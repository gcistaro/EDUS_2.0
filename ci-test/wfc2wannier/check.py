#!/usr/bin/env python3
"""Checks the output of EDUS_wfc2wannier on the synthetic cases produced by generate.py.
Usage: check.py <case directory>   (exit code 1 on failure)"""
import json
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate import LAT, OMEGA, REC, gauss_coef, miller_sphere  # noqa: E402


def read_record(f):
    (n,) = struct.unpack("<i", f.read(4))
    data = f.read(n)
    (m,) = struct.unpack("<i", f.read(4))
    assert n == m
    return data


def read_wfc(fname):
    with open(fname, "rb") as f:
        ik, x, y, z, ispin, gamma_only, scalef = struct.unpack("<i3diid", read_record(f))
        ngw, igwx, npol, nbnd = struct.unpack("<4i", read_record(f))
        B = np.frombuffer(read_record(f), "<f8").reshape(3, 3)
        M = np.frombuffer(read_record(f), "<i4").reshape(igwx, 3)
        evc = np.array([np.frombuffer(read_record(f), "<c16") for _ in range(nbnd)])
    return dict(ik=ik, xk=np.array([x, y, z]), npol=npol, nbnd=nbnd, B=B, M=M, evc=evc, gamma_only=gamma_only)


def read_pair_density(fname):
    with open(fname, "rb") as f:
        nb, n1, n2, n3, nR, nk, m1, m2, m3 = struct.unpack("<9i", read_record(f))
        wl = np.frombuffer(read_record(f), "<i4")
        qa = np.frombuffer(read_record(f), "<f8")
        R = np.frombuffer(read_record(f), "<i4").reshape(nR, 3)
        rho = np.zeros((nR, nb, nb, n1, n2, n3), complex)
        for iR in range(nR):
            for a in range(nb):
                for b in range(nb):
                    rho[iR, a, b] = np.frombuffer(read_record(f), "<c16").reshape(n3, n2, n1).transpose(2, 1, 0)
    return dict(wl=wl, q=qa[:3], a=qa[3:].reshape(3, 3), R=R, nk=nk, mp=(m1, m2, m3), grid=(n1, n2, n3), rho=rho)


def u_analytic(kc, tau, sigma, chi, grid):
    """u_k(r) = 1/sqrt(Omega) sum_G c(G) e^{iG.r} of the analytic Bloch sum at any k, shape (npol, N1, N2, N3)"""
    M = miller_sphere(kc, sigma=sigma)
    buf = np.zeros(grid, complex)
    np.add.at(buf, (M[:, 0] % grid[0], M[:, 1] % grid[1], M[:, 2] % grid[2]), gauss_coef(kc, M, tau, sigma))
    u = np.fft.ifftn(buf) * np.prod(grid) / np.sqrt(OMEGA)
    return u[None] if chi is None else np.array([chi[0] * u, chi[1] * u])


def check_pair_density(d, ref, kpts, chi):
    """rho^{q,R}_{nn'}(r) = sum_k e^{ik.R} u*_{nk}(r) u_{n'k+q}(r) from the analytic u_{k+q} (no folding of k+q)"""
    pd = json.load(open(f"{d}/input.json")).get("pair_density")
    if pd is None:
        return True
    ok = True
    sigma, taus = ref["sigma"], ref["taus"]
    for iq in range(len(pd.get("q_list", [[0, 0, 0]]))):
        out = read_pair_density(f"{d}/out/pair_density_q{iq + 1}.dat")
        grid, q, wl = out["grid"], out["q"], out["wl"] - 1
        chis = [None if chi is None else chi[n] for n in wl]
        expected = np.zeros_like(out["rho"])
        for kc in kpts:
            uk = [u_analytic(kc, taus[n], sigma, c, grid) for n, c in zip(wl, chis)]
            ukq = [u_analytic(kc + q, taus[n], sigma, c, grid) for n, c in zip(wl, chis)]
            for iR, R in enumerate(out["R"]):
                ph = np.exp(2j * np.pi * kc @ R)
                for a in range(len(wl)):
                    for b in range(len(wl)):
                        expected[iR, a, b] += ph * np.sum(np.conj(uk[a]) * ukq[b], axis=0)
        err = np.max(np.abs(out["rho"] - expected)) / np.max(np.abs(expected))
        # q = 0, R = 0 (modulo the k grid): int_cell rho^{0,0}_{nn} = N_k (normalization of u_nk).
        # Only the diagonal: the synthetic Gaussian Wannier functions overlap, so the off-diagonal elements and
        # R != 0 are not delta-like here (they are covered by the comparison with the analytic reference above)
        e_int = 0.0
        if np.allclose(q, 0):
            integral = out["rho"].sum(axis=(3, 4, 5)) * OMEGA / np.prod(grid) / out["nk"]
            for iR, R in enumerate(out["R"]):
                if all(R[x] % out["mp"][x] == 0 for x in range(3)):
                    e_int = max(e_int, np.max(np.abs(np.diag(integral[iR]) - 1)))
        print(f"  pair density q = {np.round(q, 4)}: grid {grid}, max rel. error {err:.2e}"
              + (f", |int rho_nn/N_k - 1| = {e_int:.1e}" if np.allclose(q, 0) else ""))
        # u.mat has 10 decimal digits
        ok &= err < 1.0e-8 and e_int < 1.0e-8
    return ok


def read_coulomb(fname):
    """KCW format: 'R1 R2 R3' then 'n m ispin Re Im' lines -> {R: U[n, m]} (Rydberg, ispin = 1)"""
    U, R = {}, None
    for line in open(fname):
        s = line.split()
        if len(s) == 3:
            R = tuple(int(x) for x in s)
            U[R] = {}
        elif len(s) == 5 and s[2] == "1":
            U[R][(int(s[0]) - 1, int(s[1]) - 1)] = float(s[3]) + 1j * float(s[4])
    return U


def check_coulomb(d, ref, kpts, chi):
    """bare U_nm(R) = Omega/N_k^3 sum_q e^{-iq.R} sum_{G, q+G != 0, |q+G|^2 <= ecut} rho*_n v rho_m, with rho^q_nn(G)
    from the analytic Bloch sums (no folding of k+q), in Rydberg"""
    cp = json.load(open(f"{d}/input.json")).get("coulomb")
    if cp is None:
        return True
    sigma, taus = ref["sigma"], ref["taus"]
    nk, nw = len(kpts), len(taus)
    mp = [int(round(1 / min(x for x in np.unique(kpts[:, i]) if x > 0))) if np.any(kpts[:, i] > 0) else 1
          for i in range(3)]
    ecut = cp["ecut_bare"]
    # G vectors of the sphere and an FFT grid with no aliasing on them
    nmax = [int(np.sqrt(ecut) * np.linalg.norm(LAT[i]) / (2 * np.pi)) + 2 for i in range(3)]
    mu = np.max(np.abs(np.concatenate([miller_sphere(kc + q, sigma=sigma) for kc in kpts for q in [np.zeros(3), np.ones(3)]])), axis=0)
    grid = tuple(int(2 * mu[i] + 2 * nmax[i] + 3) for i in range(3))
    chis = [None if chi is None else chi[n] for n in range(nw)]
    qs = np.array([[i / mp[0], j / mp[1], k / mp[2]] for i in range(mp[0]) for j in range(mp[1]) for k in range(mp[2])])
    Vq = np.zeros((len(qs), nw, nw), complex)
    for iq, q in enumerate(qs):
        rhoG = []
        for n in range(nw):
            rho = 0
            for kc in kpts:
                uk = u_analytic(kc, taus[n], sigma, chis[n], grid)
                ukq = u_analytic(kc + q, taus[n], sigma, chis[n], grid)
                rho = rho + np.sum(np.conj(uk) * ukq, axis=0)
            rhoG.append(np.fft.fftn(rho) / np.prod(grid))
        r = [np.arange(-n, n + 1) for n in nmax]
        G = np.array(np.meshgrid(*r, indexing="ij")).reshape(3, -1).T
        Q = (q + G) @ REC
        Q2 = np.sum(Q * Q, axis=1)
        keep = (Q2 <= ecut) & (Q2 > 1e-12)
        G, Q2 = G[keep], Q2[keep]
        idx = (G[:, 0] % grid[0], G[:, 1] % grid[1], G[:, 2] % grid[2])
        rg = np.array([x[idx] for x in rhoG])
        Vq[iq] = np.conj(rg) @ (rg * 4 * np.pi / Q2).T
    U = read_coulomb(f"{d}/out/barecoulomb.txt")
    err, umax = 0.0, 0.0
    for R, Um in U.items():
        ph = np.exp(-2j * np.pi * qs @ np.array(R))
        expected = 2 * OMEGA / nk**3 * np.tensordot(ph, Vq, axes=1)
        for (n, m), v in Um.items():
            err = max(err, abs(v - expected[n, m]))
            umax = max(umax, abs(expected[n, m]))
    print(f"  bare Coulomb U_nm(R): {len(U)} R, max |U| = {umax:.4f} Ry, max error {err:.2e} Ry")
    return err < 1.0e-8 * max(umax, 1.0)


def main(d):
    ref = json.load(open(f"{d}/reference.json"))
    kpts = np.array(ref["kpts"])
    taus = ref["taus"]
    chi = None if ref["spinors"] is None else [np.array([c[0] + 1j * c[1] for c in sp]) for sp in ref["spinors"]]
    ok = True

    # 1) Bloch functions in the Wannier gauge == analytic Bloch sums of the Gaussians
    err = 0.0
    for ik, kc in enumerate(kpts):
        wf = read_wfc(f"{d}/out/wfc{ik + 1}.dat")
        assert wf["nbnd"] == len(taus) and not wf["gamma_only"]
        err = max(err, np.max(np.abs(wf["xk"] - kc @ wf["B"])))
        for n, t in enumerate(taus):
            c = gauss_coef(kc, wf["M"], t, ref["sigma"])
            if chi is not None:
                c = np.concatenate([chi[n][0] * c, chi[n][1] * c])
            err = max(err, np.max(np.abs(wf["evc"][n] - c)))
    print(f"  Bloch functions: max error {err:.2e}")
    ok &= err < 1.0e-9

    # 2) real-space Wannier functions: norm, centre, spread
    data = np.loadtxt(f"{d}/out/wannier_centres.txt", ndmin=2)
    for row in data:
        n = int(row[0]) - 1
        e_norm = abs(row[1] - 1.0)
        e_c = np.max(np.abs(row[2:5] - np.array(ref["centres_ang"][n])))
        e_s = abs(row[5] - ref["spread_ang2"])
        print(f"  WF {n + 1}: |norm-1| = {e_norm:.1e}  |centre error| = {e_c:.1e} Ang  |spread error| = {e_s:.1e} Ang^2")
        ok &= e_norm < 1.0e-6 and e_c < 1.0e-6 and e_s < 1.0e-5

    # 3) xsf files exist and are complete
    for n in range(len(taus)):
        txt = open(f"{d}/out/wannier_{n + 1:05d}.xsf").read()
        ok &= "END_BLOCK_DATAGRID_3D" in txt

    # 4) pair densities
    ok &= check_pair_density(d, ref, kpts, chi)

    # 5) bare Coulomb matrix elements
    ok &= check_coulomb(d, ref, kpts, chi)
    return ok


if __name__ == "__main__":
    good = main(sys.argv[1])
    print("  ->", "PASSED" if good else "FAILED")
    sys.exit(0 if good else 1)
