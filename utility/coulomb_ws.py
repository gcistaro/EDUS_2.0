#!/usr/bin/env python3
"""
Wigner-Seitz symmetrization of KCW Coulomb files (barecoulomb*.txt, screencoulomb*.txt) for EDUS.

KCW computes V_mn(R) and W_mn(R) on an Nsc x Nsc x Nsc supercell, so they are periodic in R with period Nsc
(R and R + Nsc*L are the same point) and the file lists one representative per class, R_i in [-Nsc/2, Nsc/2 - 1].
EDUS (src/ModelCoulomb/ModelCoulomb.cpp) puts each value exactly at that R and zero elsewhere: R = -Nsc/2 is
present while +Nsc/2 is not, so the interaction is not invariant under inversion nor under the point group.
With the screened exchange this gives e.g. a transverse current along y,z for a pump polarized along x in a
cubic crystal.

As wannier90 does for H(R), every entry (R, m, n) is assigned to the image R + Nsc*L (L in {-1,0,1}^3) that
minimizes the distance between the centers, |tau_m - (R + Nsc*L).A - tau_n| (same convention as
ModelCoulomb.cpp: center m in cell 0, center n in cell R); if several images are degenerate (within --tol)
the value is split equally among them. The sum over R of every (m,n) element is unchanged and the output has
W_nm(-R) for every W_mn(R). Lattice vectors and Wannier centers (diagonal of r_mn(R=0)) are read from the
wannier90 _tb.dat file used by EDUS. Same file format in output (it can have more R vectors than the input).

Usage: python3 coulomb_ws.py model_tb.dat screencoulomb.txt screencoulomb_ws.txt [--tol 1e-4]
"""
import argparse
import itertools

import numpy as np


def read_tb(fname):
    """Lattice vectors (rows) and Wannier centers from a wannier90 _tb.dat file (both in Angstrom)."""
    lines = open(fname).read().split("\n")
    A = np.array([[float(x) for x in lines[i].split()] for i in (1, 2, 3)])
    nw, nr = int(lines[4]), int(lines[5])
    i, ndeg = 6, 0
    while ndeg < nr:
        ndeg += len(lines[i].split())
        i += 1
    for ncol in (1, 3):           # skip the H(R) blocks, then read the r(R) blocks
        for _ in range(nr):
            while not lines[i].strip():
                i += 1
            R = [int(x) for x in lines[i].split()]
            i += 1
            block = [lines[i + j].split() for j in range(nw * nw)]
            i += nw * nw
            if ncol == 3 and R == [0, 0, 0]:
                tau = np.zeros((nw, 3))
                for s in block:
                    m, n = int(s[0]) - 1, int(s[1]) - 1
                    if m == n:
                        tau[m] = [float(s[2 + 2 * c]) for c in range(3)]
                return A, tau
    raise RuntimeError("r(R=0) not found in " + fname)


def read_coulomb(fname):
    """Blocks 'R1 R2 R3' followed by lines 'm n s Re Im' (s = 1, 2 duplicated). Returns {R: {(m, n, s): (Re, Im)}}."""
    data, R = {}, None
    for line in open(fname):
        s = line.split()
        if len(s) == 3:
            R = tuple(int(x) for x in s)
            data[R] = {}
        elif len(s) == 5:
            data[R][(int(s[0]), int(s[1]), int(s[2]))] = (float(s[3]), float(s[4]))
    return data


def symmetrize(A, tau, data, tol):
    Rs = list(data)
    nsc = max(max(r) for r in Rs) - min(min(r) for r in Rs) + 1
    out = {}
    for R, block in data.items():
        imgs = [tuple(int(x) for x in np.array(R) + nsc * np.array(L)) for L in itertools.product((-1, 0, 1), repeat=3)]
        for (m, n, s), (re, im) in block.items():
            d = [np.linalg.norm(tau[m - 1] - np.array(Ri) @ A - tau[n - 1]) for Ri in imgs]
            best = [Ri for Ri, di in zip(imgs, d) if di < min(d) + tol]
            for Ri in best:
                acc = out.setdefault(Ri, {}).setdefault((m, n, s), [0., 0.])
                acc[0] += re / len(best)
                acc[1] += im / len(best)
    keys = sorted({k for b in data.values() for k in b})
    for b in out.values():          # every R block has all the (m, n, s) lines, as EDUS expects
        for k in keys:
            b.setdefault(k, [0., 0.])
    return nsc, keys, out


def write_coulomb(fname, keys, out):
    with open(fname, "w") as f:
        for R in sorted(out):
            f.write("%7d%7d%7d\n" % R)
            for k in keys:
                f.write("%12d%12d%12d  %.17g   %.17g\n" % (*k, out[R][k][0], out[R][k][1]))


def main():
    p = argparse.ArgumentParser(description="Wigner-Seitz symmetrization of KCW Coulomb files for EDUS")
    p.add_argument("tb", help="wannier90 _tb.dat file of the model (lattice and Wannier centers)")
    p.add_argument("input", help="KCW Coulomb file (barecoulomb*.txt or screencoulomb*.txt)")
    p.add_argument("output", help="symmetrized Coulomb file")
    p.add_argument("--tol", type=float, default=1e-4, help="degeneracy tolerance on the distances (Angstrom)")
    a = p.parse_args()
    A, tau = read_tb(a.tb)
    data = read_coulomb(a.input)
    nsc, keys, out = symmetrize(A, tau, data, a.tol)
    write_coulomb(a.output, keys, out)
    tot_in = sum(v[0] for b in data.values() for v in b.values())
    tot_out = sum(v[0] for b in out.values() for v in b.values())
    asym = max(abs(out[R][(m, n, s)][0] - out.get(tuple(-x for x in R), {}).get((n, m, s), [np.inf])[0])
               for R in out for (m, n, s) in out[R])
    print(f"{a.input}: supercell {nsc}^3, R {len(data)} -> {len(out)}, sum {tot_in:.10g} -> {tot_out:.10g}, "
          f"max|W_mn(R) - W_nm(-R)| = {asym:.2e}")


if __name__ == "__main__":
    main()
