#!/usr/bin/env python3
"""Joins the txt outputs of a run restarted from a checkpoint in a new directory with the ones of the
parent run: the lines of the parent before the restart time, followed by the lines of the restart.

    python3 restart_stitch.py PARENT/Output RESTART/Output JOINED_DIR

Only the standard library is used (it runs with the python3 of the cluster)."""
import os
import sys

FILES = ["Time.txt", "Laser.txt", "Laser_A.txt", "Population.txt", "Population_wannier.txt",
         "Velocity.txt", "Energy.txt"]


def data_lines(path):
    """Header (lines starting with '#') and data lines of a txt output"""
    header, data = [], []
    with open(path) as f:
        for line in f:
            (header if line.startswith("#") else data).append(line)
    return header, data


def main(parent, restart, joined):
    t_parent = [float(x) for x in data_lines(os.path.join(parent, "Time.txt"))[1]]
    t_restart = [float(x) for x in data_lines(os.path.join(restart, "Time.txt"))[1]]
    if not t_restart:
        sys.exit("no lines in the Time.txt of the restart")
    # times are printed with 6 digits: half of the print interval as tolerance
    tol = 0.5 * (t_restart[1] - t_restart[0]) if len(t_restart) > 1 else 1.e-6
    nkeep = sum(1 for t in t_parent if t < t_restart[0] - tol)
    if nkeep == len(t_parent):
        print(f"WARNING: the parent ends at t={t_parent[-1]} a.u., before the restart ({t_restart[0]} a.u.)")
    os.makedirs(joined, exist_ok=True)
    for name in FILES:
        p, r = os.path.join(parent, name), os.path.join(restart, name)
        if not (os.path.exists(p) and os.path.exists(r)):
            continue
        hp, dp = data_lines(p)
        _, dr = data_lines(r)
        with open(os.path.join(joined, name), "w") as f:
            f.writelines(hp + dp[:nkeep] + dr)
    print(f"{nkeep} lines from {parent} (t < {t_restart[0]} a.u.), {len(t_restart)} from {restart} -> {joined}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
