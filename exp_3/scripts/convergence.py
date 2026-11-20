"""Task 3: L2 comparison of the F100 simulations against the semi-analytic model.

For each snapshot and variable q, the relative L2 norm on the (uniform) grid points is
    L2 = sqrt( sum_i (q_sim,i - q_sam,i)^2 ) / sqrt( sum_i q_sam,i^2 ),
and the reported value is the mean of L2 over the 5 snapshots.
"""
import glob
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from sam.model import WindBubble  # noqa: E402
from sam import units as u        # noqa: E402

HERE = os.path.join(os.path.dirname(__file__), "..")
SIM = os.path.join(HERE, "sim_outputs")
OUT = os.path.join(HERE, "convergence")
VARS = {"density": 2, "velocity": 4, "pressure": 5, "temperature": 6}   # sim csv columns
SAMKEY = {"density": "rho", "velocity": "v", "pressure": "P", "temperature": "T"}


def l2(a, b):
    return np.sqrt(np.sum((a - b) ** 2)) / np.sqrt(np.sum(b ** 2))


def main(r_excl_pc=0.0, suffix=""):
    os.makedirs(OUT, exist_ok=True)
    sam = WindBubble.from_astro(1e-4, 1e38, 100.0, R_inj_AU=100.0, T_amb=10.0)
    Ns = sorted({int(f.split("_N")[1][:5]) for f in glob.glob(os.path.join(SIM, "F100_sim_*_snap*.csv"))})
    res = {v: np.zeros((len(Ns), 2)) for v in VARS}
    detail = []
    for j, N in enumerate(Ns):
        times = np.loadtxt(os.path.join(SIM, f"F100_sim_N{N:05d}_times.txt"), usecols=1)
        cache = {}
        for c, tag in enumerate(("fo", "so")):
            acc = {v: [] for v in VARS}
            for k, t in enumerate(times, start=1):
                d = np.loadtxt(os.path.join(SIM, f"F100_sim_{tag}_N{N:05d}_snap{k:03d}.csv"),
                               delimiter=",")
                if k not in cache:
                    cache[k] = sam.profile(d[:, 0], t)
                p = cache[k]
                m = d[:, 1] >= r_excl_pc
                for v, col in VARS.items():
                    e = l2(d[m, col], p[SAMKEY[v]][m])
                    acc[v].append(e)
                    detail.append((N, tag, k, v, e))
            for v in VARS:
                res[v][j, c] = np.mean(acc[v])
    for v in VARS:
        with open(os.path.join(OUT, f"{v}_L2_norm{suffix}.txt"), "w") as f:
            f.write(f"# r >= {r_excl_pc} pc only\n" if r_excl_pc else "")
            f.write("# N  mean_relative_L2_first_order  mean_relative_L2_second_order\n")
            for j, N in enumerate(Ns):
                f.write(f"{N:6d} {res[v][j,0]:.6e} {res[v][j,1]:.6e}\n")
        print(v, "\n", np.column_stack([Ns, res[v]]))
    with open(os.path.join(OUT, f"per_snapshot_L2{suffix}.txt"), "w") as f:
        f.write("# N order snap variable relative_L2\n")
        for row in detail:
            f.write("%d %s %03d %s %.6e\n" % row)


if __name__ == "__main__":
    main()
    # supplementary: exclude the unresolved injection/free-wind core
    main(r_excl_pc=0.1, suffix="_excl_r_lt_0.1pc")
