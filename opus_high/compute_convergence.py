"""Task 3: L1 relative-error norms of the last snapshot vs. the SAM
(evaluated with the simulation's R_inj = 0.5 pc), for r < R2_SAM."""
import os
import numpy as np

from sam.model import WindBubble
from sam.units import PC, KYR, KMS, K_B

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(HERE, "sim_outputs")
OUT = os.path.join(HERE, "convergence")
COLS = {"density": 1, "velocity": 2, "pressure": 3, "temperature": 4}


def main():
    os.makedirs(OUT, exist_ok=True)
    bubble = WindBubble(1e-4, 1e38, 100.0, R_inj_cm=0.5 * PC)
    res = {k: [] for k in COLS}
    for N in (64, 128, 256):
        t_kyr = np.loadtxt(os.path.join(SIM, "F100_sim_N%05d_times.txt" % N))[-1]
        t = t_kyr * KYR
        st = bubble.state(t)
        row = {k: [N] for k in COLS}
        for vv in ("fo", "so"):
            d = np.loadtxt(os.path.join(SIM, "F100_sim_%s_N%05d_snap005.csv" % (vv, N)),
                           delimiter=",", skiprows=1)
            r = d[:, 0]
            m = r < st["R2"] / PC
            rho, v, P, T = bubble.profile(r[m] * PC, t, st)
            sam = {"density": rho, "velocity": v / KMS, "pressure": P / K_B, "temperature": T}
            for k, c in COLS.items():
                row[k].append(np.mean(np.abs(d[m, c] - sam[k]) / np.abs(sam[k])))
        for k in COLS:
            res[k].append(row[k])
    for k in COLS:
        with open(os.path.join(OUT, "%s_L1_norm.txt" % k), "w") as f:
            for N, fo, so in res[k]:
                f.write("%d %.8e %.8e\n" % (N, fo, so))
        print(k, ["N=%d fo=%.4f so=%.4f" % tuple(x) for x in res[k]])


if __name__ == "__main__":
    main()
