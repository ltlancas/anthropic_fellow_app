import numpy as np, sam
b = sam.Bubble(1e-4, 1e38, 100); t = 100 * sam.KYR; R2 = b.R2(t) / sam.PC
cols = {"density": 1, "velocity": 2, "pressure": 3, "temperature": 4}
res = {k: [] for k in cols}
for N in (64, 128, 256):
    row = {k: [N] for k in cols}
    for vv in ("fo", "so"):
        d = np.loadtxt(f"sim_outputs/F100_sim_{vv}_N{N:05d}_snap005.csv", delimiter=",", skiprows=1)
        m = d[:, 0] < R2; r = d[m, 0] * sam.PC
        rho, v, p, T = b.profile(r, t)
        ref = {"density": rho, "velocity": v / 1e5, "pressure": p / sam.KB, "temperature": T}
        for k, c in cols.items():
            row[k].append(np.mean(np.abs(d[m, c] - ref[k]) / np.abs(ref[k])))
    for k in cols: res[k].append(row[k])
for k in cols:
    np.savetxt(f"convergence/{k}_L1_norm.txt", res[k], fmt=["%d", "%.6e", "%.6e"])
    print(k, res[k])
