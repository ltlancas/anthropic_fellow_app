"""Task 2: run the F100 wind bubble with the hydro code.

Usage: .venv/bin/python scripts/run_F100.py --N 1024 --order 2
"""
import argparse
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from hydro.solver import Grid, Simulation  # noqa: E402
from sam.model import WindBubble           # noqa: E402
from sam import units as u                 # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "sim_outputs")
T_END_KYR = 100.0
N_SNAP = 5
T_AMB = 10.0


def setup(N, order, cfl=0.4):
    mdot, L, nH = u.msun_per_yr_to_cgs(1e-4), 1e38, 100.0
    sam = WindBubble(mdot, L, nH, R_inj=100 * u.AU, T_amb=T_AMB)
    r_max = 1.5 * sam.R2(T_END_KYR * u.kyr)          # domain sized from the SAM
    grid = Grid(N, r_max)
    rho0 = u.nH_to_rho(nH)
    W0 = np.vstack([np.full(N, rho0), np.zeros(N), np.full(N, u.pressure_from_T(rho0, T_AMB))])
    sim = Simulation(grid, W0, order=order, cfl=cfl, mdot=mdot, L=L, R_inj=100 * u.AU)
    return sim


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--order", type=int, choices=(1, 2), required=True)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    sim = setup(a.N, a.order)
    tag = {1: "fo", 2: "so"}[a.order]
    times = np.linspace(T_END_KYR / N_SNAP, T_END_KYR, N_SNAP) * u.kyr
    t0 = time.time()
    for k, t in enumerate(times, start=1):
        sim.run_to(t)
        W = sim.primitives()
        r = sim.grid.rc
        data = np.column_stack([r, r / u.pc, W[0], u.rho_to_nH(W[0]), W[1], W[2],
                                u.temperature(W[0], W[2])])
        fn = os.path.join(OUT, f"F100_sim_{tag}_N{a.N:05d}_snap{k:03d}.csv")
        np.savetxt(fn, data, delimiter=",", fmt="%.8e",
                   header=f"t_s={sim.t:.8e} t_kyr={sim.t/u.kyr:.6f} order={a.order} N={a.N}\n"
                          "r_cm,r_pc,rho_gcm3,nH_cm3,v_cms,P_dyncm2,T_K", comments="# ")
        print(f"[{tag} N={a.N}] snap {k} t={sim.t/u.kyr:.2f} kyr steps={sim.nsteps} "
              f"wall={time.time()-t0:.1f}s", flush=True)
    with open(os.path.join(OUT, f"F100_sim_N{a.N:05d}_times.txt"), "w") as f:
        f.write("# snapshot  t_s  t_kyr\n")
        for k, t in enumerate(times, start=1):
            f.write(f"{k:03d} {t:.8e} {t/u.kyr:.6f}\n")


if __name__ == "__main__":
    main()
