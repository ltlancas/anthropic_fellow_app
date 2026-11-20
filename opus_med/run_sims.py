"""Task 2 driver: F100 with R_inj = 0.5 pc, orders 1 & 2, N = 64, 128, 256."""
import os, time
import numpy as np
from hydro import Grid, Hydro, injection_source
from sam import WindBubble, pc, kyr, km, k_B, m_p, M_sun, yr, MU_H, MU, GAMMA

MDOT, LW, NH, T_AMB = 1e-4, 1e38, 100.0, 100.0
R_INJ = 0.5 * pc
T_SNAP = [20.0, 40.0, 60.0, 80.0, 100.0]           # kyr
R_MAX = 1.5 * WindBubble(MDOT, LW, NH, R_inj=R_INJ).R2(T_SNAP[-1] * kyr)


def run(N, order):
    g = Grid(N, R_MAX)
    h = Hydro(g, gamma=GAMMA, order=order, cfl=0.4,
              source=injection_source(R_INJ, MDOT * M_sun / yr, LW))
    rho0 = MU_H * m_p * NH
    h.set_initial(rho0, 0.0, rho0 / (MU * m_p) * k_B * T_AMB)
    vv = {1: "fo", 2: "so"}[order]
    nsteps = 0
    for k, ts in enumerate(T_SNAP, start=1):
        nsteps += h.evolve(ts * kyr)
        rho, v, P = h.primitives()
        T = P * MU * m_p / (rho * k_B)
        np.savetxt(f"sim_outputs/F100_sim_{vv}_N{N:05d}_snap{k:03d}.csv",
                   np.column_stack([g.rc / pc, rho, v / km, P / k_B, T]),
                   delimiter=",", header="r,rho,v,P,T", comments="", fmt="%.8e")
    return nsteps


if __name__ == "__main__":
    os.makedirs("sim_outputs", exist_ok=True)
    print(f"R_max = {R_MAX / pc:.4f} pc")
    for N in (64, 128, 256):
        with open(f"sim_outputs/F100_sim_N{N:05d}_times.txt", "w") as f:
            f.write("\n".join(f"{t:.6e}" for t in T_SNAP) + "\n")
        for order in (1, 2):
            t0 = time.time()
            n = run(N, order)
            print(f"N={N} order={order}: {n} steps, {time.time() - t0:.1f} s")
