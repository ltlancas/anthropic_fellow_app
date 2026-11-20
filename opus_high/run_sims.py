"""Task 2: run the F100 wind bubble with R_inj = 0.5 pc at first and second
order, N = 64, 128, 256.  Code units: pc, kyr, Msun."""
import os
import time
import numpy as np

from hydro.solver import Hydro
from sam.units import CODE, MSUN, YR, K_B, KMS, nH_to_rho, pressure_from_T, temperature

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sim_outputs")

MDOT = 1e-4 * MSUN / YR       # g/s
LW = 1e38                     # erg/s
NH = 100.0
T_AMB = 100.0                 # K (same as SAM)
R_INJ = 0.5                   # pc
R_MAX = 10.0                  # pc; SAM gives R2(100 kyr) = 7.62 pc -> 1.31 R2
T_OUT = [20.0, 40.0, 60.0, 80.0, 100.0]   # kyr


def write_snapshot(h, k, order, N):
    rho, v, P = h.W
    rho_c = rho * CODE.density
    P_c = P * CODE.pressure
    data = np.column_stack([h.r, rho_c, v * CODE.velocity / KMS, P_c / K_B,
                            temperature(P_c, rho_c)])
    name = "F100_sim_%s_N%05d_snap%03d.csv" % ("fo" if order == 1 else "so", N, k + 1)
    np.savetxt(os.path.join(OUT, name), data, delimiter=",", header="r,rho,v,P,T",
               comments="", fmt="%.8e")


def run(order, N):
    rho0 = nH_to_rho(NH)
    P0 = pressure_from_T(rho0, T_AMB)
    h = Hydro(N, R_MAX, order=order, geometry="spherical", cfl=0.4,
              Mdot=CODE.to_code(MDOT, "mdot"), L=CODE.to_code(LW, "luminosity"),
              R_inj=R_INJ, rho_floor=1e-10 * CODE.to_code(rho0, "density"),
              P_floor=1e-6 * CODE.to_code(P0, "pressure"))
    h.set_primitive(CODE.to_code(rho0, "density"), 0.0, CODE.to_code(P0, "pressure"))
    t0 = time.time()
    h.evolve(T_OUT[-1], T_OUT, callback=lambda hh, k: write_snapshot(hh, k, order, N))
    print("order=%d N=%d steps=%d wall=%.1fs" % (order, N, h.nsteps, time.time() - t0))


def main():
    os.makedirs(OUT, exist_ok=True)
    for N in (64, 128, 256):
        with open(os.path.join(OUT, "F100_sim_N%05d_times.txt" % N), "w") as f:
            f.write("\n".join("%.6f" % t for t in T_OUT) + "\n")
        for order in (1, 2):
            run(order, N)


if __name__ == "__main__":
    main()
