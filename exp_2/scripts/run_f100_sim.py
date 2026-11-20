"""Run the 1D hydro simulation of a wind bubble (default: F100).

Usage: python scripts/run_f100_sim.py --order {1,2} --n N [--outdir DIR]
                                      [--problem F100|F10|F1|S100|S10|S1]

The domain is set from the semi-analytic model: r_out = 1.25 R2(t_end), and
the inner boundary r_in = 0.1 pc lies inside the free-wind region
(R1 > 0.1 pc for t >~ 0.1 kyr).  Ghost cells at r < r_in are filled with the
exact steady Chevalier & Clegg (1985) wind produced by uniform injection of
Mdot and L_w inside 100 AU, which is equivalent to volumetric injection
because the wind is steady and supersonic (all characteristics enter the
domain) at r_in.
"""
import argparse
import os
import time

import numpy as np

from hydro1d import Hydro1D, CodeUnits
from windbubble import WindBubble, units as U
from windbubble.model import CC85Wind

# --- problem parameters (Mdot [Msun/yr], L_w [erg/s], n_H [cm^-3]) --------
PROBLEMS = {"F100": (1e-4, 1e38, 100.0), "F10": (1e-4, 1e38, 10.0),
            "F1": (1e-4, 1e38, 1.0), "S100": (1e-2, 1e38, 100.0),
            "S10": (1e-2, 1e38, 10.0), "S1": (1e-2, 1e38, 1.0)}
MDOT, LW, NH = PROBLEMS["F100"]
T_END_KYR, R_INJ_AU, T0 = 100.0, 100.0, 100.0
N_SNAP = 5
R_IN_PC = 0.1
R_OUT_FAC = 1.25


def snapshot_times():
    return T_END_KYR * np.arange(1, N_SNAP + 1) / N_SNAP     # kyr


def domain():
    b = WindBubble(MDOT, LW, NH, r_inj_au=R_INJ_AU, T0=T0).solve(T_END_KYR * U.KYR)
    R2 = b.state(T_END_KYR * U.KYR)["R2"] / U.PC
    return R_IN_PC, R_OUT_FAC * R2


def write_snapshot(fname, h, cu):
    rho, u, P = h.primitive()
    rho_c = cu.to_cgs(rho, "density")
    P_c = cu.to_cgs(P, "pressure")
    u_c = cu.to_cgs(u, "velocity")
    T = U.temperature(P_c, rho_c)
    hdr = ("# t = %.10e kyr\n" % h.t +
           "radius_pc,density_g_cm3,velocity_km_s,pressure_dyn_cm2,temperature_K")
    np.savetxt(fname, np.column_stack([h.r, rho_c, u_c / U.KM, P_c, T]),
               delimiter=",", header=hdr, comments="", fmt="%.10e")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--order", type=int, choices=(1, 2), required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--outdir", default="output/sim")
    ap.add_argument("--cfl", type=float, default=0.4)
    ap.add_argument("--problem", default="F100", choices=sorted(PROBLEMS))
    a = ap.parse_args()
    global MDOT, LW, NH
    MDOT, LW, NH = PROBLEMS[a.problem]
    os.makedirs(a.outdir, exist_ok=True)

    cu = CodeUnits()
    r_in, r_out = domain()
    wind = CC85Wind(U.mdot_cgs(MDOT), LW, R_INJ_AU * U.AU)
    rho0 = U.rho_from_nH(NH)
    P0 = U.pressure_from_T(T0, rho0)

    def inner_bc(r, t):
        rho, u, P = wind.state(r * cu.length)
        return (cu.to_code(rho, "density"), cu.to_code(u, "velocity"),
                cu.to_code(P, "pressure"))

    h = Hydro1D(r_in, r_out, a.n, order=a.order, cfl=a.cfl, limiter="mc",
                inner_bc=inner_bc, outer_bc="outflow",
                rho_floor=1e-12, P_floor=1e-22)
    h.set_primitive(cu.to_code(rho0, "density"), 0.0, cu.to_code(P0, "pressure"))
    # the wind front needs ~ r_in / v_w to reach the inner boundary
    h.t = r_in * cu.length / wind.vinf / cu.time

    tag = "fo" if a.order == 1 else "so"
    times = snapshot_times()
    t0 = time.time()
    for k, ts in enumerate(times, start=1):
        h.evolve(ts)
        fn = os.path.join(a.outdir, "%s_sim_%s_N%05d_snap%03d.csv" % (a.problem, tag, a.n, k))
        write_snapshot(fn, h, cu)
        print("%s N=%d snap %d t=%.3f kyr steps=%d wall=%.1fs" %
              (tag, a.n, k, h.t, h.nstep, time.time() - t0), flush=True)
    with open(os.path.join(a.outdir, "%s_sim_N%05d_times.txt" % (a.problem, a.n)), "w") as f:
        f.write("# snapshot  time_kyr\n")
        for k, ts in enumerate(times, start=1):
            f.write("%03d %.6f\n" % (k, ts))


if __name__ == "__main__":
    main()
