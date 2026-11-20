"""Solve the semi-analytic model for the six bubbles and write outputs.

For each experiment <name> writes (in the output directory, default '.'):
  <name>.csv : radius_pc, density_g_cm3, pressure_dyn_cm2, temperature_K,
               velocity_km_s   for 0 <= r <= 1.5 R2 at t_out
  <name>.txt : R_bubble [pc] (forward-shock radius R2) and P_hot/k_B
               [K cm^-3] (pressure of the shocked-wind gas at the contact)
"""
import argparse
import os

import numpy as np

from windbubble import WindBubble, units as U

EXPERIMENTS = [
    # name,  Mdot [Msun/yr], L_w [erg/s], n_H [cm^-3], t_out [kyr]
    ("F100", 1e-4, 1e38, 100.0, 100.0),
    ("F10", 1e-4, 1e38, 10.0, 100.0),
    ("F1", 1e-4, 1e38, 1.0, 100.0),
    ("S100", 1e-2, 1e38, 100.0, 100.0),
    ("S10", 1e-2, 1e38, 10.0, 100.0),
    ("S1", 1e-2, 1e38, 1.0, 100.0),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    for name, mdot, Lw, nH, tk in EXPERIMENTS:
        t = tk * U.KYR
        b = WindBubble(mdot, Lw, nH, r_inj_au=100.0).solve(t)
        r = b.output_grid(t, rmax_fac=1.5)
        p = b.profile(r, t)
        s = p["state"]
        hdr = ("# %s: Mdot=%g Msun/yr, L_w=%g erg/s, n_H=%g cm^-3, t=%g kyr, "
               "r_inj=100 AU, mu_H=1.4, mu=%.4f\n" % (name, mdot, Lw, nH, tk, U.MU) +
               "radius_pc,density_g_cm3,pressure_dyn_cm2,temperature_K,velocity_km_s")
        np.savetxt(os.path.join(a.outdir, name + ".csv"),
                   np.column_stack([r / U.PC, p["rho"], p["P"], p["T"], p["u"] / U.KM]),
                   delimiter=",", header=hdr, comments="", fmt="%.10e")
        with open(os.path.join(a.outdir, name + ".txt"), "w") as f:
            f.write("# R_bubble_pc  P_hot_over_k_K_cm3\n")
            f.write("%.6e %.6e\n" % (s["R2"] / U.PC, s["Pb"] / U.K_B))
        rows.append((name, s["R1"] / U.PC, s["Rc"] / U.PC, s["R2"] / U.PC,
                     s["Pb"] / U.K_B, s["T_hot"], s["q"], s["E_tot"] / s["E_inj"],
                     p["state"].get("vol_rescale", np.nan)))
        print("%-5s R1=%.3f Rc=%.3f R2=%.3f pc  P/k=%.4e  T_hot=%.3e  q=%.3f "
              "E/Lt=%.4f volfac=%.3f" % rows[-1], flush=True)


if __name__ == "__main__":
    main()
