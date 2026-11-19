"""
Fiducial run: L_w = 1e38 erg/s, Mdot = 1e-2 Msun/yr (default), n_H = 100 cm^-3,
mu_H = 1.4, t = 1e5 yr.  Writes density and pressure (plus T and v) versus
radius out to 1.5 R2 into output/.

    .venv/bin/python scripts/run_fiducial.py                # Mdot = 1e-2 Msun/yr
    .venv/bin/python scripts/run_fiducial.py --mdot 1e-4    # Mdot = 1e-4 Msun/yr

Outputs for non-default Mdot carry a suffix, e.g. *_t1e5yr_mdot1e-04.csv.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from wind_bubble import WindBubble, units as u

OUT = Path(__file__).resolve().parents[1] / "output"
OUT.mkdir(exist_ok=True)

# ---- inputs in astrophysical units -> CGS ---------------------------------
L_w = 1e38 * u.erg_s
ap = argparse.ArgumentParser()
ap.add_argument("--mdot", type=float, default=1e-2, help="mass-loss rate [Msun/yr]")
args = ap.parse_args()
Mdot = args.mdot * u.Msun / u.yr
tag = "_t1e5yr" + ("" if args.mdot == 1e-2 else f"_mdot{args.mdot:.0e}")
n_H = 100.0 * u.per_cm3
mu_H = 1.4
t = 1e5 * u.yr

wb = WindBubble(L_w, Mdot, n_H, mu_H=mu_H).solve(t)
prof = wb.profiles(t, n_r=3000, r_max_factor=1.5)
st, diag = prof["state"], prof["diagnostics"]

# ---- table -----------------------------------------------------------------
cols = np.column_stack([
    prof["r"] / u.pc, prof["r"], prof["rho"], prof["n_H"], prof["P"],
    prof["P"] / u.k_B, prof["T"], prof["v"] / u.km_s, prof["region"],
])
header = ("r_pc,r_cm,rho_g_cm3,n_H_cm3,P_dyn_cm2,P_over_k_K_cm3,T_K,v_km_s,"
          "region(0=free_wind,1=shocked_wind,2=shell,3=ambient)")
np.savetxt(OUT / f"profiles{tag}.csv", cols, delimiter=",", header=header,
           fmt=["%.6e"] * 8 + ["%d"], comments="")

summary = {
    "inputs": {"L_w_erg_s": L_w, "Mdot_Msun_yr": Mdot / (u.Msun / u.yr),
               "n_H_cm3": n_H, "mu_H": mu_H, "t_yr": t / u.yr},
    "derived": {"v_w_km_s": wb.v_w / u.km_s, "rho0_g_cm3": wb.rho0,
                "t_swept_yr": wb.t_sw / u.yr, "P0_dyn_cm2": wb.P0},
    "state": {"R1_pc": st.R1 / u.pc, "Rc_pc": st.Rc / u.pc, "R2_pc": st.R2 / u.pc,
              "shell_thickness_pc": (st.R2 - st.Rc) / u.pc,
              "V_shell_km_s": st.V / u.km_s, "dR1dt_km_s": st.dR1dt / u.km_s,
              "P_bubble_dyn_cm2": st.P_b, "P_bubble_over_k_K_cm3": st.P_b / u.k_B,
              "E_bubble_erg": st.E_b, "M_shocked_wind_Msun": st.M_b / u.Msun,
              "M_shell_Msun": st.M_s / u.Msun},
    "weaver1977_ratio": wb.weaver_comparison(t),
    "diagnostics": diag,
}
summary = json.loads(json.dumps(summary, default=float))
(OUT / f"summary{tag}.json").write_text(json.dumps(summary, indent=2))

for k, v in summary["state"].items():
    print(f"{k:28s} {v:.4g}")
print("ratio to Weaver+77:", {k: round(v, 3) for k, v in summary["weaver1977_ratio"].items()})
print("L_cool/L_w (diagnostic):", round(diag["L_cool_over_L_w"], 2))

# ---- figure ----------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit("matplotlib not installed; table written, figure skipped")

ink, muted, line = "#1f2328", "#6e7781", "#2a6fdb"
r = prof["r"] / u.pc
panels = [("n_H", r"$n_{\rm H}$ [cm$^{-3}$]", prof["n_H"], True),
          ("P", r"$P$ [dyn cm$^{-2}$]", prof["P"], True),
          ("T", r"$T$ [K]", prof["T"], True),
          ("v", r"$v$ [km s$^{-1}$]", prof["v"] / u.km_s, False)]
fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, constrained_layout=True)
for ax, (_, lab, y, logy) in zip(axes.flat, panels):
    ax.plot(r, y, color=line, lw=2)
    if logy:
        ax.set_yscale("log")
    for R, name in ((st.R1, "$R_1$"), (st.Rc, "$R_c$"), (st.R2, "$R_2$")):
        ax.axvline(R / u.pc, color=muted, lw=0.8, ls="--")
    ax.set_ylabel(lab, color=ink)
    ax.grid(alpha=0.25, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
ax0 = axes[0, 0]
ymax = ax0.get_ylim()[1]
for R, name, ha in ((st.R1, "$R_1$", "center"), (st.Rc, "$R_c$ ", "right"),
                    (st.R2, " $R_2$", "left")):
    ax0.text(R / u.pc, ymax, name, ha=ha, va="bottom", color=muted, fontsize=9)
for ax in axes[1]:
    ax.set_xlabel("r [pc]", color=ink)
    ax.set_xlim(0, r.max())
fig.suptitle(r"Wind bubble at $t=10^5$ yr  ($L_w=10^{38}$ erg/s, "
             rf"$\dot M=10^{{{np.log10(args.mdot):.0f}}}\,M_\odot$/yr, $n_{{\rm H}}=100$ cm$^{{-3}}$)",
             color=ink, fontsize=11)
fig.savefig(OUT / f"profiles{tag}.png", dpi=150)
print("wrote", OUT)
