"""Task 1: evaluate the semi-analytic model for the six experiments."""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from sam.model import WindBubble  # noqa: E402
from sam import units as u        # noqa: E402

EXPERIMENTS = {
    "F100": (1e-4, 1e38, 100.0, 100.0),
    "F10": (1e-4, 1e38, 10.0, 100.0),
    "F1": (1e-4, 1e38, 1.0, 100.0),
    "S100": (1e-2, 1e38, 100.0, 100.0),
    "S10": (1e-2, 1e38, 10.0, 100.0),
    "S1": (1e-2, 1e38, 1.0, 100.0),
}

OUT = os.path.join(os.path.dirname(__file__), "..", "sam_outputs")


def main(n=4000):
    os.makedirs(OUT, exist_ok=True)
    for name, (mdot, L, nH, t_kyr) in EXPERIMENTS.items():
        b = WindBubble.from_astro(mdot, L, nH, R_inj_AU=100.0)
        t = t_kyr * u.kyr
        s = b.summary(t)
        r = np.linspace(0, 1.5 * s["R2"], n + 1)[1:]
        p = b.profile(r, t)
        hdr = ("r_pc,r_cm,rho_gcm3,nH_cm3,P_dyncm2,P_over_k_Kcm3,T_K,v_kms\n"
               f"# {name}: Mdot={mdot} Msun/yr L={L} erg/s nH={nH} t={t_kyr} kyr; "
               f"R1={s['R1']/u.pc:.5g} Rc={s['Rc']/u.pc:.5g} R2={s['R2']/u.pc:.5g} pc; "
               f"t/t_ch={s['t_over_tch']:.4g}; self-similar valid={s['valid']}")
        data = np.column_stack([r / u.pc, r, p["rho"], u.rho_to_nH(p["rho"]), p["P"],
                                p["P"] / u.k_B, p["T"], p["v"] / u.km])
        with open(os.path.join(OUT, f"{name}.csv"), "w") as f:
            f.write(hdr.split("\n")[1] + "\n" + hdr.split("\n")[0] + "\n")
            np.savetxt(f, data, delimiter=",", fmt="%.8e")
        with open(os.path.join(OUT, f"{name}.txt"), "w") as f:
            f.write(f"{s['R2']/u.pc:.6e} {s['P_b']/u.k_B:.6e}\n")
        print(f"{name}: R2={s['R2']/u.pc:.4f} pc  P/k={s['P_b']/u.k_B:.4e} K/cm^3  "
              f"R1={s['R1']/u.pc:.3f} Rc={s['Rc']/u.pc:.3f}  t/tch={s['t_over_tch']:.3g} valid={s['valid']}")


if __name__ == "__main__":
    main()
