"""Task 1: evaluate the semi-analytic model for F100, F10, F1."""
import os
import numpy as np

from sam.model import WindBubble
from sam.units import PC, KYR, KMS, K_B

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sam_outputs")

EXPERIMENTS = {  # name: (Mdot [Msun/yr], L [erg/s], nH [cm^-3], t_out [kyr])
    "F100": (1e-4, 1e38, 100.0, 100.0),
    "F10": (1e-4, 1e38, 10.0, 100.0),
    "F1": (1e-4, 1e38, 1.0, 100.0),
}


def main():
    os.makedirs(OUT, exist_ok=True)
    for name, (Mdot, L, nH, tk) in EXPERIMENTS.items():
        b = WindBubble(Mdot, L, nH)
        t = tk * KYR
        st = b.state(t)
        rmax = 1.5 * st["R2"]
        r = np.unique(np.concatenate([
            np.geomspace(1e-3 * b.R_inj, rmax, 600),
            np.linspace(0, rmax, 3001)[1:],
        ]))
        rho, v, P, T = b.profile(r, t, st)
        data = np.column_stack([r / PC, rho, v / KMS, P / K_B, T])
        np.savetxt(os.path.join(OUT, name + ".csv"), data, delimiter=",",
                   header="r,rho,v,P,T", comments="", fmt="%.8e")
        with open(os.path.join(OUT, name + ".txt"), "w") as f:
            f.write("%.8e %.8e\n" % (st["R2"] / PC, st["Pb"] / K_B))
        print("%s: xi2=%.4f R1=%.3f pc Rc=%.3f pc R2=%.3f pc  P_b/k=%.4e K/cm3"
              % (name, st["xi2"], st["R1"] / PC, st["Rc"] / PC, st["R2"] / PC,
                 st["Pb"] / K_B))


if __name__ == "__main__":
    main()
