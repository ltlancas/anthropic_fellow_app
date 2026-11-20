"""Task 1 driver: produce sam_outputs/ for F100, F10, F1."""
import os
import numpy as np
from sam import WindBubble, AU

EXPS = {"F100": 100.0, "F10": 10.0, "F1": 1.0}
MDOT, LW, TOUT = 1e-4, 1e38, 100.0

os.makedirs("sam_outputs", exist_ok=True)
for name, nH in EXPS.items():
    b = WindBubble(MDOT, LW, nH, R_inj=100 * AU)
    s = b.summary(TOUT)
    R2 = s["R2_pc"]
    r = np.unique(np.concatenate([
        np.logspace(-6, np.log10(R2 * 1.5), 600),
        np.linspace(R2 * 1e-3, R2 * 1.5, 1500)]))
    rr, rho, v, Pk, T = b.profiles_physical_units(r, TOUT)
    np.savetxt(f"sam_outputs/{name}.csv", np.column_stack([rr, rho, v, Pk, T]),
               delimiter=",", header="r,rho,v,P,T", comments="", fmt="%.8e")
    with open(f"sam_outputs/{name}.txt", "w") as f:
        f.write(f"{R2:.8e} {s['Pb_over_k']:.8e}\n")
    print(name, s)
