import numpy as np, sam
for name, nH in [("F100", 100), ("F10", 10), ("F1", 1)]:
    b = sam.Bubble(1e-4, 1e38, nH); t = 100 * sam.KYR
    R2 = b.R2(t); r = np.linspace(R2 / 2000, 1.5 * R2, 2000)
    rho, v, p, T = b.profile(r, t)
    np.savetxt(f"sam_outputs/{name}.csv", np.c_[r / sam.PC, rho, v / 1e5, p / sam.KB, T],
               delimiter=",", header="r,rho,v,P,T", comments="")
    open(f"sam_outputs/{name}.txt", "w").write(f"{R2/sam.PC:.6e} {b.P_hot(t)/sam.KB:.6e}\n")
    print(name, open(f"sam_outputs/{name}.txt").read().strip(), b.R1(t)/sam.PC)
