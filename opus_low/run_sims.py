import numpy as np, sam, hydro, time
b = sam.Bubble(1e-4, 1e38, 100); T0 = 100.0
rmax = 1.5 * b.R2(100 * sam.KYR)  # domain from SAM outer radius
P0 = b.rho0 * sam.KB * T0 / (sam.MU * sam.MP)
for N in (64, 128, 256):
    times = [20, 40, 60, 80, 100]
    np.savetxt(f"sim_outputs/F100_sim_N{N:05d}_times.txt", times, fmt="%g")
    for o, vv in ((1, "fo"), (2, "so")):
        t0 = time.time()
        h = hydro.Hydro(N, rmax, b.rho0, P0, order=o, Mdot=b.Mdot, L=b.L, r_inj=0.5 * sam.PC)
        for k, tk in enumerate(times, 1):
            h.evolve_to(tk * sam.KYR)
            rho, v, P = h.prim(h.U); T = P * sam.MU * sam.MP / (rho * sam.KB)
            np.savetxt(f"sim_outputs/F100_sim_{vv}_N{N:05d}_snap{k:03d}.csv",
                       np.c_[h.r / sam.PC, rho, v / 1e5, P / sam.KB, T], delimiter=",",
                       header="r,rho,v,P,T", comments="")
        print(N, vv, "%.1fs" % (time.time() - t0))
