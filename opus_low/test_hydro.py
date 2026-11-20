"""Quick tests (<10 s): uniform state stays uniform; mass/energy conservation of injection."""
import numpy as np, time, hydro
t0 = time.time()
for o in (1, 2):
    h = hydro.Hydro(64, 1.0, 1.0, 1.0, order=o); h.evolve_to(0.5)
    assert np.allclose(h.prim(h.U)[0], 1) and np.allclose(h.prim(h.U)[1], 0, atol=1e-10)
    h = hydro.Hydro(128, 1.0, 1.0, 1.0, order=o, Mdot=1.0, L=10.0, r_inj=0.1); h.evolve_to(0.1)
    m = (h.U[0] * h.V).sum(); m0 = h.V.sum()
    assert abs(m - m0 - 0.1) / 0.1 < 1e-6, m - m0
    e = (h.U[2] * h.V).sum() - h.V.sum() / (2/3)
    assert abs(e - 1.0) < 1e-6, e
print("tests passed in %.2f s" % (time.time() - t0))
