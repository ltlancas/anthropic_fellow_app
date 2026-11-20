"""Fast (<10 s) tests for the hydro code.  Run: venv/bin/python tests/test_hydro.py"""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from hydro import Grid, Hydro, injection_source

G = 1.4


def sod(order, N=400):
    g = Grid(N, 1.0, geometry="planar")
    h = Hydro(g, gamma=G, order=order, inner_bc="outflow", outer_bc="outflow", cfl=0.5)
    left = g.rc < 0.5
    h.set_initial(np.where(left, 1.0, 0.125), 0.0, np.where(left, 1.0, 0.1))
    h.evolve(0.2)
    return g, h.primitives()


def test_sod():
    # exact star state: p* = 0.30313, u* = 0.92745, rho*L = 0.42632, rho*R = 0.26557
    for order in (1, 2):
        g, (rho, v, P) = sod(order)
        plateau = (g.rc > 0.60) & (g.rc < 0.66)   # between contact (0.685) and rarefaction tail
        assert abs(P[plateau].mean() - 0.30313) < 3e-3, P[plateau].mean()
        assert abs(v[plateau].mean() - 0.92745) < 1e-2
        assert abs(rho[plateau].mean() - 0.42632) < 1e-2
        sh = (g.rc > 0.72) & (g.rc < 0.82)
        assert abs(rho[sh].mean() - 0.26557) < 1e-2
    # second order should be more accurate than first order in the L1 sense
    def exact_rho(x):
        # analytic Sod density at t=0.2
        t = 0.2; cl = np.sqrt(G); x0 = 0.5
        ps, us = 0.30313, 0.92745
        rl = 1.0 * (ps / 1.0)**(1 / G); csl = np.sqrt(G * ps / rl)
        sh = 1.75216
        r = np.empty_like(x)
        xi = (x - x0) / t
        r[xi < -cl] = 1.0
        fan = (xi >= -cl) & (xi < us - csl)
        u = 2 / (G + 1) * (cl + xi[fan]); c = cl - (G - 1) / 2 * u
        r[fan] = (c / cl)**(2 / (G - 1))
        r[(xi >= us - csl) & (xi < us)] = rl
        r[(xi >= us) & (xi < sh)] = 0.26557
        r[xi >= sh] = 0.125
        return r
    errs = []
    for order in (1, 2):
        g, (rho, v, P) = sod(order, 200)
        errs.append(np.mean(np.abs(rho - exact_rho(g.rc))))
    assert errs[1] < 0.7 * errs[0], errs


def test_uniform_spherical():
    for order in (1, 2):
        g = Grid(64, 1.0)
        h = Hydro(g, order=order)
        h.set_initial(1.0, 0.0, 1.0)
        h.evolve(0.5)
        rho, v, P = h.primitives()
        assert np.max(np.abs(rho - 1)) < 1e-12 and np.max(np.abs(v)) < 1e-12
        assert np.max(np.abs(P - 1)) < 1e-12


def test_conservation_spherical():
    for order in (1, 2):
        g = Grid(100, 1.0)
        h = Hydro(g, order=order, outer_bc="reflect")
        h.set_initial(1.0, 0.0, np.where(g.rc < 0.2, 100.0, 1.0))
        M0 = np.sum(h.U[0] * g.vol); E0 = np.sum(h.U[2] * g.vol)
        h.evolve(0.1)
        M1 = np.sum(h.U[0] * g.vol); E1 = np.sum(h.U[2] * g.vol)
        assert abs(M1 / M0 - 1) < 1e-12 and abs(E1 / E0 - 1) < 1e-12


def test_injection_budget():
    g = Grid(80, 1.0)
    Mdot, L, Rinj = 1.0, 1.0, 0.13   # Rinj not aligned with faces
    h = Hydro(g, order=2, outer_bc="reflect", source=injection_source(Rinj, Mdot, L))
    h.set_initial(1.0, 0.0, 1.0)
    M0 = 4 * np.pi * np.sum(h.U[0] * g.vol); E0 = 4 * np.pi * np.sum(h.U[2] * g.vol)
    h.evolve(0.05)
    M1 = 4 * np.pi * np.sum(h.U[0] * g.vol); E1 = 4 * np.pi * np.sum(h.U[2] * g.vol)
    assert abs((M1 - M0) / (Mdot * 0.05) - 1) < 1e-9
    assert abs((E1 - E0) / (L * 0.05) - 1) < 1e-9


def test_wind_smoke():
    """Short physical-units wind bubble: stays finite and positive."""
    pc, kyr, mp, kB = 3.0857e18, 3.15576e10, 1.6726e-24, 1.380649e-16
    g = Grid(48, 4 * pc)
    rho0 = 1.4 * mp * 100
    h = Hydro(g, order=2, source=injection_source(0.5 * pc, 1e-4 * 1.989e33 / 3.15576e7, 1e38))
    h.set_initial(rho0, 0.0, 100 * kB * 100)
    h.evolve(5 * kyr)
    rho, v, P = h.primitives()
    assert np.all(np.isfinite(rho)) and np.all(rho > 0) and np.all(P > 0)
    assert np.max(v) > 1e7   # wind has been launched



if __name__ == "__main__":
    t0 = time.time()
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f(); print("PASS", name)
    print(f"all tests passed in {time.time() - t0:.2f} s")
