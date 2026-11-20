"""Fast unit tests for hydro.solver (run: venv/bin/python tests/test_hydro.py)."""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from hydro.solver import Hydro  # noqa: E402

# exact Sod (rhoL=1,pL=1,rhoR=0.125,pR=0.1, gamma=1.4) star values
SOD = dict(p=0.30313, u=0.92745, rhoL=0.42632, rhoR=0.26557, S=1.75216)


def _sod(order, N=400, t=0.2):
    h = Hydro(N, 1.0, gamma=1.4, order=order, geometry="cartesian",
              inner_bc="outflow", cfl=0.4)
    left = h.r < 0.5
    h.set_primitive(np.where(left, 1.0, 0.125), 0.0, np.where(left, 1.0, 0.1))
    h.evolve(t)
    return h


def test_sod():
    for order, tol in ((1, 0.03), (2, 0.01)):
        h = _sod(order)
        rho, u, p = h.W
        x = h.r
        t = 0.2
        reg_R = (x > 0.5 + SOD["u"] * t + 0.03) & (x < 0.5 + SOD["S"] * t - 0.03)
        # rarefaction tail moves at u* - c*_L = -0.070
        reg_L = (x > 0.5 - 0.070 * t + 0.03) & (x < 0.5 + SOD["u"] * t - 0.03)
        for arr, val, reg in ((rho, SOD["rhoR"], reg_R), (rho, SOD["rhoL"], reg_L),
                              (p, SOD["p"], reg_R | reg_L), (u, SOD["u"], reg_R | reg_L)):
            err = np.abs(arr[reg].mean() - val) / val
            assert err < tol, (order, val, err)
        # shock position: last cell with density above the mid-jump value
        xs = x[np.where(rho > 0.5 * (SOD["rhoR"] + 0.125))[0].max()]
        assert abs(xs - (0.5 + SOD["S"] * t)) < 0.01, (order, xs)


def test_second_order_more_accurate():
    e = []
    for order in (1, 2):
        h = _sod(order, N=200)
        rho = h.W[0]
        x = h.r
        reg = (x > 0.5 + SOD["u"] * 0.2 + 0.02) & (x < 0.5 + SOD["S"] * 0.2 - 0.02)
        e.append(np.mean(np.abs(rho[reg] - SOD["rhoR"])))
    assert e[1] < e[0], e


def test_static_spherical():
    for order in (1, 2):
        h = Hydro(64, 1.0, order=order)
        h.set_primitive(1.0, 0.0, 1.0)
        for _ in range(50):
            h.step(h.timestep())
        rho, v, p = h.W
        assert np.max(np.abs(v)) < 1e-12 and np.max(np.abs(p - 1)) < 1e-12


def test_injection_conservation():
    Mdot, L = 2.0, 50.0
    for order in (1, 2):
        h = Hydro(128, 1.0, order=order, Mdot=Mdot, L=L, R_inj=0.1)
        h.set_primitive(1.0, 0.0, 0.01)
        M0, _, E0 = h.totals()
        h.evolve(0.02)
        M, _, E = h.totals()
        assert abs(h.W[1][-1]) < 1e-12           # nothing reached boundary
        assert abs(M - (M0 + Mdot * h.t)) / M < 1e-10
        assert abs(E - (E0 + L * h.t)) / E < 1e-10


if __name__ == "__main__":
    t0 = time.time()
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
    print("all tests passed in %.2f s" % (time.time() - t0))
