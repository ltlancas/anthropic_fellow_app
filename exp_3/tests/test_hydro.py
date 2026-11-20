"""Fast (<10 s) verification tests for the hydro code. Run: .venv/bin/python -m pytest tests -q"""
import os
import sys
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from hydro.solver import Grid, Simulation  # noqa: E402

g = 5.0 / 3.0


def exact_riemann(x, t, WL, WR, x0=0.5, gam=1.4):
    rL, uL, pL = WL; rR, uR, pR = WR
    cL, cR = np.sqrt(gam * pL / rL), np.sqrt(gam * pR / rR)

    def fk(p, rk, pk, ck):
        if p > pk:
            A, B = 2 / ((gam + 1) * rk), (gam - 1) / (gam + 1) * pk
            return (p - pk) * np.sqrt(A / (p + B))
        return 2 * ck / (gam - 1) * ((p / pk) ** ((gam - 1) / (2 * gam)) - 1)

    ps = brentq(lambda p: fk(p, rL, pL, cL) + fk(p, rR, pR, cR) + uR - uL, 1e-10, 1e3)
    us = 0.5 * (uL + uR) + 0.5 * (fk(ps, rR, pR, cR) - fk(ps, rL, pL, cL))
    out = np.zeros((3, x.size))
    for i, xi in enumerate(x):
        s = (xi - x0) / t
        if s < us:  # left side; Sod: rarefaction on left
            rs = rL * (ps / pL) ** (1 / gam)
            cs = cL * (ps / pL) ** ((gam - 1) / (2 * gam))
            if s < uL - cL:
                out[:, i] = WL
            elif s > us - cs:
                out[:, i] = rs, us, ps
            else:
                uu = 2 / (gam + 1) * (cL + (gam - 1) / 2 * uL + s)
                cc = uu - s
                rr = rL * (cc / cL) ** (2 / (gam - 1))
                out[:, i] = rr, uu, pL * (rr / rL) ** gam
        else:  # right side; Sod: shock
            rs = rR * ((ps / pR + (gam - 1) / (gam + 1)) / ((gam - 1) / (gam + 1) * ps / pR + 1))
            S = uR + cR * np.sqrt((gam + 1) / (2 * gam) * ps / pR + (gam - 1) / (2 * gam))
            out[:, i] = (rs, us, ps) if s < S else WR
    return out


def _sod(N, order):
    grid = Grid(N, 1.0, geometry="cartesian")
    x = grid.rc
    W0 = np.where(x < 0.5, [[1.0], [0.0], [1.0]], [[0.125], [0.0], [0.1]])
    sim = Simulation(grid, W0, gamma=1.4, order=order, cfl=0.4)
    sim.run_to(0.2)
    ex = exact_riemann(x, 0.2, (1.0, 0.0, 1.0), (0.125, 0.0, 0.1))
    return np.mean(np.abs(sim.primitives()[0] - ex[0]))


def test_sod_converges_and_second_order_better():
    e1 = [_sod(N, 1) for N in (100, 200, 400)]
    e2 = [_sod(N, 2) for N in (100, 200, 400)]
    assert e1[2] < e1[1] < e1[0]
    assert e2[2] < e2[1] < e2[0]
    assert e2[1] < e1[1]
    assert e2[2] < 0.01


def test_static_uniform_sphere_is_preserved():
    grid = Grid(64, 1.0)
    W0 = np.vstack([np.ones(64), np.zeros(64), np.ones(64)])
    for order in (1, 2):
        sim = Simulation(grid, W0, gamma=g, order=order)
        sim.run_to(0.5)
        W = sim.primitives()
        assert np.max(np.abs(W[1])) < 1e-12
        assert np.max(np.abs(W[0] - 1)) < 1e-12


def test_sedov_taylor_radius():
    N = 400
    grid = Grid(N, 1.0)
    E0, rho0, t = 1.0, 1.0, 0.05
    W0 = np.vstack([np.full(N, rho0), np.zeros(N), np.full(N, 1e-8)])
    nb = 4
    W0[2, :nb] = (g - 1) * E0 / grid.vol[:nb].sum()
    for order in (1, 2):
        sim = Simulation(grid, W0, gamma=g, order=order)
        sim.run_to(t)
        W = sim.primitives()
        Rs = grid.rc[np.argmax(W[0])]
        R_exact = 1.15167 * (E0 * t ** 2 / rho0) ** 0.2
        assert abs(Rs / R_exact - 1) < 0.04, (order, Rs, R_exact)
        # total energy conserved (nothing reaches the boundary)
        Etot = np.sum(sim.U[2] * grid.vol)
        assert abs(Etot / E0 - 1) < 1e-6


def test_injection_conserves_mass_and_energy():
    N = 200
    grid = Grid(N, 1.0)
    W0 = np.vstack([np.full(N, 1.0), np.zeros(N), np.full(N, 1e-3)])
    mdot, L = 0.3, 2.0
    for order in (1, 2):
        sim = Simulation(grid, W0, gamma=g, order=order, mdot=mdot, L=L, R_inj=0.05)
        M0 = np.sum(sim.U[0] * grid.vol); E0 = np.sum(sim.U[2] * grid.vol)
        sim.run_to(0.05)
        M1 = np.sum(sim.U[0] * grid.vol); E1 = np.sum(sim.U[2] * grid.vol)
        assert abs((M1 - M0) / (mdot * 0.05) - 1) < 1e-8
        assert abs((E1 - E0) / (L * 0.05) - 1) < 1e-8
