"""Validation tests for the wind-bubble model.  Run with:  .venv/bin/pytest -q"""
import numpy as np
import pytest

from wind_bubble import WindBubble, units as u, weaver

T_FID = 1e5 * u.yr


@pytest.fixture(scope="module")
def fiducial():
    return WindBubble(1e38, 1e-2 * u.Msun / u.yr, 100.0).solve(T_FID)


def test_units():
    assert np.isclose(u.mean_mass_per_particle(1.4, state="ionized"), 1.4 / 2.3)
    assert np.isclose(u.mean_mass_per_particle(1.4, state="atomic"), 1.4 / 1.1)
    wb = WindBubble(1e38, 1e-2 * u.Msun / u.yr, 100.0)
    assert np.isclose(wb.rho0, 140.0 * u.m_p)
    assert np.isclose(wb.v_w, np.sqrt(2e38 / (1e-2 * u.Msun / u.yr)))


def test_weaver_limit():
    """Light, fast wind at t >> t_sw must reproduce Weaver et al. (1977)."""
    wb = WindBubble(1e36, 1e-7 * u.Msun / u.yr, 100.0).solve(1e4 * u.yr)
    c = wb.weaver_comparison(1e4 * u.yr)
    assert abs(c["Rc"] - 1) < 1e-3
    assert abs(c["V"] - 1) < 1e-3
    assert abs(c["P_b"] - 1) < 1e-3
    assert abs(c["E_b"] - 1) < 1e-2
    # R1 differs from Weaver only by the stagnation-pressure factor
    assert abs(c["R1"] - (16.0 / 15.0) ** 1.25) < 5e-3


def test_initial_condition_insensitivity():
    ref = WindBubble(1e38, 1e-2 * u.Msun / u.yr, 100.0).solve(T_FID).state(T_FID)
    for frac, r1f in [(1e-3, 0.5), (1e-5, 0.9)]:
        st = WindBubble(1e38, 1e-2 * u.Msun / u.yr, 100.0).solve(
            T_FID, t_start_frac=frac, R1_frac0=r1f).state(T_FID)
        for k in ("R1", "Rc", "R2", "P_b", "V"):
            assert abs(getattr(st, k) / getattr(ref, k) - 1) < 1e-4, k


def test_global_conservation(fiducial):
    wb = fiducial
    st = wb.state(T_FID)
    # energy: L t = KE(free wind) + E_b + E_kin,shell + work on ambient (P0 term)
    E_fw = 0.5 * wb.Mdot_w * st.R1 / wb.v_w * wb.v_w**2
    E_shell = 0.5 * (4 * np.pi / 3 * wb.rho0 * st.Rc**3) * st.V**2
    # radiated energy of the isothermal shell is not tracked -> E_tot <= L t
    assert E_fw + st.E_b + E_shell < wb.L_w * T_FID
    # mass: injected = free wind + shocked wind
    assert np.isclose(st.M_b + wb.Mdot_w * st.R1 / wb.v_w, wb.Mdot_w * T_FID, rtol=1e-12)


def test_profile_consistency(fiducial):
    p = fiducial.profiles(T_FID)
    d = p["diagnostics"]
    # Lagrangian reconstruction closes: velocity at Rc equals shell velocity,
    # and the bubble volume is filled by shocked wind
    assert abs(d["v_at_Rc_over_V"] - 1) < 1e-8
    assert d["V_ic_fraction"] < 1e-6
    # shocked-wind mass from the profile equals M_b
    st = p["state"]
    m = p["region"] == 1
    r, rho = p["r"][m], p["rho"][m]
    lag = fiducial._shocked_wind_lagrangian(T_FID / fiducial.t_sw)
    rl = lag["r"] * fiducial.l_sc
    rhol = lag["rho"] * fiducial.m_sc / fiducial.l_sc**3
    M = np.trapezoid(4 * np.pi * rl**2 * rhol, rl)
    assert abs(M / st.M_b - 1) < 1e-3
    # shell contains all swept-up mass
    m2 = p["region"] == 2
    r2 = np.linspace(st.Rc, st.R2, 2001)
    rho2 = np.interp(r2, p["r"][m2], p["rho"][m2])
    Ms = np.trapezoid(4 * np.pi * r2**2 * rho2, r2)
    assert abs(Ms / st.M_s - 1) < 1e-2
    # ordering and pressure equilibrium across contact discontinuity
    assert 0 < st.R1 < st.Rc < st.R2
    i = np.searchsorted(p["r"], st.Rc)
    assert np.isclose(p["P"][i - 1], p["P"][i], rtol=1e-6)
    # ambient medium beyond R2
    assert np.allclose(p["n_H"][p["r"] > st.R2 * 1.001], 100.0)


def test_reverse_shock_jump(fiducial):
    """Mass flux across R1 is conserved and post-shock gas is subsonic."""
    st = fiducial.state(T_FID)
    p = fiducial.profiles(T_FID, r=[st.R1 * (1 - 1e-9), st.R1 * (1 + 1e-9)])
    rho_w, rho_b = p["rho"]
    v_w, v_b = p["v"]
    assert np.isclose(rho_w * (v_w - st.dR1dt), rho_b * (v_b - st.dR1dt), rtol=1e-6)
    cs = np.sqrt(fiducial.gamma * st.P_b / rho_b)
    assert abs(v_b - st.dR1dt) < cs
