# Model assumptions

This file lists every assumption behind the semi-analytic wind-bubble model in
`wind_bubble/`. Section 7 states where these assumptions break down for the
fiducial parameters requested (L_w = 1e38 erg/s, Mdot = 1e-2 Msun/yr,
n_H = 100 cm^-3, t = 1e5 yr).

## Literature basis

- Castor, McCray & Weaver 1975, ApJL 200, L107: thin-shell, energy-driven bubble.
- Weaver, McCray, Castor, Shapiro & Moore 1977, ApJ 218, 377: the four-zone
  structure and the self-similar solution R2 = 0.763 (L t^3/rho0)^(1/5) and
  P = 7/(3850 pi)^(2/5) L^(2/5) rho0^(3/5) t^(-4/5). This solution is
  implemented in `wind_bubble/weaver.py` and is the late-time limit of this model.
- Koo & McKee 1992, ApJ 388, 93: energy-driven versus momentum-driven bubbles,
  including the role of the swept-up mass and the reverse shock.
- Draine 2011, *Physics of the Interstellar and Intergalactic Medium*, eq. 34.4:
  the cooling-function fit for 1e5 K < T < 2e7 K. Above 2e7 K it is joined
  continuously to a free-free branch, Λ ∝ T^(1/2). It is used only for a
  diagnostic.

How this model differs from Weaver et al. (1977): it does **not** assume
R1 << R2, and it does **not** assume the wind mass is negligible next to the
swept-up mass. Both are wrong for the fiducial wind. Its speed is
v_w ≈ 178 km/s, the swept-up mass equals the injected mass at
t_sw ≈ 1.1e4 yr, so t = 1e5 yr is only about 9 t_sw, and R1 ≈ 0.77 Rc.
The model integrates the full time evolution from the free-expansion phase
onward, so it is not self-similar.

## 1. Geometry and inputs

1. The flow is strictly spherically symmetric and one-dimensional.
2. Mass and energy are injected at r = 0 at constant rates, starting at t = 0.
3. The injected energy is entirely the wind's kinetic energy, so
   v_w = sqrt(2 L_w / Mdot). The free wind is launched with negligible thermal energy.
4. The ambient medium is uniform, static and infinite, with
   rho0 = mu_H m_p n_H, where mu_H = 1.4. Gravity, magnetic fields,
   turbulence and the star's own motion are ignored.

## 2. Composition and units

5. The gas is H plus He only, with helium abundance x_He = (mu_H - 1)/4 = 0.1
   by number (consistent with mu_H = 1.4).
6. The free wind, shocked wind and shell are fully ionized, with mean mass per
   particle mu = 1.4/2.3 ≈ 0.609 m_p. The ambient medium is neutral atomic,
   with mu = 1.4/1.1 ≈ 1.27 m_p. The ambient state can be changed through the
   `ambient_state` argument.
7. All internal calculations are in CGS. `wind_bubble/units.py` defines the
   constants (CODATA 2018 / IAU 2015) and the conversions; a year is the Julian
   year. Multiplying by a unit converts into CGS; dividing converts out.
8. All gas has an adiabatic index of γ = 5/3.

## 3. Zone structure

9. **Free wind (r < R1).** The wind is cold, highly supersonic and ballistic:
   rho = Mdot/(4π r² v_w) and v = v_w. Its temperature T_wind (default 1e4 K,
   i.e. photoionized) is used only to report a pressure and temperature in this
   zone. It does not affect the dynamics because the reverse shock is assumed strong.
10. **Reverse shock (R1).** The shock is strong and adiabatic, with Rankine–Hugoniot
    compression ratio (γ+1)/(γ-1) = 4 in the shock frame.
11. **Shocked wind (R1 < r < Rc).**
    - The bubble is **isobaric**: its pressure P_b(t) is uniform. This holds
      because its sound speed exceeds its flow speeds; for the fiducial case
      c_s ≈ 90 km/s against flow speeds of 40–60 km/s.
    - Its energy is all thermal, E_b = P_b V_b/(γ-1). The bulk kinetic energy
      of the shocked wind is neglected. The error is O(V²/v_w²), which is
      about 5% for the fiducial case.
    - Each fluid element is adiabatic after it is shocked. There is **no
      radiative cooling** and **no thermal conduction or evaporation** from the
      shell. Weaver et al. (1977) include conduction; it is left out here.
12. **Reverse-shock condition.** This is an energy-consistent version of
    Weaver's condition P_b = ¾ ρ_w v_w².
    - Freshly shocked gas has the post-shock entropy from the jump conditions.
      It is then brought isentropically to the bubble pressure, so its residual
      post-shock kinetic energy becomes thermal pressure. This is the
      stagnation-pressure correction, and it makes P_b ≈ 0.88 ρ_w w² rather
      than 0.75 ρ_w w².
    - The shock velocity is then fixed by requiring that the wind energy flux
      across the *moving* shock, L_w(1 − Ṙ1/v_w), equals the enthalpy flux of
      the newly shocked gas plus the work P_b 4πR1² Ṙ1.
    - With this condition, the Lagrangian reconstruction of the interior
      exactly fills the bubble volume and matches the shell velocity at Rc.
      Both are verified to better than 1e-8 in the tests. With Weaver's
      condition the reconstruction would miss about 20% of the volume.
    - In the Weaver limit the result is identical except that R1 is larger by
      a factor (16/15)^(5/4) ≈ 1.084.
13. **Contact discontinuity (Rc).** The pressure is continuous across Rc, and
    there is no mixing between shocked wind and shell gas.
14. **Swept-up shell (Rc < r < R2).**
    - The shell is thin and **isothermal** at T_shell = 1e4 K, so it is
      radiative and photoionized. Its dynamics are those of a thin shell at Rc
      whose mass is all the ambient gas swept up, M_s = (4π/3) ρ0 Rc³.
    - For the profiles, the shell's pressure varies linearly from P_b at Rc to
      the isothermal post-shock value at R2, and its density is ρ = P/c_s².
    - R2 is chosen so that the shell holds exactly (4π/3) ρ0 R2³. The gas
      velocity is interpolated linearly from V at Rc to the post-shock value at R2.
    - The model refuses to build profiles if P_b < ρ0 c_s², because the
      isothermal shell would then not be denser than the ambient gas.
15. **Ambient gas (r > R2).** The gas is undisturbed at T_ambient = 100 K. Its
    pressure P0 is included in the shell's equation of motion.

## 4. Governing equations

With y = (Rc, V, R1, E_b):

- Rc changes as dRc/dt = V.
- The shell's momentum obeys d(M_s V)/dt = 4π Rc² (P_b − P0), with M_s = (4π/3) ρ0 Rc³.
- The bubble's energy obeys dE_b/dt = L_w (1 − Ṙ1/v_w) − 4π Rc² P_b V.
- R1 changes as dR1/dt = v_w − w(P_b/ρ_w(R1)), where w is the solution of the
  shock condition in item 12.

These are integrated with an implicit Radau solver, in units where
v_w = Mdot = 1 and the time unit is t_sw. Integration starts at t0 = 1e-4 t_sw
in the free-expansion phase, with Rc = v_w t0, V = v_w, R1 = 0.9 Rc and a
pressure consistent with Ṙ1 = 0.9 v_w. The results do not depend on these
starting values: varying t0 over 1e-5 to 1e-3 t_sw and R1/Rc over 0.5 to 0.9
changes all global quantities at 1e5 yr by less than 1e-4 (a tested property).

## 5. Interior profiles

- **Density.** Gas shocked at time t′ has density
  ρ(t) = ρ_eff(t′) [P_b(t)/P_b(t′)]^(1/γ) at the current time. These elements
  are stacked outward from R1, with the most recently shocked gas innermost,
  to give r(t′). Temperature follows from T = μ m_p P_b/(ρ k).
- **Velocity.** The velocity comes from the continuity equation for an isobaric,
  adiabatic flow:
  v(r) = [R1² Ṙ1 + Ṁ_sh/(4π ρ_eff) − (Ṗ_b/(3γ P_b))(r³ − R1³)] / r².
- **Cold layer at Rc.** A thin, cold, dense layer appears just inside Rc. It is
  wind that was shocked weakly during the early free-expansion phase (t ≲ t_sw),
  and it contains about 8% of the shocked-wind mass. This is what the adiabatic
  model predicts. In reality this gas would cool and merge with the shell.
- **Time dependence.** Profiles are available at any time t up to the
  integration end, through `WindBubble.profiles(t)`. Global histories are
  available through `WindBubble.history(t_array)`.

## 6. Numerics

- The dense ODE output is resampled on about 40,000 shock times for the
  Lagrangian integral, which uses the trapezoid rule.
- The ODE tolerances are rtol = 1e-9 and atol = 1e-14 in the dimensionless units.

## 7. Caveats for the fiducial parameters (please read)

- **The adiabatic assumption is badly violated.**
  - The shocked wind is only at T ≈ 1.5–3.7e5 K, near the peak of the cooling
    curve, and has n_H ≈ 25–60 cm⁻³. Its post-shock cooling time is about
    1e3 yr, which is far shorter than 1e5 yr.
  - The a-posteriori cooling luminosity of the adiabatic solution is about
    67 L_w, counting only gas at T > 1e5 K.
  - A real bubble with these inputs would therefore lose its hot interior to
    radiation. Its shocked wind would collapse into a cold, thin layer, and the
    bubble would be **momentum-driven**.
  - Because the wind is so heavily mass-loaded and t is only about 9 t_sw,
    the two regimes have not yet separated in *size*:
    - A fully radiative thin shell that conserves momentum, with the wind mass
      included (d[(M_s + Mdot t)V]/dt = Mdot v_w), reaches 6.56 pc at 1e5 yr.
      The adiabatic model gives Rc = 6.40 pc.
    - The adiabatic shell plus bubble carries about 1.8e5 Msun km/s of
      momentum, close to the total wind momentum Mdot v_w t.
  - The shell radius and velocity are therefore fairly robust. The interior
    state is not: the hot-bubble pressure, temperature and density profiles
    would change drastically.
  - The profiles produced here are the correct *energy-conserving*
    (adiabatic-bubble) solution that was requested. Interpret the interior
    with this in mind.
- Photoionization by the driving star(s), which would create an HII region,
  raise the ambient pressure and set the shell temperature, is not modeled
  except through the fixed T_shell.
- Thermal conduction is omitted. In the Weaver et al. (1977) picture it would
  evaporate shell gas into the bubble, raising the density and lowering the
  temperature near Rc, which would make cooling even stronger.

## 8. Case with Mdot = 1e-4 Msun/yr (other inputs unchanged)

Run this case with `scripts/run_fiducial.py --mdot 1e-4`. Its outputs have the
suffix `_mdot1e-04`.

- **Weaver regime.** The wind speed is v_w ≈ 1780 km/s and t_sw ≈ 34 yr, so
  t ≈ 3000 t_sw. The bubble is firmly in the Weaver et al. (1977) regime:
  - R1 ≈ 0.26 Rc.
  - Rc, V and P_b agree with the self-similar solution to better than 1e-3.
  - E_b is 0.98 of (5/11) L_w t.
- **Adiabatic assumption holds.** The shocked wind is at T ≈ 1e7–5e7 K with
  n_H ≈ 0.2–0.3 cm⁻³. Its post-shock cooling time is about 2e8 yr, and the
  a-posteriori cooling luminosity is about 5e-4 L_w.
- **Conduction would matter here.** Weaver et al. (1977) include thermal
  conduction, which this model omits. With it, the inner edge of the shell
  would evaporate into the bubble. That would raise the bubble's density by
  roughly an order of magnitude near Rc and lower its temperature there.
  The global dynamics would be essentially unchanged.
