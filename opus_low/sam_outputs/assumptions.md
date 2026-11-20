# Semi-analytic model: assumptions and derivation

Model follows the adiabatic (energy-conserving) phase of Weaver et al. (1977, ApJ 218, 377).

## Assumptions
- Spherical symmetry; ideal gas with gamma = 5/3; no cooling, conduction, gravity or magnetic fields.
- Uniform ambient density rho0 = 1.4 m_p n_H; ambient pressure neglected (strong outer shock).
- Mass and energy injected steadily (Mdot, L) inside r_inj; wind terminal speed v_w = sqrt(2L/Mdot) (all injected energy becomes wind kinetic energy). Inside r_inj the wind velocity is taken to rise linearly to v_w.
- Self-similar evolution, t >> time for wind to cross r_inj; shocked wind mass is negligible compared to swept-up mass.
- Temperature T = P mu m_p/(rho k_B) with mu = 1.4 (mean mass per particle taken equal to the mass per H nucleus, as specified).

## Derivation
Energy: E_th(bubble) + E_kin(shell) = (5/11) L t + (15/77) L t for the self-similar solution; with momentum equation of the thin shell
d/dt[(4/3)pi R^3 rho0 dR/dt] = 4 pi R^2 P and P = (gamma-1) E_th / ((4/3) pi R^3), assuming R ∝ t^{3/5}, gives
- R2(t) = (125/(154 pi))^{1/5} (L t^3/rho0)^{1/5},
- V2 = dR2/dt = (3/5) R2/t,
- P_hot = 7/(3850 pi)^{2/5} L^{2/5} rho0^{3/5} t^{-4/5}.

Structure (inner to outer):
1. Free wind, r < R1: rho = Mdot/(4 pi r^2 v_w), v = v_w, negligible pressure.
2. Reverse shock R1: strong-shock jump gives post-shock P = (3/4) rho_w v_w^2; setting this equal to P_hot gives R1 = sqrt(3 Mdot v_w/(16 pi P_hot)).
3. Shocked wind R1 < r < Rc: isobaric at P_hot, density 4 rho_w(R1) (post-shock value, held uniform), velocity linear from v_w/4 to the contact speed.
4. Shell Rc < r < R2: strong-shock jump conditions, rho = 4 rho0, v = (3/4)V2, P = (3/4) rho0 V2^2. Mass conservation (swept mass in density 4 rho0) gives Rc = (3/4)^{1/3} R2.
5. Ambient r > R2: rho0, v = 0, P = 0.

The reported hot-gas pressure is P_hot; the reported outer radius is R2.
