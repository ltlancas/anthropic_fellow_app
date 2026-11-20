# Semi-analytic model for adiabatic wind-blown bubbles: assumptions & derivations

Code: `sam/model.py` (model), `sam/units.py` (constants/unit conversion, cgs internally),
`scripts/run_sam.py` (produces this directory). Run with `.venv/bin/python scripts/run_sam.py`.

## Assumptions

1. **Spherical symmetry, ideal gas, γ = 5/3**. No radiative cooling, no thermal conduction,
   no gravity, no magnetic fields (fully adiabatic except at shocks).
2. **Composition**: ρ = 1.4 m_p n_H everywhere. As instructed ("all gas … mean molecular weight
   1.4 m_p"), the same μ = 1.4 is used in the equation of state, T = μ m_p P/(ρ k_B), so
   P/k_B = n_H T. (A physically ionised gas would have μ ≈ 0.6; temperatures here scale with μ.)
3. **Wind source**: constant Ṁ_w and L_w deposited uniformly (per unit volume, as mass and
   thermal energy, at rest) inside R_inj = 100 AU. Terminal speed v_∞ = (2L_w/Ṁ_w)^{1/2}
   (1782 km/s for F runs, 178 km/s for S runs).
4. **Ambient medium**: uniform ρ_0, at rest, cold for the dynamics (strong outer shock). Only for
   filling the ambient pressure/temperature outside R_2 in the output files a temperature
   T_amb = 10 K is used.
5. **Energy-driven self-similar regime**: the wind mass and the energy stored in the free-wind
   region are negligible compared with those of the bubble, i.e. t ≫ t_ch where
   R_ch = (Ṁ_w/(ρ_0 v_∞))^{1/2}, t_ch = R_ch/v_∞ (time when swept-up mass ≈ wind mass), and
   R_1 ≪ R_c. Then R_2 ∝ (L t³/ρ_0)^{1/5} exactly.
6. **Reverse-shock layer is quasi-steady** (R_1 moves at (2/5)R_1/t ≪ v_∞), and the reverse shock
   is treated as stationary in the jump conditions.
7. Time t is measured from wind turn-on, with the bubble taken to be self-similar from t = 0.

### Validity of the six experiments (t = 100 kyr)

| run | t/t_ch | R_1/R_c | status |
|---|---|---|---|
| F100 | 1447 | 0.26 | valid |
| F10 | 457 | 0.33 | valid |
| F1 | 145 | 0.42 | valid |
| S100 | 4.6 | 0.83 | marginal (finite wind mass non-negligible) |
| S10 | 1.45 | 1.05 | **invalid** |
| S1 | 0.46 | 1.32 | **invalid** |

For S10 and S1 the self-similar bubble pressure cannot confine the slow, massive wind inside the
contact discontinuity. Here the code falls back to placing the termination shock at R_c, with no
hot shocked-wind region: free wind for r < R_c, then the self-similar shell. In the
wind-dominated limit (t ≪ t_ch) the true solution approaches R_2 ≈ v_∞ t (≈ 18 pc at 100 kyr),
so the reported R_2 and pressure for S10/S1 (and to a lesser degree S100) should be considered
order-of-magnitude only. A model with finite wind mass, needed for t ≲ t_ch, was not built.

## Derivations

### 1. Free wind (r < R_1): Chevalier & Clegg (1985)
With uniform injection q = Ṁ/(4πR_inj³/3) and energy per unit mass v_∞²/2, steady flow gives
mass flux ρ v r² = q r³/3 (r < R_inj) or Ṁ/4π (r > R_inj), and the Bernoulli relation
v²/2 + c²/(γ−1) = v_∞²/2 everywhere. Transonic at r = R_inj, the Mach number M(r) satisfies

- r < R_inj: (r/R_inj) = [(3γ+1/M²)/(1+3γ)]^{−(3γ+1)/(5γ+1)} [(γ−1+2/M²)/(1+γ)]^{(γ+1)/(2(5γ+1))}
- r > R_inj: M [((γ−1)M²+2)/(γ+1)]^{−(γ+1)/(2(γ−1))} = (R_inj/r)²

(the second follows from isentropy ρ ∝ c^{2/(γ−1)}, ρ v r² = const and Bernoulli).
Then c² = v_∞²/(M² + 2/(γ−1)), v = Mc, ρ from the mass flux, P = ρc²/γ.

### 2. Shocked ambient shell (R_c < r < R_2): similarity solution
Write ξ = r/R_2, v = Ṙ_2 U(ξ), ρ = ρ_0 G(ξ), P = ρ_0 Ṙ_2² Π(ξ), with R_2 ∝ t^{3/5} so that
δ ≡ R_2 R̈_2/Ṙ_2² = −2/3. With W = U − ξ and c² = γΠ/G, the continuity, momentum and entropy
equations become

- W G' + G(U' + 2U/ξ) = 0
- δ U + W U' = −Π'/G
- 2δ + W Π'/Π − γ W G'/G = 0

which solve to U' = [2c²U/ξ + 2δΠ/G − δUW]/(W² − c²), G' = −G(U'+2U/ξ)/W,
Π' = Π[−γ(U'+2U/ξ) − 2δ]/W. Strong-shock conditions at ξ = 1: U = 2/(γ+1), G = (γ+1)/(γ−1),
Π = 2/(γ+1). The equations are integrated inward (DOP853, rtol 1e−11) until U = ξ, which is the
contact discontinuity (the piston): **ξ_c = 0.85839, Π_c = 0.58966**.

### 3. Global energy and the radius constant
The hot bubble (R_1 < r < R_c) is at uniform pressure P_b = Π_c ρ_0 Ṙ_2² with negligible kinetic
energy, and R_1³ ≪ R_c³. Energy conservation, L t = E_shell + E_bubble:

L t = 4π ρ_0 R_2³ Ṙ_2² J,  J = ∫_{ξ_c}^1 ξ²(G U²/2 + Π/(γ−1)) dξ + ξ_c³ Π_c/(3(γ−1)).

With R_2 = A (L/ρ_0)^{1/5} t^{3/5}, R_2³Ṙ_2² = (9/25) A⁵ t, so A = [25/(36π J)]^{1/5}.
Numerically **A = 0.88367** (cf. 0.88 in Weaver et al. 1977, eq. 21). As a check, the bubble
fraction of the energy is E_b/(L t) = 0.454545 = 5/11, which is the exact result from
dE_b/dt = L − P_b dV_b/dt with E_b = (3/2)P_b V_b, P_b ∝ t^{−4/5}, V_b ∝ t^{9/5}.

Outputs: R_bubble ≡ R_2 (outer shock radius), P_hot/k_B = Π_c ρ_0 Ṙ_2²/k_B, Ṙ_2 = 3R_2/(5t).

### 4. Reverse shock and quasi-steady layer behind it
The pre-shock state at R_1 is CC85. The Rankine–Hugoniot jump for a stationary shock with
upstream Mach M_1 gives (ρ_s, v_s, P_s). Behind it, the flow is subsonic (M_s² ≈ 1/5), and
on a scale ~R_1 it is steady and adiabatic: ρ v r² = const, entropy K = P_s/ρ_s^γ constant, and
v²/2 + γP/((γ−1)ρ) = const. It therefore decelerates and approaches the stagnation pressure
P_st = P_s (1 + (γ−1)M_s²/2)^{γ/(γ−1)} (≈ 1.175 P_s). **R_1 is found from P_st(R_1) = P_b.**

### 5. Hot bubble interior (Lagrangian/entropy-history solution)
Away from R_1 the pressure is uniform, P_b(t) ∝ t^{−4/5}. A wind parcel shocked at time t_s
keeps its entropy, so its density now is ρ = ρ_st(t_s) (P_b(t)/P_b(t_s))^{1/γ}, where
ρ_st ∝ P_s/v_∞² is the stagnation density. With x ≡ t_s/t and strong-shock
ρ_s = 16 P_s/(3v_∞²), P_s = P_b/f, f = P_st/P_s = (16/15)^{5/2}:

ρ(x) = (5 P_b/v_∞²) x^{−8/25}.

The volume occupied by gas shocked after t_s is ∫ Ṁ dt'/ρ ∝ (1 − x^{33/25}). Normalising to the
bubble volume ΔV = (4π/3)(R_c³ − R_1³) gives the mapping
x(r) = [1 − (4π/3)(r³ − R_1³)/ΔV]^{25/33}. Consistency check: ΔV = (25/33) Ṁ t ρ_st^{−1} is
equivalent to E_b = 5Lt/11, the same value as the shell solution gives, so the two pieces agree.
Differentiating the Lagrangian volume coordinate at fixed t_s (V_1 ∝ t^{6/5}, ΔV ∝ t^{9/5}) gives

4π r² v = dV_1/dt + (ΔV/t) [ (9/5)(1 − x^{33/25}) + (33/25) x^{33/25} ],

which equals 4πR_c² Ṙ_c at the contact, matching the shell.

**Composite solution** in R_1 ≤ r < R_c: P = P_layer(r), ρ = ρ_layer(r) x^{−8/25},
v = v_layer(r) + [v_outer(r) − v_outer(x→1)]. The (integrable) density divergence at R_c
(x → 0) comes from the earliest, lowest-entropy shocked wind. Physically this gas was shocked
before the flow became self-similar, so it is not reliable there.

T is computed from item 2 of the assumptions in every region.
