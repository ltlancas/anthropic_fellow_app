# Semi-analytic model: assumptions and derivations

Code: `sam/` (python package; `sam/bubble.py` holds the physics, `sam/units.py` the CGS constants).
Driver: `run_sam.py`.

## Assumptions
1. **Spherical symmetry, ideal gas, γ = 5/3, fully adiabatic.** No radiative cooling, thermal conduction, gravity, magnetic fields or mixing. That means the swept-up shell is *not* thin. It is a hot, adiabatic layer of shocked ISM.
2. **Composition.** ρ = 1.4 m_p n_H. For the temperature I take the same 1.4 m_p as the mean mass per particle, so T = P · 1.4 m_p / (ρ k_B) = P / (n_H k_B). The task only specifies "1.4 m_p per H nucleus", so this is an interpretation. A fully ionised H+He gas would have μ ≈ 0.6 instead.
3. **Wind.** Mass (Ṁ_w) and energy (L_w) are injected uniformly in a sphere of radius R_inj (100 AU for Task 1) with no momentum injection. The terminal wind speed is v_∞ = (2 L_w/Ṁ_w)^{1/2} = 1782 km/s.
4. **Ambient medium.** Uniform, at rest, T_amb = 100 K (P/k = 10⁴ K cm⁻³ for n_H = 100). The task leaves the ambient temperature unspecified. The ambient pressure is neglected in the dynamics (strong-shock limit), since P_amb/P_bubble ≲ 10⁻³. The ambient values are only used to fill r > R₂ in the output files.
5. **Energy budget.** L_w t = (thermal energy of the shocked wind) + (kinetic + thermal energy of the shell). The energy in the free-wind region (≈ L_w R₁/v_∞ ≲ 10⁻³ L_w t) and the volume of the free-wind region inside the bubble (R₁³/R_c³ ≈ 1.5%) are neglected when fixing the bubble pressure. Neglecting them is what makes the outer solution self-similar.
6. **Shocked wind.** Its sound speed (≈ 1000 km/s) is far larger than the expansion speed (≈ 45 km/s), so its pressure is taken as spatially uniform, P_b(t) (Weaver et al. 1977).

## Structure (inside → out)
| Region | Radii | Solution |
|---|---|---|
| (a) free wind | r < R₁ | Chevalier & Clegg (1985) steady wind with uniform injection |
| (b) shocked wind | R₁ < r < R_c | uniform P_b; Weaver velocity field; density from adiabatic Lagrangian history |
| (c) shocked ISM shell | R_c < r < R₂ | self-similar adiabatic solution (numerically integrated ODEs) |
| (d) ambient | r > R₂ | ρ₀, v = 0, T_amb |

## (c) Self-similar shell and outer shock
Dimensional analysis with constant L and uniform ρ₀ gives R₂ = ξ₀ (L t³/ρ₀)^{1/5}. The shock speed is U = Ṙ₂ = (3/5) R₂/t, so U ∝ t^{-2/5} and dU/dt = −(2/3) U²/R₂.
Write λ = r/R₂, v = U V(λ), ρ = ρ₀ G(λ), P = ρ₀ U² Π(λ). Substituting into the Euler equations (∂/∂t at fixed r = −λ (U/R₂) d/dλ, plus the explicit time dependence of U) gives, with w = V − λ:

- continuity: w G' + G V' = −2GV/λ
- momentum: w V' + Π'/G = (2/3) V
- entropy, D(Pρ^{-γ})/Dt = 0: w Π'/Π − γ w G'/G = 4/3

The strong-shock jump conditions at λ = 1 are G = (γ+1)/(γ−1) = 4, V = Π = 2/(γ+1) = 3/4.
I integrate the system inward with `solve_ivp` (rtol 10⁻¹¹), solving the 3×3 linear system at each step. The integration stops at the contact discontinuity, where the fluid moves with the self-similar frame: V(λ_c) = λ_c. Result: **λ_c = R_c/R₂ = 0.8584**, **Π_c = 0.5897**.

Energy normalisation: the bubble thermal energy is E_b = (3/2) P_b (4π/3) R_c³ with P_b = ρ₀ U² Π_c, and E_shell = 4π R₂³ ρ₀ U² ∫_{λ_c}^{1} [½ G V² + Π/(γ−1)] λ² dλ. Then
L t = E_b + E_shell ≡ K ρ₀ U² R₂³ = (9/25) K ξ₀⁵ L t, so ξ₀ = (25/(9K))^{1/5}.

Numerically **ξ₀ = 0.8837**, which agrees with Weaver et al. (1977), R₂ = 0.88 (L t³/ρ₀)^{1/5}. The energy split is 45.45% bubble thermal (exactly 5/11, as in Weaver et al.), 29.7% shell thermal and 24.9% shell kinetic. This is a check on the integration.

**Hot-gas pressure** (reported in `<exp>.txt`): P_b = ρ₀ U² Π_c = 0.5897 ρ₀ (0.6 R₂/t)².

## (a) Free wind: Chevalier & Clegg (1985)
Inside R_inj: q_m = Ṁ/V_inj and q_E = L/V_inj. Outside R_inj: an adiabatic steady wind. The Mach number M(r) comes from the CC85 implicit relations, solved with `brentq` (subsonic inside, M = 1 at R_inj, supersonic outside):

- r < R_inj: [(3γ+1/M²)/(1+3γ)]^{−(3γ+1)/(5γ+1)} [(γ−1+2/M²)/(1+γ)]^{(γ+1)/(2(5γ+1))} = r/R_inj
- r > R_inj: M^{2/(γ−1)} [(γ−1+2/M²)/(1+γ)]^{(γ+1)/(2(γ−1))} = (r/R_inj)²

Bernoulli gives v²/2 + c²/(γ−1) = v_∞²/2, so v = v_∞ [1 + 2/((γ−1)M²)]^{−1/2}. Mass conservation gives ρv = Ṁ r/(4πR_inj³) inside R_inj and ρv = Ṁ/(4πr²) outside, and P = ρ c²/γ with c = v/M.

## Reverse shock R₁
R₁ is the radius where the post-shock pressure of the free wind equals P_b. The full Rankine–Hugoniot conditions are applied in the shock frame, with shock speed Ṙ₁ = (2/5) R₁/t. This follows from ρ_w v_∞² ∝ R₁⁻² balancing P_b ∝ t^{−4/5}. The root is found with `brentq`.

## (b) Shocked wind
**Pressure.** Uniform, P = P_b(t) ∝ t^{−4/5}.

**Velocity** (Weaver et al. 1977, eq. 6). The flow is adiabatic, so D ln ρ/Dt = (3/5) d ln P_b/dt = −12/(25t). Continuity then gives (1/r²) ∂(r²v)/∂r = 12/(25t), so v = (4/25) r/t + C/r². The boundary condition v(R_c) = Ṙ_c = (3/5)R_c/t fixes C = (11/25) R_c³/t:

  v(r,t) = 11 R_c³/(25 t r²) + 4 r/(25 t).

**Density.** With x = r/R_c(t), pathlines satisfy dx/d ln t = (11/25)(1 − x³)/x², which integrates to (1 − x³) t^{33/25} = const along a pathline. For a point (r, t), I solve for the time t' at which that fluid element crossed the reverse shock, i.e. x(t') = R₁(t')/R_c(t'). I then take the exact post-shock density ρ_s(t') from the jump conditions. Adiabatic compression at uniform pressure then gives

  ρ(r,t) = ρ_s(t') (t/t')^{−12/25}.

**Accuracy.** The quasi-static region (b) is approximate. Its velocity field does not exactly match the post-shock velocity at R₁, and the R₁ volume is ignored in the energy budget. As a result, the total mass in (a)+(b) at 100 kyr for F100 is 1.07 × Ṁt (7% excess).

## Results at t = 100 kyr
| Exp | R₁ [pc] | R_c [pc] | R₂ [pc] | P_b/k [K cm⁻³] |
|---|---|---|---|---|
| F100 | 1.589 | 6.549 | 7.629 | 2.004e7 |
| F10 | 3.160 | 10.380 | 12.092 | 5.033e6 |
| F1 | 6.262 | 16.451 | 19.164 | 1.264e6 |

## Output grid
The radii are the union of 600 log-spaced points (10⁻⁶ pc to 1.5 R₂) and 1500 linear points (10⁻³ R₂ to 1.5 R₂).

## References
- Castor, McCray & Weaver 1975, ApJ 200, L107
- Weaver, McCray, Castor, Shapiro & Moore 1977, ApJ 218, 377
- Chevalier & Clegg 1985, Nature 317, 44
