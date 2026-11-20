# Semi-analytic model (SAM) of an adiabatic wind-blown bubble: assumptions and derivations

Code: `sam/` (`units.py`, `cc85.py`, `shell.py`, `model.py`); driver `run_sam.py`.
Run with `venv/bin/python run_sam.py` from the `opus_high/` directory.

## Global assumptions

1. Spherical symmetry, ideal gas with $\gamma = 5/3$, **no radiative cooling, no thermal
   conduction, no gravity, no magnetic fields**. That is, the system is adiabatic everywhere apart from shocks.
2. The gas has mass $\mu_H m_p = 1.4\,m_p$ per H nucleus, so $\rho = 1.4 m_p n_{\rm H}$. I also use
   $\mu = 1.4$ as the mean mass per particle in the temperature, $T = \mu m_p P/(\rho k_B)$. This
   is equivalent to $P/k_B = n_{\rm H} T$. (The task gives only one molecular weight, so I use the same value
   for both. If you want a different particle $\mu$, rescale $T$ by $\mu/1.4$.)
3. The ambient medium is uniform, static, and has temperature $T_{\rm amb} = 100$ K (not given in the task; I picked
   it). Its pressure is negligible for the dynamics: the outer shock has Mach number $\gtrsim 40$. The ambient pressure is used only
   for the $r > R_2$ part of the output profiles.
4. The wind injects mass $\dot M_w$ and energy $L_w$ (thermal plus kinetic) uniformly per unit volume
   inside $R_{\rm inj}$. Injection starts at $t=0$ and stays constant after that. The terminal wind speed is
   $v_\infty = \sqrt{2L_w/\dot M_w} = 1782$ km/s.
5. Units: the SAM works in cgs. `sam/units.py` provides the constants and a `CodeUnits` class
   ([L] = pc, [t] = kyr, [M] = M$_\odot$), which the hydro code uses.

## Structure (inner to outer)

### 1. Free wind, $r < R_1$: Chevalier & Clegg (1985)

The flow-crossing time of this region ($R_1/v_\infty \sim 10^3$ yr) is short compared with the age, so I take the flow to be
steady. With uniform injection rates $q = \dot M/V_{\rm inj}$ and $Q = L/V_{\rm inj}$, the steady equations integrate
to

- mass: $\rho u = q r/3$ for $r<R$, and $\rho u = \dot M/(4\pi r^2)$ for $r>R$;
- energy (Bernoulli): $u^2/2 + c^2/(\gamma-1) = Q/q = v_\infty^2/2$ everywhere, so
  $u = v_\infty \mathcal{M}/\sqrt{\mathcal{M}^2 + 2/(\gamma-1)}$.

Combining these with momentum gives the CC85 implicit relations for the Mach number. The flow becomes sonic exactly at $R_{\rm inj}$:

$$\Big(\tfrac{3\gamma + 1/\mathcal M^2}{1+3\gamma}\Big)^{-\frac{3\gamma+1}{5\gamma+1}}
\Big(\tfrac{\gamma-1+2/\mathcal M^2}{1+\gamma}\Big)^{\frac{\gamma+1}{2(5\gamma+1)}} = \frac{r}{R}\quad (r<R),\qquad
\mathcal M^{\frac{2}{\gamma-1}}\Big(\tfrac{\gamma-1+2/\mathcal M^2}{1+\gamma}\Big)^{\frac{\gamma+1}{2(\gamma-1)}}=\Big(\frac rR\Big)^2\quad(r>R).$$

I solve these with `brentq`. Then $P = \rho c^2/\gamma$ with $c = u/\mathcal M$.

### 2. Inner (wind) shock, $R_1$

The wind meets a hot bubble of uniform pressure $P_b(t)$. I find $R_1$ by requiring that the Rankine–Hugoniot post-shock
pressure of the CC85 wind equals $P_b$, using the full jump conditions with the finite upstream Mach number. The shock
velocity is $\dot R_1 = 0.4 R_1/t$: ram-pressure balance $\rho_w v_\infty^2 \propto R_1^{-2} \propto P_b
\propto t^{-4/5}$ gives $R_1 \propto t^{2/5}$. The upstream velocity in the shock frame is $u_1 = v_w - \dot R_1$, and

$$\frac{\rho_2}{\rho_1} = \frac{(\gamma+1)\mathcal M_1^2}{(\gamma-1)\mathcal M_1^2+2},\qquad
\frac{P_2}{P_1} = \frac{2\gamma\mathcal M_1^2-(\gamma-1)}{\gamma+1}.$$

### 3. Hot shocked-wind bubble, $R_1 < r < R_c$

- **Pressure.** The sound speed is $\gg$ the flow speeds, so the pressure is uniform: $P = P_b(t)$.
- **Velocity** (Weaver et al. 1977). If the pressure is uniform and the flow is adiabatic, then $D\ln P/Dt = -\gamma\nabla\cdot v$.
  Using $P_b\propto t^{-4/5}$, this gives $\frac{1}{r^2}\partial_r(r^2 v) = \frac{4}{5\gamma t} = \frac{12}{25t}$.
  Integrating with the contact condition $v(R_c) = \dot R_c = 3R_c/(5t)$ gives
  $$ v = \frac{4}{25}\frac{r}{t} + \frac{11}{25}\frac{R_c^3}{t\,r^2}. $$
- **Density** comes from each fluid element's entropy history. In $x = r/R_c(t)$ the velocity field gives
  $dx/d\ln t = \frac{11}{25}(x^{-2}-x)$, so $1-x^3 = (1-x_0^3)(t'/t)^{33/25}$ for an element that crossed
  the inner shock at time $t'$, where $x_0 = R_1(t')/R_c(t')$. I use the power-law scalings
  $R_1\propto t^{2/5}$, $R_c\propto t^{3/5}$, and $\rho_2 \propto R_1^{-2} \propto t^{-4/5}$. With these, $x_0^3(t') = x_0^3(t)\,s^{-3/5}$, where $s=t'/t$, and
  $$ 1 - x^3 = s^{33/25} - x_0^3(t)\, s^{18/25}. $$
  I solve this for $s$ by bisection. Adiabatic compression from $P_b(t')$ to $P_b(t)$ then gives
  $$\rho(r,t) = \rho_2(t')\,[P_b(t)/P_b(t')]^{3/5} = \rho_2(t)\, s^{-8/25}.$$
- This hot-bubble model is approximate. The velocity just behind $R_1$ ($\approx 180$ km/s for F100) does not
  match the Rankine–Hugoniot post-shock velocity ($\approx 450$ km/s). In reality there is a short region where the post-shock flow
  decelerates and the pressure is not exactly uniform. The model neglects it.

### 4. Shocked ambient shell, $R_c < r < R_2$: self-similar adiabatic solution

The swept-up mass, $\propto t^{9/5}$, quickly dominates the wind mass. The total energy grows as $E = L_w t$. Dimensional analysis then
gives $R_2 = \xi_2 (L_w/\rho_0)^{1/5} t^{3/5}$, so $\alpha \equiv d\ln R_2/d\ln t = 3/5$. I write
$\lambda = r/R_2$, $v = \dot R_2 V(\lambda)$, $\rho = \rho_0 G(\lambda)$, $P = \rho_0\dot R_2^2\Pi(\lambda)$, and use
$\partial_t|_r = -(\lambda\dot R_2/R_2)\partial_\lambda$ and $\ddot R_2 = \frac{\alpha-1}{\alpha}\dot R_2^2/R_2$. The Euler
equations (continuity, momentum, and entropy conservation $D(P\rho^{-\gamma})/Dt=0$) then become, with $D\equiv V-\lambda$ and
$c_1 = (1-\alpha)/\alpha$:

$$ D G' + G\,(V' + 2V/\lambda) = 0,\qquad D V' + \Pi'/G = c_1 V,\qquad
D\,\Pi'/\Pi - \gamma D\, G'/G = 2c_1 .$$

Solving for $V'$:

$$ V' = \frac{c_1 V D - (\Pi/G)(2c_1 - 2\gamma V/\lambda)}{D^2 - \gamma\Pi/G},\quad
G' = -\frac{G(V'+2V/\lambda)}{D},\quad \Pi' = \frac{\Pi}{D}\left(2c_1-\gamma V' - \frac{2\gamma V}{\lambda}\right).$$

The boundary conditions at $\lambda = 1$ are the strong-shock jumps: $G = 4$, $V = \Pi = 3/4$. I integrate inward with LSODA until the
contact discontinuity, where $V=\lambda$ (fluid moves with the similarity variable). Results:

- $\lambda_c = R_c/R_2 = 0.8584$ and $\Pi_c = 0.5897$, so $P_b = 0.5897\,\rho_0\dot R_2^2$.
- Mass check: $3\int G\lambda^2 d\lambda = 1.00000$, i.e. all swept-up mass is in the shell.
- Shell energy integral: $I_{\rm sh} = \int_{\lambda_c}^1 \lambda^2(\tfrac12 GV^2 + \Pi/(\gamma-1))\,d\lambda = 0.2238$.

### 5. Energy normalisation, $\xi_2$

In the point-source limit ($R_1\to0$), the hot bubble's energy is purely thermal, $\frac{4\pi}{3}R_c^3 P_b/(\gamma-1)$. Then

$$ L_w t = 4\pi\rho_0 \alpha^2 \xi_2^5 \frac{L_w}{\rho_0} t \left[I_{\rm sh} + \frac{\lambda_c^3\Pi_c}{3(\gamma-1)}\right]
\;\Rightarrow\; \xi_2 = 0.8837,$$

which agrees with the adiabatic value of 0.88 in Weaver et al. (1977).

**Finite-$R_1$ correction (used by default).** At the output time, I solve for $\xi_2$ so that
$L_w t = E_{\rm shell} + E_{\rm bubble} + E_{\rm free\,wind}$. Here $E_{\rm bubble}$ is the thermal plus kinetic energy between
$R_1$ and $R_c$, from the profiles above, and $E_{\rm free\,wind}$ is the energy integrated over the CC85 profile out to $R_1$.
The shell shape ($\lambda_c$, $\Pi_c$, profiles) is kept self-similar. The correction is small: $\xi_2 = 0.8829$ (F100),
0.8821 (F10), 0.8805 (F1).

### Outputs

`sam_outputs/<name>.csv` gives the profiles for $0 < r \le 1.5 R_2$ (log- plus linear-spaced radii, which resolve the
100 AU injection region). `sam_outputs/<name>.txt` gives $R_2$ [pc] and $P_b/k_B$ [K cm$^{-3}$].

| Run | $R_1$ [pc] | $R_c$ [pc] | $R_2$ [pc] | $P_b/k_B$ [K cm$^{-3}$] |
|---|---|---|---|---|
| F100 | 1.591 | 6.543 | 7.623 | $2.000\times10^7$ |
| F10 | 3.166 | 10.36 | 12.07 | $5.015\times10^6$ |
| F1 | 6.284 | 16.39 | 19.10 | $1.255\times10^6$ |

### Summary of approximations / caveats

- Steady free wind (CC85).
- Hot bubble at uniform pressure. The velocity near $R_1$ does not match the post-shock value, and the density uses power-law
  histories for $R_1$, $R_c$, $P_b$.
- The shell is exactly self-similar, $\propto t^{3/5}$. The ambient pressure and the finite wind mass are neglected in the shell.
- Ambient $T = 100$ K, and $\mu = 1.4$ for temperatures.
- When the SAM is compared with the simulation (Task 3), it is evaluated with the simulation's injection radius, 0.5 pc.

References: Chevalier & Clegg 1985, Nature 317, 44; Weaver, McCray, Castor, Shapiro & Moore 1977, ApJ 218, 377;
Castor, McCray & Weaver 1975, ApJ 200, L107.
