# Wind Blown Bubble Grader

## Task Specifications for Grader

A previous instance of Claude has attempted to create a semi-analytic model for the evolution of wind-blown bubbles as a function of time and radius under the assumption of constant injection of mass($\dot{M}_w$) and energy ($L_w$) within a given radius ($r_{\rm fb}$), into a background with a given constant number density of hydrogen nuclei ($n_{\rm H}$) and mean molecular weight of the gas (per hydrogen nucleus, $\mu_{\rm H}$). Below I list a set of experiments that the model should have carried out with corresponding names, parameters (wind mass loss, energy input, and injection radius), and the time at which it was asked to output the profiles.

| Experiment Name | $\dot{M}_w$ [$M_{\odot}/{\rm yr}$] | $L_w$ [${\rm erg}\,{\rm s}^{-1}$]  | $n_{\rm H}$ [${\rm cm}^{-3}$] | $t_{\rm out}$ [${\rm kyr}$]
| --- | --- | --- | --- | --- |
| F100 | $10^{-4}$ | $10^{38}$ | $100$ | $100$ |
| F10 | $10^{-4}$ | $10^{38}$ | $10$ | $100$ |
| F1 | $10^{-4}$ | $10^{38}$ | $1$ | $100$ |
| S100 | $10^{-2}$ | $10^{38}$ | $100$ | $100$ |
| S10 | $10^{-2}$ | $10^{38}$ | $10$ | $100$ |
| S1 | $10^{-2}$ | $10^{38}$ | $1$ | $100$ |

For each of the experiments the injection radius is assumed to be $100\, {\rm AU}$ and the mean molecular weight of the gas is assumed to be $1.4 m_p$ where $m_p$ is the proton mass.

## How to Grade

In the `grader_tools/` directory there is a python file (`ad_wind_spec.py`) that specifies a class `AdiabaticWind()` that takes as arguments

| Argument | Meaning |
| --- | --- |
| `Lwind` | $L_w$ or enery input rate |
| `Mdotw` | $\dot{M}w$ or mass loss rate|
| `rfb` | $r_{\rm fb}$ or injection radius |
| `rho0` | $\rho$ or background mass density |

It has methods `radius`, `pressure`
