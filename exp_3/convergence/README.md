# Convergence results (F100, sim vs SAM)

Metric: relative L2 = ||q_sim − q_SAM||₂ / ||q_SAM||₂ over all grid points (uniform grid), averaged
over the 5 snapshots (20–100 kyr). Columns: N, first order, second order.

- `*_L2_norm.txt`: required outputs, using the whole domain.
- `*_L2_norm_excl_r_lt_0.1pc.txt`: supplementary, r ≥ 0.1 pc only. Added after the results
  were seen, because the pressure and density norms over the whole domain are dominated (>90% of
  the squared error) by the first one or two cells. There, the 100 AU injection region is smaller
  than one cell (Δr ≥ 290 AU), so the simulation's central free wind cannot match the resolved
  CC85 profile of the SAM. That full-domain error falls steadily with N (pressure: 161 → 8).
- `per_snapshot_L2*.txt`: per-snapshot values.

Interpretation: outside the core, errors are 3–20% and are dominated by discontinuities: the
reverse shock (~1.7 pc), the contact (~6.5 pc), where the SAM density has an integrable spike,
and the outer shock. For discontinuous solutions L2 converges slowly, at best ~N^{-1/2}. There
is also a floor from SAM model error (self-similar assumptions, steady reverse-shock layer). So
monotone first- or second-order convergence rates should not be expected. Second order gives
sharper discontinuities, so a small offset in shock position costs more under L2. That is why
its velocity and temperature norms come out larger than first order here.

Spot check at 100 kyr, N=2048, second order: outer shock at 7.636 pc (SAM 7.629 pc) and bubble
P/k = 1.999e7 K cm⁻³ (SAM 2.004e7).
