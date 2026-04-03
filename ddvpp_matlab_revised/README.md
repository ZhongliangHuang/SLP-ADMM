# DDVPP IEEE 118 MATLAB case package

## Included features

This package implements a complete case-study scaffold for the DDVPP framework:

1. Wang-style nodal frequency state-space model on generator buses with full-bus reconstruction
2. worst disturbance bus scan
3. critical oscillatory mode selection with one mode per conjugate pair
4. MAC-based modal tracking across SLP outer iterations
5. analytic eigenvalue real-part sensitivity with respect to IBR virtual inertia and damping
6. adaptive trust-region update using the ratio between actual and predicted modal push
7. inner ADMM with residual balancing
8. diagnostics for modal normalization, finite-difference sensitivity error, QSS, and RoCoF consistency
9. figure export, including spatial RoCoFmax distribution

## Files

- `run_ddvpp_case118_demo.m` main entry
- `setup_ddvpp_case118_ddvpp.m` deterministic case initialization
- `ddvpp_slp_admm_case118.m` nested outer/inner solver
- `ddvpp_diagnostics.m` self-checks that expose likely implementation errors
- `ddvpp_*.m` helper functions

## How to run

1. Add MATPOWER to the MATLAB path.
2. Put this folder on the MATLAB path.
3. Run:

```matlab
results = run_ddvpp_case118_demo();
```

4. Outputs are saved to:

```matlab
./ddvpp_outputs/
```

## Output figures

- `fig1_frequency_trajectories.png`
- `fig2_spatial_nadir.png`
- `fig3_eigenvalues.png`
- `fig4_trust_region.png`
- `fig5_admm_balance.png`
- `fig6_final_allocation.png`
- `fig7_rocof_spatial.png`
- `summary.txt`

## Key corrections in this version

- The state matrix now uses `omega0 * I` in the angle-frequency coupling block.
- `J`, `L`, `F`, `Dtilde`, and `N` follow the Wang model structure consistently.
- Frequency output is reconstructed to all buses via the frequency divider.
- The modal envelope uses `|(u^H B d) * v_omega / lambda|`, which is consistent with the spectral formula.
- QSS requirement is converted to a damping deficit, not a frequency-error deficit.
- RoCoF constraints include the conversion from pu/s to Hz/s.
- Simulation uses a block exponential and does not require `pinv(A)`.
- Diagnostics compare analytic sensitivities against finite differences.
