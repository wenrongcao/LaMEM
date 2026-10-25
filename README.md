Figures and verification inputs for the `temperature_dependent_kappa` pull request
(branch `feature/temperature-dependent-kappa`). Not meant to be merged.

- `verification/steady.dat` (units = geo) and `verification/steady_si.dat` (units = si):
  1D column, 10 km, 0 / 1000 C, run to steady state (40 Myr, 82 steps).
- `verification/plot_verification.jl`: exact steady-state solution and both figures.
  Expects `<dir>/val/` (steady.dat output) and `<dir>/sphere_td/`, `<dir>/sphere_ck/`
  (t43 input with `nstep_out = 1`, `dt_out = 0.02`; the `_ck` copy uses
  `temperature_dependent_kappa = 0` and `k = 2.5`). Run `julia plot_verification.jl <dir>`
  with CairoMakie and LaMEM.jl installed.
