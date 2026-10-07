# Conical-flow Python modernization

This is an independent, modular migration of the **angular cone equations** in the archived Project 1 Fortran source. Original sources remain unchanged. The inviscid solver has an independent Taylor–Maccoll reference; the added viscous local-station model has constitutive, analytic-operator and Couette-heating verification, but full axisymmetric viscous-cone benchmark validation remains open.

[AUSM and viscous comparison report](results-v2/report.html) · [Viscous equations and source audit](VISCOUS.md) · [Extension verification](results-v2/verification.json)

The runnable core uses Python and NumPy; Matplotlib generates scientific figures. No packaged CFD solver is used. Comments, documentation and generated reports are in English.

The runnable default is **HLLC/CWENO3**, which reached the residual gate on all six case/grid combinations in the common-CFL study. MUSCL-MC and CWENO5 remain available for comparison, but some shock-containing runs fail to reach a steady state; the report preserves those failures. Smooth fifth-order accuracy alone does not establish shock-case robustness.

## Expanded algorithm and shock-tube edition

The expanded shared core has **15 pointwise Euler fluxes**, **20 FV reconstructions**, and a standalone **JST** central stencil. [Algorithm inventory, sources, scope and executed cone gates](ADDITIONAL_METHODS.md) distinguish new methods from existing aliases and retained failed cases. The cone notebook includes **46 additional runs** alongside the original 68 recorded cases.

A **[separate shock-tube Colab notebook](../shock-tube/README.md)** tests the same shared algorithms with exact transient Riemann references, full combination coverage, conservation budgets and mesh-error comparisons.

## Open the complete 2026 edition

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Ehsan-Roohi/Ehsan-Roohi/blob/main/courses/cfd-2/modern-python/conical/Conical_Flow_Modern_2026.ipynb)

**[One Colab notebook](Conical_Flow_Modern_2026.ipynb)** displays the original 68-case study and 46 additional cases. Expanded comparisons include all 15 pointwise fluxes in both Euler cases and both viscous thermal walls, separate inviscid JST curves, and the added FV reconstructions. The original six-reconstruction, higher-order viscous/grid and DG comparisons remain available. Setup is collapsed, with no encoded source block. Recorded runs display by default; an optional fresh comparison selects all 15 pointwise fluxes by default.

[Published result report with figures and acceptance tables](RESULTS.md) · [DG companion and its limits](DG.md) · [Download source, results and vector figures](../../downloads/conical-modern-2026.zip) · [Notebook execution evidence](colab-validation.json)

The original coursework is retained as an archive. The Python edition is a corrected independent implementation and extension; it does not certify the original executable. Smooth verification and successful steady runs are distinguished from experimental or nonconverged cases.

## Run and compare

```bash
python verify.py
python study.py --cells 60 120 240
python stability_check.py
python render_report.py
python solver.py --mach 7.95 --outer 22 --cells 240 --flux hllc --reconstruction cweno3 --steps 40000 --output results/my-run
```

Run from this directory on a normal host. The restricted desktop environment required launching from the authorized Courses root, with script/output paths relative to that root.

[Measured comparison report](results/report.html) · [Run metadata](results/comparison.json) · [Verification](results/verification.json)

## Equations and nondimensionalization

Let `U = [rho, rho*vr, rho*vtheta, rho*E]`, with velocity normalized by upstream speed, density by upstream density, and pressure by `rho_inf*U_inf**2`. Thus `p_inf = 1/(gamma*M_inf**2)`. This matches the selected Fortran source's ideal-gas nondimensionalization.

The pseudo-time angular balance law is:

```text
dU/dtau + dF/dtheta = S
F = [rho*vt, rho*vr*vt, rho*vt**2+p, (rho*E+p)*vt]
S = [-rho*(2*vr+vt*cot(theta)),
     -rho*(2*vr**2-vt**2+vr*vt*cot(theta)),
     -rho*(3*vr*vt+vt**2*cot(theta)),
     -(rho*E+p)*(2*vr+vt*cot(theta))]
```

This is a conical/axisymmetric angular reduction, not planar wedge flow. Finite volumes store angular cell averages. Each face has one numerical flux shared by its adjacent cells. Source averages use three-point Gaussian quadrature and the reconstructed cell polynomial.

## Algorithms

| Component | Implemented choices |
| --- | --- |
| Numerical flux | Roe with acoustic entropy smoothing, Van Leer, HLLC, HLLE, Rusanov, AUSM, AUSM+, AUSM+-up, AUSM+-up2, SLAU2, entropy-conservative + LF |
| Reconstruction | Piecewise constant, MUSCL-MC, CWENO3/5; characteristic CWENO3/5; nonuniform primitive-minmod for viscous runs |
| Time stepping | CFL-based SSPRK3 toward a steady state |
| Geometric source | Same conical inviscid terms as the source, integrated through the reconstruction |
| Wall | Reflected angular momentum in ghosts; zero mass/energy wall flux and pressure force |
| Outer boundary | Fixed upstream Cartesian-uniform state transformed to the angular basis |
| Safeguards | Sampled-state admissibility scaling; rejected/reduced steps; HLLC face-local HLLE fallback |
| Viscous extension | Sutherland/constant viscosity, heat conduction, stretched wall mesh, FV diffusion of order 2/4, adiabatic/isothermal no-slip walls |
| Optional steady solve | Damped pseudo-transient Newton; colored Jacobian; explicit recovery and independent SSPRK drift check |

CWENO uses one nonlinear polynomial throughout a cell, so it supplies both face states and source quadrature. The code derives interpolation matrices from **cell averages**, not nodal finite differences. CWENO3/5 use the framework of [Cravero et al.](https://arxiv.org/abs/1607.07319); mesh-dependent epsilon follows the considerations in [Cravero and Semplice](https://arxiv.org/abs/1503.00736). This independent implementation uses ordinary Jiang–Shu-style weights, not a claimed implementation of optimal CWENO-Z or TENO.

Nominal third-/fifth-order spatial accuracy is checked on a smooth periodic FV spatial operator. SSPRK3 has third-order time accuracy. Shock discontinuities, nonlinear limiting and the present wall closure can reduce full-problem convergence order.

## Source-to-module map and corrections

| Archived role | Python module | Migration decision |
| --- | --- | --- |
| INPUT, GEOM | `Config`, `ConicalSolver` in solver.py | Explicit physical configuration; cell averages and separate ghosts |
| INITIALIZE, CALC | state.py | Consistent conserved/primitive state and exact freestream averages |
| BC | solver.py | Explicit inviscid wall and fixed outer state |
| FLUX | flux.py, reconstruction.py | One complete flux vector per face; limited/high-order choices |
| SCI, SMRI, SMTI, SEI | state.py, solver.py | Documented geometric sources and consistent quadrature |
| UPDATE, CONVERGENCE | solver.py | SSPRK stages, all-equation residual gate, finite/positive checks |
| OUTPUT | study.py, render_report.py | Raw fields, JSON metadata, figures and English report |
| VISCUS | viscous.py, fv_mesh.py, steady.py | Corrected local-station viscous model; full 2-D benchmark validation remains open |

Specific issues identified in `all-codes/Second0.for`:

- Unit 5 is opened as an output file before interactive input.
- The second-order tangential left-face flux uses right-face dissipation coefficients.
- The radial and energy second-order branches write `FINVT` instead of their own flux arrays.
- The energy second-order branch uses radial dissipation coefficients.
- Unlimited reconstructed states and viscosity work near uninitialized boundary entries are unsafe.

Accordingly, the comparison is **classical versus higher-order methods on one corrected Python core**, not a bitwise reproduction of the defective Fortran executable. The original default selects uniform flow; the main study deliberately selects its inviscid wall problem.

## Cases, reference and evidence

- Historical-style case: Mach 7.95, cone half-angle 10 degrees, angular outer boundary 22 degrees.
- Independent benchmark: Mach 2.35, cone half-angle 10 degrees, outer boundary 45 degrees.
- Three physical-cell grids: 60, 120 and 240. The original N=240 includes boundary entries, so the discretizations are not identical.
- Independent Taylor–Maccoll shooting: Rankine–Hugoniot shock states, RK4 integration and weak-branch bisection. ODE refinement and wall/shock checks are recorded.

The [NASA cone page](https://www.grc.nasa.gov/www/wind/valid/cone10/cone10.html) has a discrepancy between its printed theoretical pressure/angle and this independently calculated reference. Its listed Wind-US wall pressure (~1.3740) agrees with our Taylor–Maccoll pressure (~1.37393637); the printed theoretical pressure (1.4234) is not silently adopted as truth. The generated report gives all figures and differences.

The full-profile pressure L1 error includes the shock. Wall pressure error is reported separately. Shock angle is estimated from the pressure-gradient peak; that estimate carries grid and interpolation sensitivity.

Run metadata records configurations, residuals, iterations, timing, positive-state minima, geometric flux/source balance and source hashes. Classical fluxes are compared at the same first-order reconstruction; reconstruction changes hold HLLC fixed. CPU timing is host-specific and runs are sequential.

## Limits

`stability_check.py` continues any common-CFL run that misses its residual gate at CFL=0.15, from its saved cell averages. It preserves the original run and records continuation iterations/timing separately. The six-property figures use a converged continuation where available and identify its changed CFL in the legend. The original pressure-only, cost/error and convergence figures retain the common-CFL experiment.

Figures cover pressure, density, local Mach number, ideal-gas temperature, radial velocity and tangential velocity, with shock-detail panels. PNG previews and PDF figures accompany seven-column CSV profiles for every method/grid and the reference. Temperature is derived from the ideal-gas equation of state; these are numerical predictions rather than measurements.

The original comparison is a steady inviscid cone benchmark; the added viscous model and optional characteristic reconstruction are described above. This is not a validated real-gas hypersonic model or a general entropy-stable/positivity-proven scheme. The original CWENO study uses component-wise conserved reconstruction. The reflected wall closure is not a complete fifth-order boundary treatment. HLLC face fallback counts are not separately instrumented; reconstruction-limiting and step-rejection counts are logged. The report distinguishes actual convergence from any exhausted iteration budget.
