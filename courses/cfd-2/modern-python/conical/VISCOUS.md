# AUSM and local viscous-cone extension

The 2026 Python extension adds the same independently implemented inviscid face-flux families to both Euler and viscous calculations. No external CFD solver runs the cases. This is a modernization of the archived course model, with an explicit distinction between verification of numerical components and validation of the physical cone model.

[Executed figures and run tables](results-v2/report.html) · [Verification](results-v2/verification.json) · [Exact shock-tube benchmark](results-v2/sod-verification.json)

## Algorithms and dates

| CLI name | Algorithm | Role |
| --- | --- | --- |
| `ausm` | Liou–Steffen AUSM, 1993 | Original mass/advection and pressure splitting |
| `ausm-plus` | Liou AUSM+, 1996 | Common critical sound speed, fourth/fifth-degree subsonic split polynomials |
| `ausm-up` | Liou AUSM+-up, 2006 | Pressure correction in mass flux and velocity correction in pressure flux |
| `ausm-up2` | Kitamura–Shima AUSM+-up2, 2013 | Modified pressure flux with velocity-magnitude-dependent dissipation |
| `slau2` | Kitamura–Shima SLAU2, 2013 | Low-dissipation mass flux and shock-oriented pressure flux |
| `ec-lf` | Chandrashekar EC flux, 2013, plus LF dissipation | Logarithmic means; independently checked face entropy identity/inequality |

The remaining Roe, Van Leer, HLLC, HLLE and Rusanov choices remain available. These algorithms are established methods newly added to the course code, not methods invented in 2026. The [2013 pressure-flux paper](https://doi.org/10.1016/j.jcp.2013.02.046) specifically motivates AUSM+-up2 and SLAU2 for hypersonic heating. All-speed HLLC extensions such as [Kitamura–Shima HLLCL, 2020](https://doi.org/10.1002/fld.4782) are relevant further reading, but are not claimed as implemented here.

AUSM+ variants use total-enthalpy critical sound speed and canonical signed normal-velocity denominators. AUSM+-up uses Kp = 0.25, Ku = 0.75, sigma = 1 and beta = 1/8. The reference cutoff is the supplied upstream Mach, capped at 1 by the formula; the present supersonic cases therefore use fa = 1. A low-Mach multidimensional benchmark has not been run. SLAU2 includes the standard counterflow and pressure terms, with no optional reduced-dissipation sensor. The pinned [SU2 v8.5.0 mathematical implementation](https://github.com/su2code/SU2/blob/v8.5.0/SU2_CFD/src/numerics/flow/convection/ausm_slau.cpp) was used for formula cross-checking; its source was not copied and its CFD solver was not run.

The [Chandrashekar flux](https://arxiv.org/abs/1209.4994) passes the two-state Tadmor identity. Adding LF dissipation gives a dissipative face entropy inequality. This does not prove global entropy stability for nonlinear reconstruction, conical sources, boundaries and SSPRK time stepping. Characteristic CWENO3/5 uses a local angular Euler eigenbasis; the optional fourth-order viscous fit uses control-volume averages, not nodal finite differences.

## Physical model and nondimensionalization

The default historical-style station uses Mach 7.95, cone angle 10°, outer angle 22°, gamma = 1.4, Re = 420,000 and Pr = 0.72. Upstream density and velocity equal 1. Upstream pressure is 1/(gamma Mach²). Tinf = 56.857153 K follows the archived stagnation temperature 775.56 K at that Mach. The ideal-gas temperature ratio is T = gamma Mach² p/rho.

This is a **fixed radial station** approximation. Primitive radial derivatives are suppressed; the spherical velocity-gradient and stress geometry remain. Physical stresses scale as mu/r, so radial geometric divergence terms remain even when primitive radial derivatives vanish. Re is local to the selected station. A real viscous cone loses exact radial similarity because Re changes along the cone. This solver cannot predict downstream boundary-layer development; it is not a full axisymmetric Navier–Stokes validation or a turbulence model.

Define mu = muhat/Re and differentiate with respect to theta. The Sutherland ratio is:

```text
muhat(T) = T^(3/2) (Tinf + 110.4) / (T Tinf + 110.4)

tau_rr = -(2/3) mu (vt' + 2 vr + vt cot(theta))
tau_tt =  (2/3) mu (2 vt' + vr - vt cot(theta))
tau_pp =  (2/3) mu (-vt' + vr + 2 vt cot(theta))
tau_rt = mu (vr' - vt)
qtheta = -mu T' / [Pr (gamma - 1) Mach²]

Fv = [0, tau_rt, tau_tt, tau_rt vr + tau_tt vt - qtheta]
Sv = [0,
      tau_rr - tau_tt - tau_pp + tau_rt cot(theta),
      2 tau_rt + (tau_tt - tau_pp) cot(theta),
      tau_rr vr + tau_rt vt + cot(theta)(tau_rt vr + tau_tt vt - qtheta)]

dU/dpseudo_time = -dF/dtheta + S + dFv/dtheta + Sv
```

The trace tau_rr + tau_tt + tau_pp is zero. A uniform Cartesian velocity rotated into the spherical basis produces zero stress; this is a verification test. Wall Cf = 2 tau_rt,w because the stress scale is rhoinf Uinf². Reported qwall is normalized by rhoinf Uinf³ and positive **into the fluid**; a cooled wall normally gives negative qwall.

## Wall and mesh

The mesh faces use an exponential mapping with stretch = 3 by default. Ghost widths mirror the first three physical cells at the cone; upstream states are exact cell averages in the rotating basis at the outer boundary.

The convective wall flux is [0, 0, pwall, 0]. Both velocity components obey no slip. Thermal conditions act in the viscous flux and source closure: qwall = 0 for adiabatic walls; Tw/Tinf is prescribed for isothermal walls. The default isothermal comparison uses Tw/Tinf = 1 (about 56.86 K), which is a mathematical cold-wall case, not a proposed laboratory temperature. It deliberately contrasts heat removal with adiabatic heating.

The corrected `fv-wall-quadratic-v2` closure fits a quadratic to the first two **primitive control-volume values**, imposing wall velocity or temperature, or zero wall temperature derivative. It supplies both the wall gradient and first-cell source polynomial. The fit has second-order derivative accuracy when supplied with accurate primitive averages, as independently tested. G4 obtains primitive averages by conservative-polynomial deconvolution and quadrature; G2 converts conservative averages directly to primitive values. That conversion can reduce wall-transport accuracy in a strongly varying density/temperature layer. The metadata's `wall_closure_order = 2` describes the formal fit operator, not a measured global boundary order for the coupled solver. The smooth interior G4 diffusion operator has measured fourth order; the complete shock-containing viscous solver is not claimed to have global fourth or fifth order.

Early exploratory outputs without the operator-version field used a first-order wall derivative. Their metadata originally labeled the wall order as 2; that label was incorrect. They remain diagnostic snapshots and are excluded from the corrected-wall comparison. Original Fortran files and the earlier inviscid results are unchanged.

## What the source audit found

In `projects/01-conical-flow/all-codes/Second0.for`, BCTYPE = 3 reflects both velocity components and copies density/energy, so the actual wall branch is no slip and adiabatic. The TCTI wall-temperature-like input is not applied in that branch. The left viscous source contains `DUTDT2(1)` inside the cell loop rather than a cell-indexed value, and uses right-face `XMUR(I)` in a left-face expression. The boundary routines also overwrite the first two old radial-velocity entries. The Python model derives shared face stresses and geometric sources directly and avoids those source-index and update-array defects. It does not certify the archived executable as correct.

## Execution and acceptance

The explicit solver uses a convective CFL bound and a parabolic viscosity/thermal diffusion bound. Sampled reconstruction states must have positive density and pressure; invalid RK stages are rejected and the time step is reduced.

The optional steady solver evaluates the **same FV residual**. It uses colored finite-difference Jacobians, pseudo-transient regularization, a damped positivity-aware line search and short explicit recoveries when necessary. The Jacobian coloring is independently checked against column-by-column differencing. The residual gate is the maximum of four equation-wise RMS values, below 2e-8. Every accepted Newton solution is additionally advanced for 20 explicit SSPRK steps and checked for drift. Small-step persistence is evidence of a consistent stationary discrete state, not proof of long-time stability.

Minmod is nonsmooth, so a finite-difference Jacobian can disagree with a random directional secant near limiter switches. Adaptive perturbations improve the final residual solve. All acceptance decisions use the actual nonlinear residual, never a linear-solver tolerance.

For the controlled viscous flux comparison, every method receives the same wall-matched coarse seed. The seed was prepared with AUSM+-up2 and is already converged for that flux. Therefore iteration counts and seeded concurrent timings are not fair speed rankings. Grid restarts use a wall-aware conservative interpolant integrated over the fine cells; the interpolant is only an initializer. Failed attempts, exact initial-checkpoint hashes and subsequent continuations/restarts remain available.

The independent checks cover equal-state/orientation flux properties, face entropy, exact Sod shock flow, control-volume conservation, stress/heat scaling, uniform Cartesian stress, constant-property Couette heating, smooth nonuniform reconstruction, interior diffusion and viscous angular operators. Taylor–Maccoll is only the inviscid reference. Full viscous-cone validation requires external physical data or a matching full axisymmetric reference and remains open.

## Reproduce

```bash
python verify_extension.py
python shock_tube.py
python extension_study.py euler --grids 60 120
python extension_study.py characteristic --grids 60
python extension_study.py steady-warm --grids 60
python extension_study.py steady-refine --grids 120 240
python extension_study.py steady-high-order --grids 60 120 240
python continue_steady.py --pattern 'fv2-steady-*-n60-*-primitive-minmod-g2.json'
python extension_report.py
```

The included `warm-adiabatic.npz` and `seed-isothermal.npz` are the common seeds. To generate them from the included explicit precursor, run `steady.py --checkpoint results-v2/pilot-primitive.npz --output results-v2/warm-adiabatic`, then `steady.py --wall isothermal --checkpoint results-v2/warm-adiabatic.npz --relax-steps 5000 --output results-v2/seed-isothermal`. All metadata, fields, figures and sources are English.
