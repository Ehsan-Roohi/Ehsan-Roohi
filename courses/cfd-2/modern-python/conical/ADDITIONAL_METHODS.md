# Additional finite-volume algorithms and shock-tube cross-checks

This extension is independently implemented from mathematical formulations. Upstream repositories are references and independent comparison candidates; their code is not republished as our authorship.

## Inventory and provenance

| Component | Added choices | Reference and scope |
|---|---|---|
| Pointwise Euler flux | Exact Godunov, Roe without entropy smoothing, Steger–Warming, semidiscrete global LF | [Python shock tube](https://github.com/chairmanmao256/Python-shock-tube), [Toro-style C++ examples](https://github.com/cangyu/Riemann-Solvers), [comparison framework](https://github.com/fhermet/euler-1d-solver) |
| MUSCL slopes | Minmod, van Leer, superbee, van Albada, alongside existing MC | Conservative cell-average linear polynomials |
| FV face reconstruction | ENO2; characteristic WENO-JS/Z3/5; characteristic WENO-JS7 | Candidate polynomials and optimal face weights derived from exact FV moments; no nodal FD substitution |
| Targeted reconstruction | TENO3/5/7 | TENO5 uses tau5 and binary stencil cutoff. TENO3/7 use the inverse-smoothness cutoff strategy visible in [CAELUM](https://github.com/navasmontilla/CAELUM/blob/main/lib/reconst.c). Fixed C_T=1e-5, exponent 6; these are FV variants, not adaptive TENO-A/LAD |
| Interface sharpening | MUSCL–THINC–BVD | Independent component-wise one-stage BVD variant with beta=1.6. It is not the published multistage P4-THINC-BVD algorithm. [BVD paper](https://arxiv.org/abs/1602.00814), [MATLAB reference collection](https://github.com/wme7/ApproximateRiemannSolvers) |
| Central stencil | JST | Pressure-sensed second/fourth artificial differences, k2=0.5/k4=0.02, SSPRK3. Standalone uniform-grid FV; not a composable Riemann solver or viscous/DG discretization. [SU2 documentation](https://su2code.github.io/docs_v7/Convective-Schemes/) |
| DG shock-tube counterpart | Modal degrees 0/1/2, all 15 shared fluxes | The cone DG core's weak form is reused without its geometric source or cone boundary. Characteristic TVB and sampled positivity scaling. [Quail Sod example](https://github.com/IhmeGroup/quail/tree/main/examples/euler/1D/sod_problem) |

The existing `hlle` key uses Davis min/max acoustic bounds; it is HLL-family flux, not Roe-averaged Einfeldt bounds. HLL from the reviewed Python repositories is not counted as a new method. The van Leer slope limiter and the van Leer flux are separate components. Counts refer to 15 pointwise fluxes, 20 FV reconstructions and one standalone JST stencil, not 300 distinct flux models.

The NASA [AUSM shock-tube code](https://github.com/nasa/shocktube) supports AUSM/AUSM+ and limited face reconstruction; those families were already present. ADERDG, arbitrary 2-D unstructured DG/FR, reacting flow and physical-viscosity shock tubes are not implemented. Classical Lax–Wendroff is a coupled space/time scheme and is not relabelled as a pointwise flux; this extension concentrates on shared FV flux/reconstruction algorithms and SSPRK3.

## Cone source and boundary treatment

WENO/TENO reconstruct faces; separate conservative CWENO polynomials supply geometric-source quadrature and center diagnostics. WENO3/TENO3 use CWENO3 source quadrature; fifth/seventh-order faces use CWENO5 source quadrature. Therefore seventh-order face reconstruction does not imply seventh-order global cone accuracy. Boundary ghost extension and limiting further reduce problem order. On stretched viscous meshes the added pointwise fluxes use the existing supported primitive-minmod/G2 operator; uniform-grid WENO/TENO face variants are not silently applied to the nonuniform mesh.

The THINC variant affects face states only; cone source quadrature retains MUSCL polynomials. It is experimental, and its failure to reach a steady cone residual is retained. Roe without entropy correction is also an educational diagnostic, not a universal default.

## Executed cone study

All new cases use N=60, the same physical parameters as the original comparison, and configuration-matched published states as initial guesses. Initializing from a different numerical method does not prescribe the converged answer. Each new method recomputes its residual and evolves its own state. Seeded iteration counts must not be ranked as fresh-start speed measurements.

46 new cone runs were executed; 31 met the recorded steady acceptance gate. Raw fields and residual histories are in [results-additional](results-additional/catalog.json).

| Physics / Mach | Flux | Reconstruction | Gate | Final residual |
|---|---|---|---|---|
| adiabatic / 7.95 | global-lf | primitive-minmod | Reached | 1.143e-08 |
| adiabatic / 7.95 | godunov | primitive-minmod | Reached | 1.705e-12 |
| adiabatic / 7.95 | roe-nc | primitive-minmod | Reached | 7.718e-11 |
| adiabatic / 7.95 | steger-warming | primitive-minmod | Reached | 4.867e-09 |
| inviscid / 2.35 | global-lf | cweno3 | Reached | 1.786e-08 |
| inviscid / 2.35 | godunov | cweno3 | Reached | 1.766e-08 |
| inviscid / 2.35 | hllc | eno2 | Not reached | 4.377e-03 |
| inviscid / 2.35 | hllc | muscl-minmod | Not reached | 3.343e-03 |
| inviscid / 2.35 | hllc | muscl-superbee | Not reached | 1.070e-02 |
| inviscid / 2.35 | hllc | muscl-thinc-bvd | Not reached | 5.011e-01 |
| inviscid / 2.35 | hllc | muscl-vanalbada | Reached | 1.870e-08 |
| inviscid / 2.35 | hllc | muscl-vanleer | Reached | 1.507e-08 |
| inviscid / 2.35 | hllc | teno3 | Not reached | 1.853e-03 |
| inviscid / 2.35 | hllc | teno5 | Reached | 1.810e-08 |
| inviscid / 2.35 | hllc | teno7 | Reached | 1.885e-08 |
| inviscid / 2.35 | hllc | weno3-js | Reached | 1.851e-08 |
| inviscid / 2.35 | hllc | weno3-z | Reached | 1.806e-08 |
| inviscid / 2.35 | hllc | weno5-js | Reached | 1.791e-08 |
| inviscid / 2.35 | hllc | weno5-z | Reached | 1.781e-08 |
| inviscid / 2.35 | hllc | weno7-js | Reached | 1.839e-08 |
| inviscid / 2.35 | jst | first | Reached | 1.777e-08 |
| inviscid / 2.35 | roe-nc | cweno3 | Reached | 1.992e-08 |
| inviscid / 2.35 | steger-warming | cweno3 | Reached | 1.924e-08 |
| inviscid / 7.95 | global-lf | cweno3 | Reached | 1.805e-08 |
| inviscid / 7.95 | godunov | cweno3 | Reached | 1.904e-08 |
| inviscid / 7.95 | hllc | eno2 | Not reached | 8.533e-01 |
| inviscid / 7.95 | hllc | muscl-minmod | Not reached | 5.250e-01 |
| inviscid / 7.95 | hllc | muscl-superbee | Not reached | 1.166e+00 |
| inviscid / 7.95 | hllc | muscl-thinc-bvd | Not reached | 1.294e+00 |
| inviscid / 7.95 | hllc | muscl-vanalbada | Not reached | 5.826e-01 |
| inviscid / 7.95 | hllc | muscl-vanleer | Not reached | 6.533e-01 |
| inviscid / 7.95 | hllc | teno3 | Not reached | 1.686e-01 |
| inviscid / 7.95 | hllc | teno5 | Not reached | 1.188e-01 |
| inviscid / 7.95 | hllc | teno7 | Not reached | 1.799e-01 |
| inviscid / 7.95 | hllc | weno3-js | Reached | 1.752e-08 |
| inviscid / 7.95 | hllc | weno3-z | Reached | 1.716e-08 |
| inviscid / 7.95 | hllc | weno5-js | Reached | 1.773e-08 |
| inviscid / 7.95 | hllc | weno5-z | Reached | 1.841e-08 |
| inviscid / 7.95 | hllc | weno7-js | Not reached | 1.555e-01 |
| inviscid / 7.95 | jst | first | Reached | 1.858e-08 |
| inviscid / 7.95 | roe-nc | cweno3 | Reached | 1.751e-08 |
| inviscid / 7.95 | steger-warming | cweno3 | Reached | 1.977e-08 |
| isothermal / 7.95 | global-lf | primitive-minmod | Reached | 6.569e-11 |
| isothermal / 7.95 | godunov | primitive-minmod | Reached | 1.236e-10 |
| isothermal / 7.95 | roe-nc | primitive-minmod | Reached | 2.937e-09 |
| isothermal / 7.95 | steger-warming | primitive-minmod | Reached | 1.052e-08 |

## Independent verification

[Executed component verification](additional-verification.json) checks all 15 fluxes on equal states, the Godunov sampler against the independent original Sod reference, a true vacuum sample, and measured smooth FV operator orders for all 20 reconstructions. WENO3 nonlinear critical-point order reduction and shock/wall limitations are shown rather than inferred away. Tests validate recorded properties, not all possible flows.

[Stored-field integrity audit](results-additional/field-integrity.json) independently checks all 46 extension fields, positive conservative means, finite histories, recorded residual gates and the explicit drift check for accepted viscous states.

Run locally from this folder:

```bash
python verify_additional.py
python additional_study.py
python shock_study.py
```

[Separate shock-tube notebook and report](../shock-tube/README.md) cover the full 15-by-20 Sod matrix, standalone JST, all fluxes with DG degrees 0/1/2, targeted three-grid studies and additional wave tests. Physical viscosity is not included in this Euler benchmark.
