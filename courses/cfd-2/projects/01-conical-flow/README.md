# Project 1 - conical flow

[Open the illustrated case page: pressure reconstruction and angular profile →](../../case-studies/conical-flow/README.md)

The source declares `PROGRAM CONICALFLOW` and radial/tangential velocity, density, pressure, temperature and energy fields. `Second0.for` offers interactive choices labeled Roe, Van Leer and H-L-R. The boundary switch labels uniform, inviscid and viscous configurations.

- [Second0.for](all-codes/Second0.for): interactive flux choice; source default `BCTYPE=1`, `SECOND=1`.
- [VAN LEER.for](all-codes/VAN%20LEER.for): defaults to `METHOD=1` (Roe), `BCTYPE=3` (viscous); the filename does not establish the active method.
- [uniform.for](all-codes/uniform.for) and [CONTOUR.for](all-codes/CONTOUR.for): additional archived variants/helper.
- [HLR final directory](hlr-final): preserved alternative source versions.
- [Original project document](../../reports/01-conical-flow/Project.doc).

Parameters, grids and output filenames are embedded in the sources. Read `Subroutine Input` and the boundary/flux routines before selecting a configuration. The original cleanup calls must be removed in a working copy; see [runtime instructions](../../docs/running-legacy-codes.md).

## Physical and numerical interpretation

This is the first numbered assignment. Its five flow-source variants store an angular grid and radial/tangential components of velocity. The `GEOM` routine supplies angular positions, `INITIALIZE` constructs an initial state, `BC` applies boundary values, and `FLUX`/`VISCUS` assemble inviscid/viscous contributions before the update and convergence routines.

One selected source, `all-codes/Second0.for`, uses 240 angular entries, a 10-degree cone angle and `XMIN=7.95`. Its default `BCTYPE=1` is a uniform-flow check rather than the wall case; `BCTYPE=2` and `3` select the inviscid and viscous wall branches. These are source settings, not an independently verified physical solution.

The historical numerical study compares flux schemes and boundary treatments using density, pressure, temperature, velocity and residual outputs. The separate Python modernization below now performs an inviscid comparison with an independently checked Taylor–Maccoll solution under matching assumptions.

## Executed checks

All six files compiled in the Linux diagnostic run. The five flow variants failed runtime checks; the contour postprocessor was not run because it needs flow outputs. The interactive `Second0.for` reopens Fortran unit 5 as `Output.plt` before `read*`, so the method selection reads the file rather than standard input and encounters EOF. Other variants triggered floating-point checks in initialization or viscosity.

These findings are recorded in the [validation report](../../docs/validation/README.md). They prevent claiming a reproduced, validated cone-flow solution.

## 2026 runnable edition

[Open the complete code in Google Colab](https://colab.research.google.com/github/Ehsan-Roohi/Ehsan-Roohi/blob/main/courses/cfd-2/modern-python/conical/Conical_Flow_Modern_2026.ipynb) · [Published results](../../modern-python/conical/RESULTS.md) · [Source and raw-result ZIP](../../downloads/conical-modern-2026.zip)

The single notebook includes 11 selectable face fluxes, classical and high-order finite-volume reconstruction, the local-station viscous model and an optional modal DG companion. It embeds the source, so students can run and modify the algorithms themselves. The DG companion is independently verified on smooth problems, but its degree-one/two cone shock runs are still exploratory; their failures are documented.

## Modern Python implementation and measured results

[Source, equations and migration guide](../../modern-python/conical/README.md) · [Measured comparison report](../../modern-python/conical/results/report.html) · [Verification evidence](../../modern-python/conical/results/verification.json)

The new implementation migrates the inviscid angular model into separate state, flux, reconstruction, solver and reference modules. Roe and Van Leer provide classical baselines on a corrected common core. HLLC, MUSCL, CWENO3 and CWENO5 extend the comparison, with consistent geometric-source quadrature and recorded residual histories. No packaged CFD solver is used.

The report covers Mach 2.35 and the historical-style Mach 7.95 case, three grids and six flux/reconstruction combinations. It gives wall-pressure errors, shock-angle estimates, full-profile errors and computation time, and explicitly marks any run that exhausts its convergence budget. Original Fortran defects are documented in the migration guide.

The [AUSM/viscous extension report](../../modern-python/conical/results-v2/report.html) adds five AUSM-family fluxes and an entropy-based comparison to both inviscid and viscous calculations. The [viscous migration guide](../../modern-python/conical/VISCOUS.md) derives the local radial-station model, audits the historical wall/source routines and explains the independent verification tests. The model now computes no-slip velocity, adiabatic/isothermal wall conditions, skin friction and heat flux. It is not a validated full axisymmetric cone solver or a certification of the archived viscous executable.
