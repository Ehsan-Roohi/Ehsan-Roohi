# Project 1 - conical flow

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

The numerical study is intended to compare flux schemes and boundary treatments using density, pressure, temperature, velocity and residual outputs. A suitable inviscid conical-flow validation would compare with an appropriate Taylor-Maccoll solution under matching assumptions; that comparison has not been completed here.

## Executed checks

All six files compiled in the Linux diagnostic run. The five flow variants failed runtime checks; the contour postprocessor was not run because it needs flow outputs. The interactive `Second0.for` reopens Fortran unit 5 as `Output.plt` before `read*`, so the method selection reads the file rather than standard input and encounters EOF. Other variants triggered floating-point checks in initialization or viscosity.

These findings are recorded in the [validation report](../../docs/validation/README.md). They prevent claiming a reproduced, validated cone-flow solution.
