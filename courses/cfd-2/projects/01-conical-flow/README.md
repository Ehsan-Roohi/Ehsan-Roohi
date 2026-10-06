# Project 1 - conical flow

The source declares `PROGRAM CONICALFLOW` and radial/tangential velocity, density, pressure, temperature and energy fields. `Second0.for` offers interactive choices labeled Roe, Van Leer and H-L-R. The boundary switch labels uniform, inviscid and viscous configurations.

- [Second0.for](all-codes/Second0.for): interactive flux choice; source default `BCTYPE=1`, `SECOND=1`.
- [VAN LEER.for](all-codes/VAN%20LEER.for): defaults to `METHOD=1` (Roe), `BCTYPE=3` (viscous); the filename does not establish the active method.
- [uniform.for](all-codes/uniform.for) and [CONTOUR.for](all-codes/CONTOUR.for): additional archived variants/helper.
- [HLR final directory](hlr-final): preserved alternative source versions.
- [Original project document](../../reports/01-conical-flow/Project.doc).

Parameters, grids and output filenames are embedded in the sources. Read `Subroutine Input` and the boundary/flux routines before selecting a configuration. The original cleanup calls must be removed in a working copy; see [runtime instructions](../../docs/running-legacy-codes.md).
