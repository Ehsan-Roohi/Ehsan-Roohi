# 01 · Conical flow

**CFD II · Dr. Mazaheri · Ehsan Roohi**

How does a cone compress a high-speed stream? This assignment studies compressible flow in conical coordinates, using an angular profile to describe the region between the cone surface and the outer flow.

![Conical-flow pressure reconstruction and original angular pressure profile](figure.png)

| Figure data | Geometry used for display | Pressure ratio in the archive |
| --- | --- | --- |
| 238 stored angular samples | 10° cone half-angle | 1.000–4.357 |

## Read the figure

Warm colors show elevated pressure near the cone; cool colors show the return toward the outer-flow pressure. The curve on the right is the original pressure profile, including its abrupt drop over a narrow angular interval. That feature is useful for studying how a numerical flux represents a compression discontinuity; the image alone does not establish shock accuracy.

The archive stores a **one-dimensional angular solution**. The colored view extends each stored pressure value along a ray under a self-similar interpretation, then mirrors the sector about the axis. The displayed radial scale is arbitrary. The 10° half-angle follows the preserved source geometry; the file's `TETA` column is an angular offset above the cone surface. No extra flow values are invented outside the stored angular range.

## What the code does

The Fortran programs update density, radial and tangential velocity, pressure, temperature and energy. Angular fluxes and geometric source terms describe the conical coordinate system. Alternative branches compare Roe, Van Leer and H-L-R fluxes, while boundary switches select uniform-flow, inviscid-wall or viscous-wall configurations.

The educational questions are how the flux scheme changes discontinuity resolution, how the wall condition affects the near-cone state, and whether the solution approaches an appropriate Taylor–Maccoll reference under matching inviscid assumptions.

## Explore this project

[Code and numerical guide](../../projects/01-conical-flow/README.md) · [Original project document](../../reports/01-conical-flow/Project.doc) · [Preserved angular output](../../results/01-conical-flow/ROE%20FINAL.plt)

**Figure status:** historical output, newly visualized. The original run configuration and convergence have not been independently established; the filename `ROE FINAL.plt` does not prove the active flux scheme. See the [execution checks](../../docs/validation/README.md) and [figure provenance](../figure-provenance.json).

---

[All three cases](../README.md) · [Next: compressible airfoil →](../compressible-airfoil/README.md)
