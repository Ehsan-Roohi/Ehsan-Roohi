# 02 · Compressible airfoil

**CFD II · Dr. Mazaheri · Ehsan Roohi**

Where does subsonic freestream flow become locally supersonic? This assignment studies compressible flow around an airfoil on an unstructured triangular mesh, with acceleration, compression and surface pressure visible in the archived fields.

![Historical airfoil Mach contour with sonic isoline and surface pressure-coefficient curves](figure.png)

| Preserved case | Mesh | Source Mach setting |
| --- | --- | --- |
| 85-medium | 1,298 nodes · 2,526 triangles | 0.85 |

## Read the figure

The colored field displays the output's stored `MACH` variable. Warm colors identify higher Mach number; the white isoline marks **M = 1**. The preserved field spans approximately **0.571–1.436**, showing locally supersonic regions even though the corresponding source sets the freestream Mach parameter to 0.85.

The right-hand panel plots the stored `CPP` pressure coefficient at upper- and lower-surface boundary nodes. Its vertical axis is inverted, as is common in airfoil pressure plots: stronger suction appears higher on the chart. Differences between the surface curves connect the flow field to aerodynamic loading, although no validated lift or drag coefficient is claimed here.

Contours use the **original 2,526 triangles**. The airfoil interior remains empty, and surface nodes are identified from boundary edges rather than from a rectangular interpolation grid.

## What the code does

The Fortran solver evolves the compressible conservation equations for mass, two momentum components and total energy. Roe-type fluxes connect neighboring triangular cells; multistage updates march toward a steady state. Thermodynamic recovery reconstructs pressure, sound speed and Mach number.

Useful comparisons include mesh refinement, discontinuity resolution, wall and farfield treatment, and conservation-residual reduction. The archive also contains finer meshes and other historical configurations, described in the project guide.

## Explore this project

[Code and mesh variants](../../projects/02-compressible-airfoil/README.md) · [85-medium source and inputs](../../projects/02-compressible-airfoil/cases/85-medium) · [Original report](../../reports/02-compressible-airfoil/Report%202.doc) · [Preserved field](../../results/02-airfoil-85-medium/output.plt)

**Figure status:** historical output, newly visualized. The source Mach setting is documented; the original output's convergence and aerodynamic accuracy have not been independently established. See the [execution checks](../../docs/validation/README.md) and [figure provenance](../figure-provenance.json).

---

[← Conical flow](../conical-flow/README.md) · [All three cases](../README.md) · [Next: viscous airfoil →](../viscous-airfoil/README.md)
