# 03 · Viscous airfoil

**CFD II · Dr. Mazaheri · Ehsan Roohi**

How should a mesh resolve the strong gradients near an airfoil wall? This assignment adds viscous and thermal fluxes to the compressible airfoil formulation and uses a mixed mesh with quadrilateral layers near the surface and triangles farther away.

![Actual viscous-airfoil mesh colored by logarithmic cell area, with a leading-edge detail](figure.png)

| Preserved mesh | Source Reynolds setting | Source Mach / angle setting |
| --- | --- | --- |
| 838 nodes · 232 quads + 1,142 triangles | 1,000 | 1.2 / 0° |

## Read the figure

Colors represent **cell area**, calculated directly from the original node coordinates and connectivity. The logarithmic scale makes the large differences between small wall-adjacent cells and larger outer cells visible. Dark colors indicate smaller cells; brighter colors indicate larger cells.

The close-up shows the quadrilateral layers around the leading edge and their transition into triangles. This geometry matters because viscous stresses and heat flux depend on gradients: wall spacing, cell shape and the consistency of reconstruction all affect the calculation.

This is a **mesh visualization**, not a velocity, pressure or temperature contour. No preserved flow-field output was located in the selected personal Roohi-code viscous directory. The actual project mesh provides a traceable illustration without substituting another person's solution or another solver's airfoil result.

## What the code does

The Fortran solver combines inviscid transport with viscous momentum and thermal flux contributions. Separate connectivity files describe quadrilateral and triangular cells. Nodal reconstruction and wall treatment supply information for the viscous terms; the preserved source sets the Reynolds parameter to 1,000 and the Mach parameter to 1.2.

Questions to investigate include sensitivity to wall-normal resolution, gradient reconstruction on mixed cells, boundary conditions, and eventual comparisons of surface pressure and skin friction with matching reference data.

## Explore this project

[Code and numerical guide](../../projects/03-viscous-airfoil/README.md) · [Viscous source and mesh inputs](../../projects/03-viscous-airfoil/roohi-code/Viscous) · [Original report](../../reports/03-viscous-airfoil/Report%203.doc)

**Figure status:** verified mesh geometry, newly visualized. Basic mesh checks pass, while the flow solver still has documented runtime and cell-centroid defects. See the [execution checks](../../docs/validation/README.md) and [figure provenance](../figure-provenance.json).

---

[← Compressible airfoil](../compressible-airfoil/README.md) · [All three cases](../README.md)
