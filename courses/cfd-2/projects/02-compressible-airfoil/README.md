# Project 2 - compressible airfoil

[Open the illustrated case page: Mach contour and surface pressure →](../../case-studies/compressible-airfoil/README.md)

The archived `AIRFOIL` programs evolve density, two velocity components, pressure and energy, with Roe-related flux routines and multistage update coefficients. They are compressible-flow codes, whereas the separately reviewed SIMPLE solver is incompressible.

[Root source variants](variants) are preserved separately from five case directories: [85-medium](cases/85-medium), [85-high](cases/85-high), [very-high](cases/very-high), [exact](cases/exact), and [uniform](cases/uniform).

The directory names are historical labels, not a guarantee of Mach number, resolution or exactness. Inspect each source's `Input` routine. For example the root `airfoilfinal.for` sets `XMIN=0.85` and an angle of 1 degree; root `airfoil4.for` sets `XMIN=0.80` and 1.25 degrees. Do not interchange these settings based on names alone.

The selected 85-medium case opens `Node.DAT`, `connectivity.DAT`, `BC.DAT` and `BE.DAT`. It and several other variants import `MSIMSL` / `MSFLIB`; a modern compiler needs a deliberate compatibility review. Some root variants reference boundary files that were not co-located with them, so root variants are not presented as ready-to-run cases.

[Project brief](../../reports/02-compressible-airfoil/Project%20II.doc) · [Report](../../reports/02-compressible-airfoil/Report%202.doc) · [Mach 2 appendix](../../reports/02-compressible-airfoil/Appendix%20M%3D2%20Results.doc) · [Selected old results](../../results/README.md).

## What this assignment is studying

This is the second numbered course project. The main airfoil codes discretize compressible conservation equations for mass, two momentum components and total energy. The neighbor and face routines connect triangular cells; Roe-type characteristic fluxes transfer conserved quantities across faces. Multistage updates march the state toward a steady solution, and thermodynamic recovery reconstructs pressure, temperature, sound speed and Mach number.

Its key numerical questions include shock resolution, sensitivity to the mesh and pseudo-time step, farfield/wall treatment, and whether all conservation residuals decrease. A historical case directory called `exact` is not proof of an analytical exact solution.

## Case sizes and settings read from the sources

| Directory | Nodes | Triangular cells | Source Mach parameter |
| --- | ---: | ---: | ---: |
| [85-medium](cases/85-medium) | 1298 | 2526 | 0.85 |
| [85-high](cases/85-high) | 4012 | 7934 | 0.85 |
| [very-high](cases/very-high) | 5284 | 10488 | 0.85 |
| [exact](cases/exact) | 409 | 776 | 1.2 |
| [uniform](cases/uniform) | 258 | 464 | 1.2 |

All five input meshes passed the recorded connectivity, signed-area and source-dimension checks. That establishes basic input consistency, not flow accuracy.

## Additional source material

The seven root variants are not seven extra course assignments. Five are airfoil iterations/configurations. `unstructured0.for` is a temperature/conduction numerical prototype, with `T`, conductivity and thermal boundary data; it is not another full airfoil-flow solver. The free-form `NUM1,CFD2,PRO2,T,uniform.f90` declares `program AliGHAFFARIYAN`; its own identifier is retained and it is not asserted to be authored by Roohi.

The root airfoil variants use dimensions such as 481/892 or 2177/4124, which do not match the five preserved case meshes above. No unrelated mesh was substituted merely to obtain an executable run.

## Executed checks

All 13 files compiled after the documented working-copy adaptations where needed. The six source files in the five co-located mesh cases were executed with bounds and floating-point checks; all stopped on invalid arithmetic or uninitialized indices. The seven root variants were built but not executed because their matching inputs were not co-located in the published selection.

See the [validation report](../../docs/validation/README.md) for the precise failure locations. New validated aerodynamic results are not claimed.
