# Code and project guide

| Project | Physical/numerical focus | Guide |
| --- | --- | --- |
| 1 | Conical compressible flow, Roe/Van Leer/H-L-R flux branches | [Conical flow](01-conical-flow/README.md) |
| 2 | Two-dimensional compressible airfoil, triangular connectivity and multistage updates | [Compressible airfoil](02-compressible-airfoil/README.md) |
| 3 | Viscous airfoil, mixed three-/four-node connectivity and Reynolds-number settings | [Viscous airfoil](03-viscous-airfoil/README.md) |
| Grid | Airfoil boundary geometry and structured grid generation | [Grid generation](04-grid-generation/README.md) |

These descriptions follow the actual source code, not just filenames. For example, `VAN LEER.for` defaults to `METHOD=1`, whose comment identifies Roe; select the intended method explicitly. Keep each case's dimensions, connectivity, boundary files and source configuration together.
