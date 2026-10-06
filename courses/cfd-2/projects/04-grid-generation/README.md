# Grid generation

[LAST GRID.FOR](my-grid/LAST%20GRID.FOR) defines airfoil and outer-circle boundaries, sets thickness `T=0.12D0`, and writes grid and error files. The archive also preserves [GRID generation.FOR](GRID%20generation.FOR) from the personal Roohi-code/grid folder.

[Historical GRID1.plt and error.plt](../../results/04-grid-generation) are original outputs. They have not been regenerated. Read [Hejranfar chapter 4](../../lectures/hejranfar/chapter-04.pdf) alongside these programs and examine grid orthogonality, stretching and cell quality before using a mesh in a flow solver.

## Role in the course

Grid generation is a supporting component associated with the Project 3/grid material. Its separate directory is an organizational choice, not evidence of a fourth numbered assignment.

Both programs build a NACA 0012 boundary and a surrounding circular farfield. Linear interpolation supplies an initial mesh, and an elliptic smoothing iteration adjusts interior coordinates. `LAST GRID.FOR` also has an orthogonality-control stage.

| Source | Structured dimensions | Recorded execution |
| --- | --- | --- |
| `GRID generation.FOR` | 51 by 51 points | Completed with exit code 0 |
| `my-grid/LAST GRID.FOR` | 25 by 25 points | Completed with relaxation input 1.0 and exit code 0 |

These are actual Linux diagnostic runs using the preserved sources with no numerical changes. The new output-file sizes are recorded in [Fortran execution evidence](../../docs/validation/fortran-validation.json). Completion confirms that the programs run; full grid-quality acceptance still requires checking generated-cell areas, Jacobians, orthogonality and stretching. The historical grids in `results/` remain identified as historical files.
