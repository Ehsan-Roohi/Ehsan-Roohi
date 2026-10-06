# Three CFD II case studies

**Dr. Mazaheri · Ehsan Roohi**

Three numbered course projects, each with a dedicated visual page, an explanation of the numerical problem, and links to the original code and report.

## 01 · Conical flow

[![Conical-flow pressure visualization](conical-flow/figure.png)](conical-flow/README.md)

An archived angular pressure profile reconstructed around a cone. [Open the conical-flow page →](conical-flow/README.md)

## 02 · Compressible airfoil

[![Compressible-airfoil Mach visualization](compressible-airfoil/figure.png)](compressible-airfoil/README.md)

Historical Mach contours, the sonic isoline and surface pressure coefficients on a triangular mesh. [Open the compressible-airfoil page →](compressible-airfoil/README.md)

## 03 · Viscous airfoil

[![Viscous-airfoil mixed-mesh visualization](viscous-airfoil/figure.png)](viscous-airfoil/README.md)

The actual mixed mesh, colored by cell area, with a close-up of the leading-edge wall layers. [Open the viscous-airfoil page →](viscous-airfoil/README.md)

## Figure provenance and reproduction

The cone and compressible-airfoil pictures visualize preserved historical output. The viscous-airfoil picture visualizes mesh geometry. None is presented as a newly validated flow solution.

Install NumPy and Matplotlib in a suitable Python environment, then run:

```text
python tools/render_case_figures.py
```

The script checks dimensions, finite numerical values, connectivity bounds and nonzero cell areas as applicable. [Figure provenance](figure-provenance.json) records input hashes, data ranges and the transformations used. [Solver validation](../docs/validation/README.md) describes the separate execution and accuracy checks.

[← Course home](../README.md) · [Lecture notes](../lectures/README.md) · [Detailed code guides](../projects/README.md)
