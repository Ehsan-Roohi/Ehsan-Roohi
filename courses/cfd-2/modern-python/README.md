# Student-written Python CFD modernization

The first implemented module is [conical flow](conical/README.md), with independently reimplemented classical fluxes and higher-order finite-volume reconstruction. [Measured results](conical/results/report.html) compare both Mach 7.95 and Mach 2.35 cones.

The [AUSM and viscous extension](conical/results-v2/report.html) adds AUSM, AUSM+, AUSM+-up, AUSM+-up2, SLAU2 and an entropy-conservative/LF flux to the common core. It includes a local viscous cone-station model, wall heat/shear comparisons, characteristic reconstruction and higher-order interior diffusion. [Equations and limitations](conical/VISCOUS.md) distinguish verified components from full viscous-cone validation.

The historical Fortran archive is preserved. Euler-airfoil and viscous-airfoil Python migrations have not yet been implemented. No external CFD solver performs the computations in this module.

[Run the all-methods Colab notebook](https://colab.research.google.com/github/Ehsan-Roohi/Ehsan-Roohi/blob/main/courses/cfd-2/modern-python/conical/Conical_Flow_Modern_2026.ipynb) · [Published results](conical/RESULTS.md) · [Optional DG companion](conical/DG.md) · [Download code and raw data](../downloads/conical-modern-2026.zip).

## Expanded flux/reconstruction and shock-tube comparison

[Expanded cone algorithms](conical/ADDITIONAL_METHODS.md) · [Separate shock-tube Colab and report](shock-tube/README.md)

15 shared Euler fluxes, 20 FV reconstructions and standalone JST are available. The cone notebook retains the original study and adds 46 cases; the separate shock-tube notebook covers 592 uniform/DG runs plus [96 uniform-grid comparisons](shock-tube/UNIFORM_GRID_RESULTS.md). The full 300-combination Sod matrix remains available. Failed cases remain visible.
