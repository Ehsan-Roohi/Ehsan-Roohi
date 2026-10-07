# Optional high-order DG companion

The finite-volume implementation remains the main modernization of the historical cone project. `dg.py` adds an independent **modal Runge–Kutta discontinuous Galerkin** discretization of exactly the same inviscid angular equations. It is a compact teaching implementation for studying a second high-order direction.

## What is implemented

- Legendre modal polynomials of degree 0, 1 or 2, with the zeroth coefficient equal to the conserved cell average.
- A weak-form operator with at least `degree + 2` Gauss points, diagonal exact polynomial mass matrices, one shared numerical face flux, and the same spherical geometric source.
- All eleven Euler fluxes available in the finite-volume code.
- SSPRK3 with a degree-dependent CFL reduction, characteristic TVB slope limiting, sampled positivity scaling that preserves the mean, and rejected-step recovery.
- A slip wall and a prescribed upstream state at the angular outer boundary.

DG is a discretization, whereas Roe/HLLC/AUSM are numerical face fluxes. They can be combined. This module is not viscous DG and is not a two-dimensional unstructured-mesh solver. Viscous calculations continue to use the finite-volume model in [VISCOUS.md](VISCOUS.md).

## Verification and executed cone runs

[Verification JSON](results-dg/verification.json) records:

| Check | Measured result |
|---|---:|
| Degree-zero reduction to the first-order FV interior/wall operator | 0 absolute difference |
| Integrated conserved-mean flux/source balance | 1.11 × 10⁻¹⁶ |
| Limiter preservation of all conserved means | 0 absolute difference |
| Last observed order, smooth periodic advection, degree 1 | 1.994 |
| Last observed order, smooth periodic advection, degree 2 | 3.000 |
| Linear Radau-FR/DG equivalence, maximum over degrees 0–2 | 4.44 × 10⁻¹⁵ |

The equivalence check implements Radau flux correction for **linear advection** and compares its nodal residual to the modal DG residual. It does not implement or certify a separate nonlinear FR cone solver.

The [six-run cone study](results-dg/comparison.json) uses Mach 2.35, cone 10°, outer boundary 45°, 40 cells, the same freestream initialization and HLLC/AUSM+-up2. Each degree has an 8,000-step maximum budget. Both degree-zero runs meet the residual and explicit-drift gates. The four degree-one/two runs **do not meet the full modal residual gate**. Their finite, positive fields are retained as exploratory results. Smooth order verification does not resolve shock limiting or prove high-order cone accuracy.

[Computed profiles, convergence and smooth-order figure](results-dg/comparison.png) · [Vector PDF](results-dg/comparison.pdf)

The residual includes higher modes. A limiter may create a nearly stationary time map while leaving a nonzero semidiscrete residual. We therefore retain the strict residual gate and do not relabel these cases as converged. Further work on shock stabilization, high-order boundary closure and steady convergence is needed before recommending DG as the default cone method.

## Run

```bash
python verify_dg.py
python dg_study.py
python dg.py --degree 2 --flux hllc --cells 40 --steps 8000 --output results-dg/my-run
```

The [single Colab notebook](Conical_Flow_Modern_2026.ipynb) embeds this module and lets students select `Physics = dg-inviscid`. It also retains all finite-volume and viscous choices.

## Primary references

- [Cockburn and Shu (2001), Runge–Kutta Discontinuous Galerkin Methods for Convection-Dominated Problems](https://helper.ipam.ucla.edu/publications/pcatut2005/pcatut_5477_reprint2.pdf).
- [Huynh (2007), A Flux Reconstruction Approach to High-Order Schemes Including Discontinuous Galerkin Methods](https://doi.org/10.2514/6.2007-4079).
- [NASA overview of the FR high-order unstructured-grid approach](https://ntrs.nasa.gov/citations/20160010066).

The implementation is independent Python/NumPy code. References describe the methods; neither a university solver nor a packaged CFD solver was copied or executed.
