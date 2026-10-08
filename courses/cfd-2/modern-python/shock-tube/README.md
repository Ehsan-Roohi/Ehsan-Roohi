# Shock-tube algorithm comparison

[Download the English lecture PDF](output/pdf/CFD_Algorithms_Lecture_2026.pdf): 33 teaching pages covering every implemented flux/reconstruction, DG, algorithm steps, advantages and limitations, worked examples, exercises with answers and measured shock-tube plots. A final reading map links the open LeVeque, Shu and Persson teaching resources that informed the explanations. [PDF source](teaching/build_lecture.py).

## Fixed uniform-grid refinement

[Six-grid report and plots](UNIFORM_GRID_RESULTS.md) compare all 15 shared fluxes and standalone JST on **80, 160, 320, 640, 1280 and 2560 fixed uniform cells**. The 96 grid/method comparisons show continuous cell-center curves, shock/contact details, physical properties, L1 errors and measured transition widths. The original 592-run study remains available below.


[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Ehsan-Roohi/Ehsan-Roohi/blob/main/courses/cfd-2/modern-python/shock-tube/Shock_Tube_All_Methods_2026.ipynb)

A separate results-first Python notebook for the numerical algorithms shared with [the cone project](../conical/README.md).
The study covers **15 shared pointwise Euler fluxes**, **20 FV reconstruction choices**, standalone **JST**, and **DG degrees 0/1/2**.
The notebook includes an English [educational algorithm guide](ALGORITHM_GUIDE.md): governing equations, update steps, every flux/reconstruction, DG, time integration, positivity controls and interpretation of convergence plots. Each method has an explanation of its mechanism and limitations.
The full **300-combination Sod matrix** is executed on 80 cells. Three-grid studies isolate fluxes (CWENO3), reconstructions (HLLC), and every DG degree/flux pair on 40, 80 and 160 elements. Lax, double rarefaction and stationary contact provide additional wave/positivity checks.

[Download source, results and figures](../../downloads/shocktube-modern-2026.zip) · [Executed report](RESULTS.md) · [All numerical metrics](run-metrics.csv) · [Notebook execution evidence](colab-validation.json) · [Algorithm and provenance audit](../conical/ADDITIONAL_METHODS.md)

Setup is collapsed with no encoded source cell. A hashed public ZIP downloads the shared editable Python modules and recorded numerical data. Recorded answers display by default; an optional fresh run can recompute one method, all fluxes, or the full matrix. No packaged CFD solver is used.

Completion means the run reached the requested physical time and passed finite/positive cell-average and boundary-budget conservation checks. It is not a guarantee of accurate wave resolution. Stopped runs retain actual times and are excluded from final-time accuracy errors. This benchmark is inviscid Euler; physical viscous transport is not tested here.

## Run from source

The numerical core is shared in [conical/shock_suite.py](../conical/shock_suite.py), not duplicated as a different implementation.

```bash
cd ../conical
python shock_suite.py --case sod --flux hllc --reconstruction weno5-z --cells 160
python shock_suite.py --case sod --flux ausm-up2 --degree 2 --cells 80
python shock_study.py
python uniform_refinement.py
python verify_additional.py
```

[Raw run records and fields](../conical/results-shocktube/catalog.json) include hashes, configuration, actual final time, errors, positivity, conservation budgets, rejected stages and reference-quadrature changes. Source code, diagrams and documentation are in English.
