# Project 2 - compressible airfoil

The archived `AIRFOIL` programs evolve density, two velocity components, pressure and energy, with Roe-related flux routines and multistage update coefficients. They are compressible-flow codes, whereas the separately reviewed SIMPLE solver is incompressible.

[Root source variants](variants) are preserved separately from five case directories: [85-medium](cases/85-medium), [85-high](cases/85-high), [very-high](cases/very-high), [exact](cases/exact), and [uniform](cases/uniform).

The directory names are historical labels, not a guarantee of Mach number, resolution or exactness. Inspect each source's `Input` routine. For example the root `airfoilfinal.for` sets `XMIN=0.85` and an angle of 1 degree; root `airfoil4.for` sets `XMIN=0.80` and 1.25 degrees. Do not interchange these settings based on names alone.

The selected 85-medium case opens `Node.DAT`, `connectivity.DAT`, `BC.DAT` and `BE.DAT`. It and several other variants import `MSIMSL` / `MSFLIB`; a modern compiler needs a deliberate compatibility review. Some root variants reference boundary files that were not co-located with them, so root variants are not presented as ready-to-run cases.

[Project brief](../../reports/02-compressible-airfoil/Project%20II.doc) · [Report](../../reports/02-compressible-airfoil/Report%202.doc) · [Mach 2 appendix](../../reports/02-compressible-airfoil/Appendix%20M%3D2%20Results.doc) · [Selected old results](../../results/README.md).
