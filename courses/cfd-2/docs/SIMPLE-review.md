# Review of jyoun35/SIMPLE

Reviewed upstream revision: [6a1d90684c9462aa48416894f147910ca20e238f](https://github.com/jyoun35/SIMPLE/tree/6a1d90684c9462aa48416894f147910ca20e238f). Review and execution date: October 6, 2026. The [validation report](validation/README.md) records the executed cases and their limits.

## What the program does

This is a teaching solver for **two-dimensional, steady, incompressible Navier-Stokes flow**. It stores horizontal velocity `u`, vertical velocity `v` and pressure `p` at cell centers on a collocated finite-volume mesh. Gmsh supplies the quadrilateral mesh; the code builds cell connectivity, face normals, face lengths and cell areas.

The discrete continuity equation requires the sum of mass fluxes through every cell's faces to vanish. The momentum equations balance convection, pressure gradients and viscous diffusion. Upwind interpolation supplies the convective terms; explicit non-orthogonal terms account for faces not perpendicular to the line joining neighboring cell centers.

SIMPLE means **Semi-Implicit Method for Pressure-Linked Equations**. Pressure has no separate time-evolution equation in this incompressible formulation. Instead, the algorithm adjusts pressure and velocity together so that the predicted flow better satisfies continuity.

## One iteration, step by step

1. **Assemble momentum operators.** The current velocity supplies convective face fluxes, while viscosity supplies diffusion. Boundary-condition objects contribute inlet, outlet and wall terms.
2. **Split diagonal and neighbor contributions.** The diagonal inverse and off-diagonal contribution provide a velocity predictor. This is a diagonal update using the previous field, rather than a full direct solve of the momentum system.
3. **Predict velocity.** The predictor generally does not satisfy continuity.
4. **Solve pressure correction.** A sparse pressure-correction equation is assembled and solved using SciPy `spsolve`.
5. **Correct pressure and velocity.** Relaxation factors damp changes to support iterative convergence.
6. **Measure continuity imbalance.** The code reports the maximum absolute cell residual, the sum of absolute cell residuals (L1), and the signed domain sum.
7. **Save fields and plots.** Checkpoints contain velocity and pressure arrays, contours, a continuity history and a surface-pressure-coefficient plot.

Rhie-Chow interpolation modifies the face velocity used in the continuity calculation. Its purpose is to avoid spurious checkerboard pressure modes on collocated grids. Non-orthogonal corrections are a separate geometric issue; both mechanisms must be consistent on skew meshes.

Sources: [main program](https://github.com/jyoun35/SIMPLE/blob/6a1d90684c9462aa48416894f147910ca20e238f/src/__main__.py), [momentum discretization](https://github.com/jyoun35/SIMPLE/blob/6a1d90684c9462aa48416894f147910ca20e238f/src/solution/discretize.py), [pressure-velocity algorithm](https://github.com/jyoun35/SIMPLE/blob/6a1d90684c9462aa48416894f147910ca20e238f/src/solution/SIMPLE.py).

## Default case and outputs

The supplied main program selects the fine NACA 2412 mesh, a speed of 34.3 m/s and an angle of attack of 3 degrees. It sets outlet pressure to 101325 Pa and viscosity to 1.79e-5, uses 1000 iterations, and saves every 10 iterations. Density is hard-coded as 1 in the reviewed formulation.

Velocity and pressure arrays are saved in `.npz` files. Figures show speed, pressure relative to the freestream, continuity imbalance, residual history and `Cp`. The `Cp` figure also plots supplied XFOIL data. Plotting a reference curve alongside a CFD result does not establish a quantitative validation.

The executed airfoil diagnostic used the upstream **very-coarse mesh and 100 iterations**, keeping the physical settings and iteration operators above. It did not execute the default fine-mesh 1000-iteration calculation. The exact uniform-flow test used the upstream rectangular mesh and 20 iterations. See [recorded evidence](validation/README.md).

## File map

| File or directory | Responsibility |
| --- | --- |
| `src/__main__.py` | Case settings, main iteration loop and saving |
| `src/mesh.py` | Gmsh import, connectivity and geometry |
| `src/solution/discretize.py` | Momentum operators and face fluxes |
| `src/solution/SIMPLE.py` | Predictor, pressure correction and solution update |
| `src/numerics/gradients.py` | Green-Gauss gradients with reconstructed nodal values |
| `src/numerics/interpolations.py` | Linear interpolation and Rhie-Chow face velocity |
| `src/numerics/non_orthogonal_correction.py` | Geometric decomposition of face fluxes |
| `src/boundaries/` | Velocity inlet, pressure outlet and slip-wall contributions |
| `src/helpers.py` | Residual evaluation, figures, pressure coefficient and checkpoints |

## Specific defects and the supporting evidence

### 1. Rhie-Chow does not preserve affine pressure on a skew configuration

The [interpolation function](https://github.com/jyoun35/SIMPLE/blob/6a1d90684c9462aa48416894f147910ca20e238f/src/numerics/interpolations.py) divides the cell pressure difference by `d dot n`, but projects the interpolated gradient onto the unit center-to-center direction. Those directions differ on a non-orthogonal grid.

The executed isolated test uses centers `P=(0,0)`, `N=(1,1)`, face normal `n=(1,0)`, exact pressure `p=x+y`, exact gradient `(1,1)`, zero input velocity and unit diagonal inverses. The original function returns face velocity `(-0.5857864376, 0)`; an affine field should leave the interpolation correction zero. The orthogonal control, with `N=(1,0)`, returns zero.

This was reproduced using the actual pinned source function. The gradient projection must be made consistent with the pressure-difference term and the chosen non-orthogonal decomposition before trusting skew-grid results.

### 2. Repeated non-orthogonal pressure corrections accumulate the right-hand side

In the [pressure-correction loop](https://github.com/jyoun35/SIMPLE/blob/6a1d90684c9462aa48416894f147910ca20e238f/src/solution/SIMPLE.py), `b_p += b_p_prime_non_ortho` retains corrections from previous passes. A fixed-point update normally uses the unchanged base right-hand side plus the current correction, rather than the cumulative sum of all preceding corrections.

This is a source-inspection finding. The default setting performs one correction; the cumulative behavior matters when the number of corrections is increased. No multi-corrector reference validation is claimed.

### 3. Boundary correction uses a stale interior-cell index

The boundary part of the same non-orthogonal loop selects boundary cell geometry using `cell`, but passes `grad_p_prime_x[cell_P]` and `grad_p_prime_y[cell_P]`. Here `cell_P` is left over from the preceding internal-face loop. It should select the gradient belonging to the current boundary cell.

This is a source-inspection finding about an index mismatch. A constant-pressure or orthogonal-grid test can hide its effect; passing the uniform-flow benchmark does not clear this defect.

### 4. The pressure-color limits are reversed

The upstream [checkpoint helper](https://github.com/jyoun35/SIMPLE/blob/6a1d90684c9462aa48416894f147910ca20e238f/src/helpers.py) declares `pmin=500` and `pmax=-500`. The intended increasing interval appears to be `[-500,500]`.

The executed function test with Matplotlib 3.11.2 saved the figure without an exception, but the colorbar changed the interval to `[-550,500]`. A universal plotting crash is not asserted. [Isolated diagnostic results](simple-diagnostics.json) record the observed behavior.

## Physical and numerical limits

- The airfoil uses [SlipWall](https://github.com/jyoun35/SIMPLE/blob/6a1d90684c9462aa48416894f147910ca20e238f/src/boundaries/SlipWall.py). It removes wall-normal penetration but does not impose no-slip. This case cannot validate a real viscous boundary layer or skin-friction drag.
- There is no compressible-flow energy equation or turbulence model in the reviewed solver.
- Density is fixed at 1; adding a density entry to the initial-condition dictionary does not change the discretization.
- The program monitors continuity, but supplies no automatic convergence stop or independent momentum-residual check.
- A small signed net flux can conceal large positive and negative cell residuals. Inspect L1 and maximum imbalance as well.
- The mesh reader reshapes 2D elements into four-node cells without validating Gmsh element types. Triangular and higher-order meshes require changes.
- Default Windows backslash paths need portable path handling on Linux.
- Checkpoints contain fields, but not a complete mesh/configuration provenance bundle or a validated restart procedure.
- The pinned tree contains no README or LICENSE. Upstream source is linked here and retained privately for testing, rather than mirrored into the course archive.

## Relationship to Mazaheri's projects

The three course projects emphasize compressible conservation equations, Roe/Van Leer/H-L-R fluxes, conical flow, airfoil flow, viscosity and mesh generation. SIMPLE supplements the **incompressible pressure-velocity coupling** part of the numerical-methods study. Read [Hejranfar chapter 6](../lectures/hejranfar/chapter-06.pdf), with chapters 4-5 for grid geometry and finite volume.

A useful validation sequence is an exact uniform state on an orthogonal mesh, affine-pressure interpolation on a skew mesh, an appropriate physical benchmark, and mesh refinement. The current evidence includes the first two checks and a finite airfoil diagnostic; it does not certify the airfoil solution.

## Reproduction

For isolated function diagnostics:

```text
python -m pip install numpy matplotlib
python tools/diagnose_simple.py --output work/simple-function-diagnostics.json
```

For solver execution, install NumPy, SciPy, Gmsh and Matplotlib, obtain the pinned upstream checkout, and follow the commands in the [validation report](validation/README.md). The harness calls the actual upstream numerical functions, records versions and source/mesh hashes, and saves field ranges and residual histories.
