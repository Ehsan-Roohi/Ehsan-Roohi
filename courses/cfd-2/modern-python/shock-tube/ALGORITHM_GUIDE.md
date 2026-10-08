# How the implemented CFD methods work

This guide explains **every selectable method in this notebook**: 15 Euler face fluxes, standalone JST, 20 finite-volume reconstructions and DG degrees 0/1/2. These are numerical algorithms for the same inviscid equations, not different physical turbulence or viscosity models. The descriptions follow the shared Python implementation; an algorithm name alone does not specify all its parameters.

## A. The common finite-volume algorithm

For planar flow, the Euler equations conserve mass, momentum and total energy:

$$U=(\rho,\rho u,E)^T,\qquad F(U)=(\rho u,\rho u^2+p,u(E+p))^T,$$
$$p=(\gamma-1)\left(E-\tfrac12\rho u^2\right),\qquad a=\sqrt{\gamma p/\rho}.$$

The code also carries transverse momentum; its velocity is zero in Sod. Each uniform cell stores a **conserved cell mean**, not a point sample. Its update is

$$L_i(U)=-\frac{\widehat F_{i+1/2}-\widehat F_{i-1/2}}{\Delta x}.$$

One time step consists of the following operations:

1. Extend boundary cell states into three ghost cells (outflow conditions).
2. Reconstruct left/right face states from neighboring cell means; control negative reconstructed density/pressure.
3. Compute one common numerical flux at each face. The same flux leaves one cell and enters its neighbor, giving discrete conservation.
4. Apply SSPRK3 stages, using a CFL step based on the largest acoustic speed. Reject and halve the step if a stage produces inadmissible cell means.

**Flux versus reconstruction:** the flux decides how waves cross a face; reconstruction decides how accurately face states represent each cell. Changing one does not change the other. JST and DG have separate spatial operators, explained below.

## B. Face-flux algorithms — all 15 choices

For left/right reconstructed states, let $F_L=F(U_L)$ and $F_R=F(U_R)$.

| Selection | Main idea and algorithm | Interpretation and limitation |
|---|---|---|
| `godunov` — exact Riemann | Solve the ideal-gas star-pressure equation by bracketing and bisection; construct shock/rarefaction/contact waves; sample at $x/t=0$ and evaluate its physical flux. Includes a vacuum branch. | Exact **local** Riemann flux; the complete spatial/time discretization still has error. More expensive than approximate solvers. |
| `roe` — Roe with entropy correction | Form square-root-density-weighted velocity/enthalpy; decompose the state jump into acoustic, contact and shear eigenwaves; subtract half the absolute-wave-speed weighted jump from the mean physical flux. Smooth small acoustic eigenvalues with a Harten-style entropy correction. | Resolves characteristic waves efficiently. The entropy correction is not a proof of full-scheme entropy stability; positivity is not automatic. |
| `roe-nc` — uncorrected Roe | Use the same Roe decomposition without the small-eigenvalue entropy correction. | A diagnostic comparison: expansion waves can produce nonphysical entropy behavior. |
| `hlle` — two-wave HLL | Bound the fan with $S_L=\min(u_L-a_L,u_R-a_R,0)$ and $S_R=\max(u_L+a_L,u_R+a_R,0)$; use the HLL intermediate flux between these bounds and the upwind physical flux outside them. | This implementation uses **Davis acoustic bounds**, despite the `hlle` key. Omitting the contact wave increases contact diffusion. |
| `hllc` — contact-restoring HLL | Estimate the two outer wave speeds and the contact speed; construct left/right star states; select the flux by the wave-speed signs. Fall back locally to HLL when star density/pressure is invalid or nonfinite. | Better contact resolution than two-wave HLL; moving contacts still spread on finite grids. |
| `rusanov` — local Lax–Friedrichs | Use $\widehat F=(F_L+F_R)/2-\alpha(U_R-U_L)/2$, with $\alpha=\max(\vert u_L\vert +a_L,\vert u_R\vert +a_R)$ at each face. | Simple acoustic-speed dissipation; generally more diffusive on slow/contact waves. |
| `global-lf` — global LF flux | Use the same central-minus-jump formula with a **single maximum speed over all faces in the current residual evaluation**. | Extra diffusion in regions whose local speeds are small. This is a semidiscrete flux, not the classical fully discrete LF update with $\Delta x/\Delta t$. |
| `vanleer` — Van Leer flux-vector splitting | Split each state's flux into positive/negative parts using Mach-number polynomials for subsonic flow; add $F^+(U_L)+F^-(U_R)$. Supersonic states contribute only in their upwind direction. | Smooth sonic splitting, with contact diffusion. This is distinct from the Van Leer **slope limiter**. |
| `steger-warming` — eigenvalue flux splitting | Split the physical-flux eigenvalues as $\lambda^\pm=(\lambda\pm\vert \lambda\vert )/2$ for each state; combine the left positive and right negative flux contributions. | Straightforward upwind splitting; can diffuse contacts. |
| `ausm` — original AUSM | Split mass advection by left/right local Mach polynomials and split pressure separately. Upwind velocity and enthalpy according to the interface mass-flux sign; add face pressure to normal momentum flux. | Separates transport and pressure; the original splitting is a baseline for its improved family. |
| `ausm-plus` — AUSM+ | Use an enthalpy-based interface sound speed and improved Mach/pressure polynomials, then the same mass-advection/pressure construction. | Improved splitting around sonic transitions; no explicit AUSM+-up low-Mach corrections. |
| `ausm-up` — AUSM+-up | Add a Mach-regulated pressure-jump correction to interface mass flux and a velocity-jump correction to interface pressure. The code uses reference Mach **0.1**. | Designed to improve pressure/velocity coupling across Mach regimes. Performance depends on the correction parameters and test. |
| `ausm-up2` — AUSM+-up2 | Keep the pressure-jump correction in the mass flux; replace the interface-pressure expression with a velocity-magnitude-based dissipation term. The implemented reference Mach remains **0.1**. | A different pressure-dissipation choice within the AUSM family; not an additional physical transport model. |
| `slau2` — SLAU2 | Use arithmetic mean sound speed, density-weighted velocity splitting and a low-Mach pressure-jump correction for mass flux; use velocity-magnitude-based interface pressure dissipation. | An all-speed splitting alternative. No optional localized artificial-diffusion sensor is included here. |
| `ec-lf` — entropy-conservative core + LF | Construct a Chandrashekar-style central flux using logarithmic means of density and $\beta=\rho/(2p)$; add local LF jump dissipation. | The central face entropy identity is checked separately. This does not certify entropy stability of the entire reconstructed, time-integrated, boundary-treated scheme. |

**HLL formula used here:**

$$\widehat F_{HLL}=\frac{S_RF_L-S_LF_R+S_LS_R(U_R-U_L)}{S_R-S_L}.$$

**AUSM family in one expression:** a signed mass flux transports upwind velocity/enthalpy, while a separately split pressure acts on normal momentum. Its variants change interface sound speed, Mach polynomials and pressure/velocity dissipation; they do not solve different governing equations.

### Standalone `jst` — Jameson–Schmidt–Turkel

JST uses a uniform-grid central stencil rather than a freely interchangeable two-state Riemann flux. Evaluate the central physical flux at the averaged state, $F((U_L+U_R)/2)$. Subtract artificial dissipation built from second and fourth state differences. A pressure-curvature sensor increases second-difference dissipation near shocks and reduces fourth-difference dissipation there:

$$\epsilon_2=0.5\,s_p,\qquad\epsilon_4=\max(0,0.02-\epsilon_2).$$

The second-difference term controls jumps; the fourth-difference term supplies gentler smooth-region damping. This is **numerical** dissipation, not molecular viscosity. JST uses its own cell-mean stencil and is not included in the 15-flux reconstruction matrix or DG options.

## C. Finite-volume reconstruction — all 20 choices

These reconstructions fit **cell-average moments**. The WENO/TENO algorithms are finite-volume face reconstructions, not finite-difference formulas applied to nodal values. Formal order refers to smooth regions with suitable boundaries and inactive limiting; discontinuities reduce measured global order.

| Selection | Main idea and algorithm | Tradeoff |
|---|---|---|
| `first` | Keep each cell's mean constant; its two face states equal that mean. | First order, inexpensive and diffusive. |
| `muscl` | Fit a linear polynomial with the monotonized-central (MC) slope: limit the centered slope by twice either one-sided difference; set it to zero when the differences disagree in sign. | Second order in smooth regions; reduces slopes around jumps/extrema. |
| `muscl-minmod` | Use the smaller magnitude of the two same-sign one-sided slopes, or zero when signs differ. | More conservative limiting, with greater diffusion. |
| `muscl-vanleer` | Use the harmonic slope $2ab/(a+b)$ for same-sign differences $a,b$, otherwise zero. | Smooth limiter response; not the Van Leer face flux. |
| `muscl-superbee` | For same-sign differences, choose the larger of the two bounded slopes $\min(2\vert a\vert,\vert b\vert)$ and $\min(\vert a\vert,2\vert b\vert)$, retaining their sign. | Sharper, more compressive profiles; not universally the most accurate choice. |
| `muscl-vanalbada` | Use the smooth slope $ab(a+b)/(a^2+b^2+\varepsilon)$ for same-sign differences, otherwise zero. | Smooth transition between limiting strengths. |
| `eno2` | Choose the one-sided linear stencil with the smaller neighboring state jump, component by component. | Second-order smooth reconstruction; abrupt stencil selection can affect profiles. |
| `cweno3` | Build a third-order optimal polynomial from cell means, two lower-order candidates and a central candidate; combine them with smoothness-dependent nonlinear weights. | One polynomial per cell supports both face evaluation and source quadrature. |
| `cweno5` | Apply the same central-WENO construction with a fifth-order optimal polynomial and wider candidates. | Higher smooth accuracy, more stencil work; limiting can reduce order. |
| `cweno3-char` | Project CWENO3 candidates into the local Euler eigenbasis; compute nonlinear weights per characteristic wave; transform back. | Limits coupling between different wave families, with extra eigenbasis work. |
| `cweno5-char` | Use characteristic weighting with the fifth-order CWENO construction. | Same purpose at a higher formal smooth order. |
| `weno3-js` | Form characteristic face values on small FV candidate stencils; blend using Jiang–Shu weights $\alpha_k=d_k/(\beta_k+\varepsilon)^2$. | Third-order smooth faces; nonlinear weights can lose accuracy near smooth critical points. |
| `weno3-z` | Use the same third-order candidates with Z weights $\alpha_k=d_k[1+(\tau/(\beta_k+\varepsilon))^2]$, where $\tau=\vert \beta_0-\beta_{last}\vert $. | Better recovery of optimal smooth weights; still limited near jumps. |
| `weno5-js` | Use fifth-order characteristic FV face candidates with Jiang–Shu weights. | Wider smooth reconstruction; critical-point and boundary effects remain. |
| `weno5-z` | Use fifth-order candidates and the Z smoothness correction. | An alternative smooth-weight recovery at fifth order. |
| `weno7-js` | Use seventh-order characteristic FV face candidates and Jiang–Shu weights. | Higher formal face order and wider stencil; no WENO7-Z option is implemented. |
| `teno3` | Normalize inverse-smoothness indicators, discard stencils below a cutoff, then renormalize the optimal linear weights over retained stencils. | An educational third-order FV cutoff variant; binary stencil selection. |
| `teno5` | Use a fifth-order TENO smoothness detector based on $\tau_5$, a binary stencil cutoff, and retained optimal weights. | Targets smooth-stencil accuracy while excluding troubled candidates. |
| `teno7` | Apply inverse-smoothness binary cutoff to seventh-order characteristic FV candidates. | An educational seventh-order variant; not TENO-A or a localized artificial-diffusion extension. |
| `muscl-thinc-bvd` | Construct MC-linear and bounded hyperbolic-tangent interface candidates; compare face-jump variation with neighbor candidates; choose the smaller variation component by component. Use interface steepness **1.6** and positivity controls. | An educational **single-stage** BVD selection; may sharpen interfaces, but is not a Riemann solver or a published multistage P4-THINC-BVD implementation. |

Here $d_k$ are optimal linear weights and $\beta_k$ measure candidate smoothness through polynomial derivatives. CWENO weights use a mesh- and component-scaled $\varepsilon$. TENO uses cutoff **$10^{-5}$** and exponent **6**. WENO/TENO face reconstructions use characteristic variables; plain CWENO and the listed linear/ENO/BVD choices work componentwise, unless `-char` is selected.

In the cone code, face-only WENO/TENO also require a separate CWENO source polynomial: CWENO3 for third-order faces and CWENO5 for fifth/seventh-order faces. Therefore seventh-order **face reconstruction** does not imply seventh-order accuracy for the entire cone solver. Planar Sod has no geometric source.

## D. Discontinuous Galerkin — degrees 0, 1 and 2

DG stores a local polynomial as Legendre coefficients rather than only one cell mean. Multiply Euler's equations by each basis function and integrate by parts: volume quadrature handles the physical flux, and the same Riemann flux couples neighboring elements at their boundaries. This notebook uses outflow boundaries and no cone geometric source.

| Option | Representation and relationship to FV | Accuracy qualification |
|---|---|---|
| DG degree 0 | One constant coefficient per element. With the same flux and time step, this reproduces first-order FV; all 15 flux checks agree to roundoff. | First-order spatial representation. |
| DG degree 1 | Mean plus one linear Legendre mode; evolve both weak-form coefficients. | Formally second order in smooth regions; shock limiting can reduce it. |
| DG degree 2 | Mean, linear and quadratic modes; evolve all three using volume quadrature and face fluxes. | Formally third order in smooth regions; degree alone does not guarantee better shock profiles. |

Each degree supports all **15 pointwise fluxes**, with JST excluded. Gauss quadrature uses at least $\max(3,p+2)$ points. A characteristic TVB slope limiter uses neighboring means; if it changes the linear mode, higher modes are zeroed. The recorded runs use **TVB parameter 0**, giving strong limiting. Sampled positivity scaling reduces higher modes toward the mean while preserving that mean. The time-step bound includes the factor $1/(2p+1)$.

## E. Time integration, admissibility and reading the plots

**SSPRK3** advances the spatial operator through three stages:

$$A=U^n+\Delta tL(U^n),$$
$$B=\tfrac34U^n+\tfrac14[A+\Delta tL(A)],$$
$$U^{n+1}=\tfrac13U^n+\tfrac23[B+\Delta tL(B)].$$

For FV, $\Delta t=\mathrm{CFL}\,\Delta x/\max(|u|+a)$, clipped to the remaining physical time; recorded uniform-grid runs use **CFL 0.25**. SSP terminology does not guarantee nonoscillatory or positive results for every reconstruction and flux. Reconstructed states are scaled toward admissible means at sampled points; invalid stage means trigger step halving, with a finite retry limit. Failed runs retain their actual stopping time.

**Fixed-grid refinement:** repeat the same algorithm on $N=80,160,320,640,1280,2560$ cells with no mesh adaptation. Decreasing $\Delta x$ reduces the physical thickness of captured shocks/contacts; it should not remove the smooth physical rarefaction fan. A finite grid still has a numerical transition. The six-grid comparison holds CWENO3 fixed for all 15 fluxes and uses JST's own stencil; it is not a six-grid run of every reconstruction/DG combination.

**What the diagnostics measure:**

- Physical profiles show density, velocity, pressure and derived quantities from computed conserved means. Connecting cell centers changes display style, not numerical accuracy.
- $L_1$ mesh errors compare numerical means against the exact Riemann solution integrated into **conserved** cell means, then converted to primitive variables. Across discontinuities, global measured order can be lower than formal smooth order.
- Shock/contact 10–90% density widths measure physical transition thickness; compare both width and width in cells. They are not a universal ranking for all wave problems.
- Time histories show transient residuals/admissibility and the SSPRK-weighted boundary-flux conservation budget. Transient RHS decay is not steady-state convergence.
- Completion, finite positive states, conservation and file hashes establish execution/data integrity. Accuracy needs the separate exact-solution errors and wave-resolution comparisons.

## F. Implementation and reference trail

The methods above are implemented in the shared, editable modules: `shock_suite.py` (FV/DG driver), `flux.py`, `additional_fluxes.py`, `ausm.py`, `entropy_flux.py`, `reconstruction.py` and `advanced_reconstruction.py`. The original cone solver supplies reusable numerical components; the planar shock-tube driver removes cone geometry. See the [implementation and provenance audit](../conical/ADDITIONAL_METHODS.md) for verification scope and detailed attribution.

Related reference implementations and papers: [NASA AUSM shock tube](https://github.com/nasa/shocktube), [Euler comparison framework](https://github.com/fhermet/euler-1d-solver), [Riemann solvers](https://github.com/cangyu/Riemann-Solvers), [CAELUM](https://github.com/navasmontilla/CAELUM), [Quail](https://github.com/IhmeGroup/quail), and [boundary variation diminishing reconstruction](https://arxiv.org/abs/1602.00814).
