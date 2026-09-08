# Research portfolio

[← Ehsan Roohi on GitHub](README.md) · [UMass research page](https://websites.umass.edu/roohie/research/) · [Google Scholar](https://scholar.google.com/citations?user=AWKLce4AAAAJ&hl=en) · [ORCID](https://orcid.org/0000-0001-5739-3210) · [YouTube channel](https://www.youtube.com/@Ehsan_Roohi)

My research connects molecular simulation, rarefied-gas dynamics, computational fluid dynamics and scientific machine learning. This page is a visual guide to selected directions; it is not a substitute for the papers, and a video or repository image is not presented as independent validation.

## Scientific machine learning for compressible and rarefied flows

The central objective is not simply to replace a solver with a neural network. It is to identify where learned representations can reduce repeated computational cost while preserving the quantities that carry physical meaning: conservation, shock location, internal-energy relaxation, boundary loads and non-equilibrium transport.

<table>
<tr>
<td width="48%" valign="top"><a href="https://www.youtube.com/watch?v=QBXaqeFriUk"><img src="https://img.youtube.com/vi/QBXaqeFriUk/hqdefault.jpg" width="100%" alt="Shock-aware artificial intelligence"><br><b>Watch: Shock-aware AI</b></a><br>A short introduction to why sharp fronts need physically meaningful evaluation rather than a single global error.</td>
<td width="48%" valign="top"><a href="https://www.youtube.com/watch?v=S6U-0hBSJ98"><img src="https://img.youtube.com/vi/S6U-0hBSJ98/hqdefault.jpg" width="100%" alt="Artificial intelligence for rarefied gas dynamics"><br><b>Watch: AI for rarefied gas dynamics</b></a><br>An overview of learned models for computationally demanding kinetic-flow analysis.</td>
</tr>
</table>

Selected papers and open artifacts:

- [Neural networks for rarefied gas dynamics: relaxation, polyatomic shocks and hypersonic-cylinder flow](https://doi.org/10.1063/5.0334590), *Physics of Fluids* **38**, 057108 (2026).
- [Physics-constrained neural collision operators for hard-sphere surrogates and ab initio angle prediction in DSMC](https://doi.org/10.1063/5.0328463), *Physics of Fluids* **38**, 057123 (2026).
- [Micro-step DeepONet: paper-linked data, code and reference results](https://github.com/Ehsan-Roohi/roohi-step-dnn-mahdavi).
- [Micro-nozzle POD and shock-aligned surrogate reproducibility](https://github.com/Ehsan-Roohi/roohi-nozzle-pod-reproducibility).
- [ShockVortexML](https://github.com/Ehsan-Roohi/ShockVortexML), ongoing work on joint shock and vortex-core identification with explicit physical audits.

## Hypersonic DSMC and shock structure

My hypersonic research connects particle collision algorithms to the prediction of shock layers and non-equilibrium flow around aerodynamic bodies. These studies complement the newer machine-learning work: the underlying DSMC calculations are a research direction in their own right.

### Hypersonic cylinder

<p align="center"><a href="https://websites.umass.edu/roohie/files/2025/05/image.jpg"><img src="https://websites.umass.edu/roohie/files/2025/05/image.jpg" width="800" alt="Hypersonic cylinder flow-field illustration from Ehsan Roohi's UMass research gallery"></a></p>

Hypersonic-cylinder flow visualization from my [UMass research gallery](https://websites.umass.edu/roohie/research/). The cylinder connects collision-sampling research with shock-layer prediction. This archived illustration is distinct from the later [2026 neural-network cylinder study](https://doi.org/10.1063/5.0334590); it is not labelled as a result of that paper.

### Hypersonic biconic geometry

<p align="center"><a href="https://websites.umass.edu/roohie/files/2025/05/image-1-1024x791-1-1-1.png"><img src="https://websites.umass.edu/roohie/files/2025/05/image-1-1024x791-1-1-1.png" width="800" alt="Shock-wave structure around a hypersonic biconic geometry"></a></p>

Shock-wave structure around a biconic body, associated on the UMass page with Goshayeshi, Roohi and Stefanov's [DSMC simulation of hypersonic flows using an improved SBT-TAS technique](https://doi.org/10.1016/j.jcp.2015.09.027), *Journal of Computational Physics* **303**, 28–44 (2015). This work links adaptive collision treatment to hypersonic-flow simulation. [Figure source](https://websites.umass.edu/roohie/research/).

### Related high-speed CFD visualization

Shock-bearing flows are useful stress tests for both numerical methods and learned models. Local front position, topology, surface quantities and conservation can reveal failures hidden by field-averaged scores.

<p align="center"><a href="https://www.youtube.com/watch?v=SQxV9lnMc1I"><img src="https://img.youtube.com/vi/SQxV9lnMc1I/hqdefault.jpg" width="720" alt="Mach 3 shock formation over a diamond airfoil"><br><b>Watch: Mach 3 shock formation over a diamond airfoil</b></a></p>

The video is a visualization of an MFC CFD schlieren calculation. It should be read as a simulation demonstration, not as a grid-convergence or experimental-validation record. Related open work includes [SU2 Diamond Airfoil Verification](https://github.com/Ehsan-Roohi/SU2-Diamond-Airfoil-Verification) and the shock/vortex project above; each repository states its own evidence scope.

## DSMC collision algorithms and non-equilibrium transport

This long-running direction develops and evaluates particle collision procedures and uses DSMC to study flows for which continuum constitutive assumptions become unreliable.

### Collision-scheme development

<p align="center"><a href="https://websites.umass.edu/roohie/files/2025/05/5-5-1024x562-1.png"><img src="https://websites.umass.edu/roohie/files/2025/05/5-5-1024x562-1.png" width="900" alt="Timeline of DSMC collision schemes including contributions by Roohi and collaborators"></a></p>

The timeline places our group's contributions within the wider development of DSMC collision schemes. My collaborative work includes intelligent SBT and adaptive-subcell methods, GBT, SSBT and SGBT. Earlier NTC, BT and SBT methods are background contributions by their respective developers, not claimed as our inventions. [Original timeline and descriptions](https://websites.umass.edu/roohie/research/).

### How collision partners are selected

<p align="center"><a href="https://websites.umass.edu/roohie/files/2025/05/6-6-1024x753-1.png"><img src="https://websites.umass.edu/roohie/files/2025/05/6-6-1024x753-1.png" width="800" alt="Comparison of collision-partner selection in BT, SBT and SSBT schemes"></a></p>

BT, SBT and SSBT use different candidate-selection procedures. The diagram illustrates those differences; collision probabilities and algorithmic details belong to the original formulations. See the [2016 review](https://doi.org/10.1016/j.physrep.2016.08.002), the [2018 GBT paper](https://doi.org/10.1016/j.jcp.2017.10.033), and the [2022 SSBT paper](https://pubs.aip.org/aip/pof/article-abstract/34/1/012010/2845513/A-symmetrized-and-simplified-Bernoulli-trial). [Figure source](https://websites.umass.edu/roohie/research/).

### Simplified generalized Bernoulli trials (SGBT)

<p align="center"><a href="https://websites.umass.edu/roohie/files/2025/05/7-7-908x1024-1.png"><img src="https://websites.umass.edu/roohie/files/2025/05/7-7-908x1024-1.png" width="640" alt="SGBT collision-selection procedure attributed to Javani and collaborators, 2024, on the UMass page"></a></p>

The UMass gallery attributes this procedure to Javani and collaborators (2024): SGBT applies a symmetrized selection procedure to a selected subset of particles. The diagram provides an algorithmic overview, not a standalone implementation specification. [Source and attribution](https://websites.umass.edu/roohie/research/).

### Sensitivity to particles per cell

<p align="center"><a href="https://websites.umass.edu/roohie/files/2025/05/8-8-1024x841-1.png"><img src="https://websites.umass.edu/roohie/files/2025/05/8-8-1024x841-1.png" width="800" alt="Heat-transfer comparison examining DSMC collision schemes at different particle counts per cell"></a></p>

The displayed heat-transfer benchmark examines sensitivity to particles per cell, a practical concern for computational cost and collision sampling. The low-particle-count behavior shown here is specific to the reported test, not a universal accuracy guarantee for arbitrary flows. [Original figure and study description](https://websites.umass.edu/roohie/research/).

Selected publications and implementations:

- [An open source, parallel DSMC code for rarefied gas flows in arbitrary geometries](https://doi.org/10.1016/j.compfluid.2010.07.014), *Computers & Fluids* (2010).
- [Collision partner selection schemes in DSMC: From micro/nano flows to hypersonic flows](https://doi.org/10.1016/j.physrep.2016.08.002), *Physics Reports* (2016).
- [A generalized form of the Bernoulli Trial collision scheme in DSMC](https://doi.org/10.1016/j.jcp.2017.10.033), *Journal of Computational Physics* (2018).
- [Generalised Bernoulli trials for OpenFOAM](https://github.com/Ehsan-Roohi/openfoam-generalised-bernoulli-trials) and [DSMC in Python](https://github.com/Ehsan-Roohi/DSMC_Python).

## Micro- and nanoscale gas flows

Rarefaction introduces velocity slip, temperature jump, thermal polarization and transport behavior that cannot be inferred safely from macroscale intuition. The work spans pressure-, shear- and thermally driven flows, porous structures, nozzles and aerodynamic configurations.

<p align="center"><a href="https://www.youtube.com/watch?v=FGuTGyQwhS0"><img src="https://img.youtube.com/vi/FGuTGyQwhS0/hqdefault.jpg" width="640" alt="Introduction to micro and nano flows"><br><b>Watch: Introduction to Micro and Nano Flows</b></a></p>

For a consolidated review, see [A comprehensive review on micro- and nano-scale gas flow effects](https://doi.org/10.1016/j.physrep.2022.10.004), *Physics Reports* **997** (2023). The 2025 book [*Advances in Direct Simulation Monte Carlo*](https://doi.org/10.1007/978-981-96-8200-3) develops the broader kinetic-theory and application context.

## Open research and teaching

[FlowMLLab](https://github.com/Ehsan-Roohi/FlowMLLab) turns selected, properly attributed research cases into reproducible lecture notes and executable audits. The course code, newly written baselines and teaching claims are kept distinct from the original paper solvers and article results. Case-level origin and reuse status are recorded in its [data provenance statement](https://github.com/Ehsan-Roohi/FlowMLLab/blob/main/DATA_PROVENANCE.md).

For collaboration or questions, contact [roohie@umass.edu](mailto:roohie@umass.edu).
