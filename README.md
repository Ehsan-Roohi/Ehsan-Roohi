# Ehsan Roohi

**Rarefied gas dynamics · Molecular simulation · Scientific machine learning**

I am a [Lecturer in Mechanical and Industrial Engineering at the University of Massachusetts Amherst](https://www.umass.edu/engineering/about/directory/ehsan-roohi). My research connects kinetic theory, computational fluid dynamics and machine learning, with applications spanning micro/nanoscale gas transport and hypersonic flows.

My work began with molecular–continuum modeling and direct simulation Monte Carlo (DSMC), and now also explores learned collision models, physics-informed neural networks and flow-field surrogates. A central question across these efforts is how to reduce computational cost while retaining the physical behavior that matters: conservation, fluctuations, transport and non-equilibrium effects.

[UMass website](https://websites.umass.edu/roohie/) · [Google Scholar](https://scholar.google.com/citations?user=AWKLce4AAAAJ&hl=en) · [ORCID](https://orcid.org/0000-0001-5739-3210) · [Email](mailto:roohie@umass.edu)

## Research areas

- **Rarefied gas dynamics and DSMC:** collision algorithms, particle methods and non-equilibrium transport.
- **Micro- and nanoscale flows:** slip and temperature jump, thermally driven transport and Knudsen pumps.
- **Hypersonic and compressible flows:** shock structure, internal-energy relaxation and gas–surface quantities.
- **Scientific machine learning:** physics-constrained collision models, PINNs, operator learning and reduced-order representations.
- **Verification and reproducibility:** numerical benchmarks, explicit data provenance and comparisons with classical methods.

## Selected research

### Molecular simulation and collision algorithms

Understanding and improving collision sampling has been a long-running part of my research, from open-source DSMC implementations to the design and assessment of collision-partner selection schemes.

- Scanlon, Roohi, White, Darbandi & Reese (2010). [An open source, parallel DSMC code for rarefied gas flows in arbitrary geometries](https://doi.org/10.1016/j.compfluid.2010.07.014). *Computers & Fluids*.
- Roohi & Stefanov (2016). [Collision partner selection schemes in DSMC: From micro/nano flows to hypersonic flows](https://doi.org/10.1016/j.physrep.2016.08.002). *Physics Reports*.
- Roohi, Stefanov, Shoja-Sani & Ejraei (2018). [A generalized form of the Bernoulli Trial collision scheme in DSMC: Derivation and evaluation](https://doi.org/10.1016/j.jcp.2017.10.033). *Journal of Computational Physics*.

### Micro/nanoflows and non-equilibrium transport

This direction examines gas-flow behavior beyond the usual continuum assumptions, including thermally driven motion and the relationship between molecular transport and macroscopic observables.

- Akhlaghi, Roohi & Stefanov (2023). [A comprehensive review on micro- and nano-scale gas flow effects: Slip-jump phenomena, Knudsen paradox, thermally-driven flows, and Knudsen pumps](https://doi.org/10.1016/j.physrep.2022.10.004). *Physics Reports*.

### Machine learning for kinetic and rarefied flows

Recent studies investigate neural representations of collision dynamics and complete flow fields. These are related but distinct tasks: reproducing a field does not, by itself, establish the stability or conservation properties of a learned collision model.

- Roohi, Shoja-Sani & Ebrahimzadeh Azghadi (2026). [Neural networks for rarefied gas dynamics: Relaxation problem, polyatomic shock waves, and hypersonic cylinder flow](https://doi.org/10.1063/5.0334590). *Physics of Fluids* **38**, 057108.
- Roohi, Shoja-sani & Stefanov (2026). [Physics constrained neural collision operators for hard sphere surrogates and ab initio angle prediction in direct simulation Monte Carlo](https://doi.org/10.1063/5.0328463). *Physics of Fluids* **38**, 057123.
- Roohi & Mahdavi (2026). [Analysis of the rarefied flow at micro-step using a DeepONet surrogate model with a physics-guided zonal loss function](https://doi.org/10.1007/s10404-026-02899-8). *Microfluidics and Nanofluidics*.

For the broader publication record, see [Google Scholar](https://scholar.google.com/citations?user=AWKLce4AAAAJ&hl=en) and [Books and Papers](https://websites.umass.edu/roohie/books-and-papers/).

## Software and open research artifacts

| Project | What you will find |
| --- | --- |
| [FlowMLLab](https://github.com/Ehsan-Roohi/FlowMLLab) | A reproducible CFD-to-SciML course: lecture notes, executable notebooks, baselines and physical validation. [Course site](https://ehsan-roohi.github.io/FlowMLLab/) |
| [Generalised Bernoulli trials for OpenFOAM](https://github.com/Ehsan-Roohi/openfoam-generalised-bernoulli-trials) | An archival DSMC collision implementation with upstream attribution and GPL licensing. |
| [Micro-step DeepONet](https://github.com/Ehsan-Roohi/roohi-step-dnn-mahdavi) | Curated data, code and reference results accompanying the micro-step study. |
| [Micro-nozzle POD reproducibility](https://github.com/Ehsan-Roohi/roohi-nozzle-pod-reproducibility) | Research artifacts for shock-centered low-rank analysis of rarefied micro-nozzle flows. |
| [DSMC in Python](https://github.com/Ehsan-Roohi/DSMC_Python) | Python notebooks for direct simulation Monte Carlo. |
| [Introduction to Compressible Flows](https://github.com/Ehsan-Roohi/Introduction-to-Compressible-Flows) | Chapter-linked computational notebooks for gas dynamics. |

Research archives and teaching implementations have different purposes. Where earlier article data are brought into FlowMLLab, the original paper and dataset are credited separately from newly written loaders, baselines and exercises. See the course's [data provenance statement](https://github.com/Ehsan-Roohi/FlowMLLab/blob/main/DATA_PROVENANCE.md) and each repository's own reuse terms and validation notes.

## Current research and preprints

- [Geometry-native machine learning reconstruction of DSMC moment fields with support monitoring](https://arxiv.org/abs/2609.01637) — 2026 preprint on reconstructing noisy DSMC moment fields. Listed as a preprint, not as a published journal article.
- [ShockVortexML](https://github.com/Ehsan-Roohi/ShockVortexML) — ongoing research on joint shock and vortex-core identification, with physical audits and documented evaluation limits.

## Books and teaching

With Hassan Akhlaghi and Stefan Stefanov, I coauthored [*Advances in Direct Simulation Monte Carlo: From Micro-Scale to Rarefied Flow Phenomena*](https://doi.org/10.1007/978-981-96-8200-3), Springer, 2025.

My [teaching](https://websites.umass.edu/roohie/teaching/) spans fluid mechanics, aerodynamics, gas dynamics, propulsion and aerospace structures. In FlowMLLab, I connect physical derivations to reproducible computations so that students can examine not only a model's predictions, but also its assumptions, numerical error and limits of validity.

I received my Ph.D. in Aerospace Engineering from Sharif University of Technology in 2010. My academic path has included Ferdowsi University of Mashhad, Xi'an Jiaotong University, the University of Maryland, Johns Hopkins University and Embry-Riddle Aeronautical University. [Academic background](https://websites.umass.edu/roohie/)

For research or teaching inquiries, please contact **[roohie@umass.edu](mailto:roohie@umass.edu)**.
