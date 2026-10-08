# Lectures and lecture notes

The personal CFD 2 projects are identified by Ehsan Roohi as coursework for Dr. Mazaheri. The lecture PDFs located for this publication are supplementary material from other named instructors. Mazaheri's own lecture files have not yet been identified in the available folders.

## 2026 Python-project teaching companion

[Modern CFD Algorithms - English lecture PDF](../modern-python/shock-tube/output/pdf/CFD_Algorithms_Lecture_2026.pdf), **33 pages**. This newly written course companion teaches all 15 implemented Euler fluxes, standalone JST, 20 FV reconstructions, DG degrees 0/1/2, SSPRK3, positivity controls and measured shock-tube refinement. It includes algorithms, advantages/limitations, worked examples and exercises with answers. The final reading map attributes the open LeVeque, Shu and Persson resources used for teaching ideas; their lecture files are linked rather than reproduced. This companion is separate from the preserved historical instructor notes below.


## Kazem Hejranfar - Sharif University of Technology

The title pages of the seven chapter PDFs explicitly identify Kazem Hejranfar and the Aerospace Engineering Department, Sharif University of Technology. Original slides are preserved without editing. Chapter topics below were identified from the opening pages.

| Chapter |Topic | Pages | File |
| --- | --- | ---: | --- |
| 1 |Governing equations and classification | 41 | [PDF](hejranfar/chapter-01.pdf) |
| 2 |Euler equations and mathematical structure | 36 | [PDF](hejranfar/chapter-02.pdf) |
| 3 |Boundary conditions | 14 | [PDF](hejranfar/chapter-03.pdf) |
| 4 |Structured grid generation | 70 | [PDF](hejranfar/chapter-04.pdf) |
| 5 |Finite-volume method | 56 | [PDF](hejranfar/chapter-05.pdf) |
| 6 |Incompressible Navier–Stokes and pressure–velocity coupling | 38 | [PDF](hejranfar/chapter-06.pdf) |
| 7 |Compact finite differences | 18 | [PDF](hejranfar/chapter-07.pdf) |

[Additional scanned chapter 5 notes](hejranfar/chapter-05-scanned.pdf) - 63 pages. This file is in the same Hejranfar source folder; its original filename is `CFD-2-Kazem-chap-5.pdf`. Its content is image-based, so no text-derived topic or authorship beyond the folder attribution is asserted.

## Peter Bakker and Bram van Leer - supplementary gas dynamics material

The source folder is `Van Leer Notes`. Visual inspection of the gas dynamics cover identifies Prof. Peter Bakker (Delft University of Technology) and Prof. Bram van Leer (University of Michigan), AE4-140, dated February 21, 2005. The separate AE 520 Fall 2011 outline identifies Bram van Leer. These are supplementary materials, not Mazaheri's lectures.

| Material | Pages | Original file |
| --- | ---: | --- |
| AE 520 Compressible Flow, Fall 2011 course outline | 3 | [courseOutline11.pdf](van-leer/courseOutline11.pdf) |
| Gas dynamics lecture notes, AE4-140 - Peter Bakker and Bram van Leer (2005) | 241 | [gasdynamics_LectureNotes.pdf](van-leer/gasdynamics_LectureNotes.pdf) |
| Scanned AE 520 exam and 2D supersonic/Burgers notes | 11 | [lecturenotes2DSSBurgers.pdf](van-leer/lecturenotes2DSSBurgers.pdf) |

The original ZIP contained those same three PDFs; the duplicate ZIP is omitted. Gas dynamics text extraction has font-encoding problems and the final PDF is scanned; its first page is an AE 520 Fall 2011 exam, despite the original Burgers filename; use the PDFs for reading.

## Suggested reading alongside the projects

- Chapters 1-3: conservation laws, Euler equations and boundary conditions before the conical-flow and compressible-airfoil projects.
- Chapter 4: mesh generation before the grid project.
- Chapter 5: finite-volume face fluxes before reviewing the airfoil solvers.
- Chapter 6: pressure-velocity coupling before studying the external SIMPLE code.
- Chapter 7: compact differences as an additional numerical-method topic; the archived solvers are not claimed to implement it.

## RWTH Aachen - Modern Simulation Software Development

[Eleven archived lecture PDFs and syllabus](aachen/README.md), Summer 2026. Dr. Lambert Theisen and Dr. Georgii Oblapenko. Includes FEM, time integration, FVM, DG, particle methods, DSMC and uncertainty quantification. These are further supplementary notes found in the local `CFD Achen` ZIP archives.
