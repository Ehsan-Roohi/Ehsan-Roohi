# Computational Fluid Dynamics II - Mazaheri
# دینامیک سیالات محاسباتی ۲ - مظاهری

Ehsan Roohi's historical CFD 2 coursework: conical-flow solvers, compressible and viscous airfoil calculations, grid generation, mesh inputs and project reports. Supplementary lecture notes by Kazem Hejranfar, Peter Bakker, Bram van Leer and RWTH Aachen are indexed with their original attribution.

این مجموعه شامل کدهای فورترن و گزارش‌های پروژه‌های درس CFD 2 مظاهری، ورودی شبکه‌ها، نمونهٔ نتایج قدیمی و جزوه‌های تکمیلی است. نام مدرس جزوه‌ها مستقل از نام درس درج شده است.

**[Lectures and lecture notes / لکچرها و جزوه‌ها](lectures/README.md)** · **[Code and project guide](projects/README.md)** · **[SIMPLE code review / بررسی کد SIMPLE](docs/SIMPLE-review.fa.md)** · **[Running the legacy codes](docs/running-legacy-codes.md)** · **[Provenance and archive policy](docs/provenance.md)**

## Materials

| Section | Contents |
| --- | --- |
| [Lectures](lectures/README.md) | Seven Hejranfar chapter PDFs, one additional scanned chapter, three Van Leer supplementary PDFs, and eleven RWTH Aachen lecture PDFs with a syllabus. Mazaheri's own lecture files have not yet been located. |
| [Project 1: conical flow](projects/01-conical-flow/README.md) | Roe, Van Leer and H-L-R method branches; primitive/conservative variables and inviscid/viscous boundary-condition variants. |
| [Project 2: compressible airfoil](projects/02-compressible-airfoil/README.md) | Legacy two-dimensional airfoil solvers, five preserved mesh cases and matching mesh/boundary files. |
| [Project 3: viscous airfoil](projects/03-viscous-airfoil/README.md) | Personal Roohi-code variants and a mixed-connectivity viscous case. |
| [Grid generation](projects/04-grid-generation/README.md) | Personal airfoil-grid programs and two selected historical grid/convergence outputs. |
| [Reports](reports/README.md) | Five original Word project documents, preserved as archival files. |
| [Historical results](results/README.md) | Four selected original output files, not newly generated validation results. |
| [SIMPLE review](docs/SIMPLE-review.fa.md) | What the Python solver does, how it relates to the course, and reproducible diagnostics for two specific defects. |

## Start here

1. Read the relevant [lecture chapter](lectures/README.md).
2. Follow the corresponding [project guide](projects/README.md).
3. Read the [legacy runtime constraints](docs/running-legacy-codes.md), then prepare an isolated working copy using `tools/prepare_run.py`.
4. Verify the archival files with `python tools/verify_archive.py`.

The original Fortran files are preserved byte for byte. Several contain Windows `CALL SYSTEM("del ...")` cleanup commands. The working-copy helper removes those commands before experimentation. No Fortran compilation or full CFD solution has been validated for this archive. Mesh cases and historical results should be treated as source material for reproduction, not a certified modern solver.

The external [jyoun35/SIMPLE](https://github.com/jyoun35/SIMPLE) repository is linked and reviewed separately. Its source is not included in this archive.
