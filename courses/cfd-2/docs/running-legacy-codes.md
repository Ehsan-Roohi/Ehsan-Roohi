# Running the legacy codes

The archive is a historical teaching collection. All 23 Fortran files have now been compiled in a Linux diagnostic run with documented working-copy adaptations. Both grid generators complete, but the course flow attempts expose runtime defects. Read the [execution and validation report](validation/README.md) before interpreting outputs.

## Prepare a separate working directory

From this course directory:

```text
python tools/prepare_run.py projects/02-compressible-airfoil/cases/85-medium work/airfoil-85-medium
```

The helper copies a selected source directory and removes active `CALL SYSTEM("del ...")` statements only from the copied `.for`/`.f90` sources. It refuses an existing destination. Original archival files remain unchanged.

## Check before compilation

- Legacy fixed-form files use tabs, continuation conventions and long lines; select compiler options after inspecting the source format.
- Several Project 2 variants import Microsoft/Compaq-era `MSIMSL` and `MSFLIB`. Determine which routines are used before replacing libraries or removing an import.
- Mesh dimensions and run parameters are compiled into many variants. Verify `NN`, `NE`, boundary counts, Mach/angle and Reynolds settings against that case's inputs.
- Filenames are case-insensitive on the original Windows environment. On Linux, normalize filenames in the working copy or update file opens consistently.
- Outputs use fixed filenames and can overwrite old data. Run each configuration in its own directory and record compiler version, source hash, mesh hashes and settings.

## Reproduce the recorded diagnostics

Use `python tools/validate_fortran.py --compiler /usr/bin/gfortran --work work/new-run --output work/results.json --timeout 30` on Linux, or use the repository CFD II validation workflow. The tool records each compilation, working-copy edit, staged input, process result and failure log. Original files remain unchanged.

## Validation expected for a modern reproduction

Check positive density/pressure, mesh validity, boundary-condition consistency, mass balance, momentum/energy residuals, and mesh refinement. Compare applicable cases to analytical/reference results. Historical residual plots alone do not establish all of these.

## Unity HPC

For any future Unity execution, discover/allocate the HPC Workspace with `ws_list`/`ws_allocate` and use the returned shared scratch path for live outputs, checkpoints and logs. Check expiry and capacity. Record a scratch-to-project mapping; archive validated results and selected recovery checkpoints to durable project storage. Do not guess a scratch path or use `/project` for frequent live writes. This publication submits no HPC job.
