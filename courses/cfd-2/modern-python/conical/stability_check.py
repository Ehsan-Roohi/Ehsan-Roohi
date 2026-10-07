"""Continue budget-exhausted runs at lower CFL without overwriting the study.

Continuation timings are separate from cold-start timings. Both attempts remain
available, so the common-CFL experiment and its failure are fully reviewable.
"""
from pathlib import Path
from dataclasses import replace
import json
import hashlib
import numpy as np
from solver import Config, ConicalSolver, save_result
from reference import taylor_maccoll
from study import metrics

ROOT = Path(__file__).resolve().parent


def main():
    study = json.loads((ROOT/"results/comparison.json").read_text(encoding="utf-8"))
    if not study["complete"]:
        raise RuntimeError("Complete the controlled study before supplementary continuations")
    results = []
    for row in study["runs"]:
        if row["converged"]:
            continue
        # This is a stability diagnostic from an already advanced state, not a
        # second cold-start grid study. Preserve a completed longer diagnostic.
        c = replace(Config(**row["config"]), cfl=.15, max_steps=10000)
        name = row["file"]+"-cfl015-continuation"
        cached_path = ROOT/"results"/(name+".json")
        if cached_path.exists() and (ROOT/"results"/(name+".npz")).exists():
            cached = json.loads(cached_path.read_text(encoding="utf-8"))
            expected = row["config"] | {"cfl": .15}
            if all(cached["config"][k] == v for k,v in expected.items() if k != "max_steps"):
                results.append(cached | {"file":name})
                print(f"REUSE completed stability diagnostic {name}", flush=True)
                continue
        print(f"CONTINUE {row['file']}: CFL {row['config']['cfl']} -> {c.cfl}", flush=True)
        solver = ConicalSolver(c)
        field = np.load(ROOT/"results"/(row["file"]+".npz"))
        solver.u = field["conserved_average"].copy()
        result = metrics(solver.run(verbose=True), taylor_maccoll(c.mach, max_step=1e-4))
        result.update(parent_run=row["file"], initialization="saved final cell averages of the common-CFL run",
                      parent_iterations=row["iterations"], parent_elapsed_seconds=row["elapsed_seconds"],
                      initial_residual=row["residual_rms_max_equation"])
        save_result(result, ROOT/"results"/name)
        results.append({k:v for k,v in result.items() if k not in ["theta", "primitive", "conserved_average", "history"]} | {"file": name})
        print(f"  converged={result['converged']} residual={result['residual_rms_max_equation']:.3e} wall_error={result['wall_pressure_error_percent']:.5f}%", flush=True)
        payload = {"purpose": "Time-step stability diagnostic; continuation cost is not a cold-start performance comparison",
                   "cfl": .15, "runs": results, "complete": False,
                   "source_sha256": {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob("*.py")}}
        (ROOT/"results/stability.json").write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")
    payload = {"purpose": "Time-step stability diagnostic; original failures retained; continuation timing reported separately; additional budgets are recorded per run, normally 10000 steps",
               "cfl": .15, "runs": results, "complete": True,
               "source_sha256": {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob("*.py")}}
    (ROOT/"results/stability.json").write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
