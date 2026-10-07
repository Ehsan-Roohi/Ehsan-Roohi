"""Run controlled flux/reconstruction/refinement comparisons and save every run."""
from pathlib import Path
from time import perf_counter
import hashlib
import json
import platform
import sys
import argparse
import numpy as np
from solver import Config, ConicalSolver, save_result
from reference import taylor_maccoll, sample_reference

ROOT = Path(__file__).resolve().parent
METHODS = [("roe", "first"), ("vanleer", "first"), ("hllc", "first"),
           ("hllc", "muscl"), ("hllc", "cweno3"), ("hllc", "cweno5")]


def metrics(result, ref):
    theta, q = result["theta"], result["primitive"]
    p = q[:, 3]*ref["gamma"]*ref["mach"]**2
    uq = sample_reference(theta, ref)
    from state import primitive
    exact_p = primitive(uq, ref["gamma"])[:, 3]*ref["gamma"]*ref["mach"]**2
    h = theta[1]-theta[0]
    gradient = abs(np.diff(p))
    k = int(np.argmax(gradient))
    shift = 0
    if 0 < k < len(gradient)-1:
        denominator = gradient[k-1]-2*gradient[k]+gradient[k+1]
        if denominator != 0:
            shift = np.clip(.5*(gradient[k-1]-gradient[k+1])/denominator, -.5, .5)
    shock = np.degrees(.5*(theta[k]+theta[k+1])+shift*h)
    # Compare only smooth regions for a separate diagnostic; shock error is retained in full L1.
    smooth = abs(theta-np.radians(ref["shock_angle_deg"])) > 4*h
    result.update(wall_pressure_error_percent=float(100*abs(result["wall_pressure_ratio"]/ref["wall_pressure_ratio"]-1)),
                  wall_mach_error_percent=float(100*abs(result["wall_mach"]/ref["wall_mach"]-1)),
                  pressure_ratio_l1=float(np.mean(abs(p-exact_p))),
                  smooth_pressure_ratio_l1=float(np.mean(abs(p[smooth]-exact_p[smooth]))),
                  shock_angle_deg=float(shock), shock_angle_error_deg=float(abs(shock-ref["shock_angle_deg"])))
    return result


def run(cells=(60, 120, 240), max_factor=140, resume=False):
    rows, refs = [], []
    total = 2*len(cells)*len(METHODS)
    count = 0
    for mach, outer in [(2.35, 45), (7.95, 22)]:
        ref = taylor_maccoll(mach, max_step=1e-4)
        refs.append({k: float(v) for k, v in ref.items() if k not in ["theta", "primitive"]})
        for n in cells:
            for flux, reconstruction in METHODS:
                count += 1
                print(f"RUN {count}/{total}: M={mach}, N={n}, {flux}/{reconstruction}", flush=True)
                c = Config(mach=mach, outer_deg=outer, cells=n, flux=flux,
                           reconstruction=reconstruction, max_steps=max_factor*n)
                name = f"m{mach:g}-n{n}-{flux}-{reconstruction}"
                start = perf_counter()
                cached_json = ROOT/"results"/(name+".json")
                cached_npz = ROOT/"results"/(name+".npz")
                use_cache = False
                if resume and cached_json.exists() and cached_npz.exists():
                    cached = json.loads(cached_json.read_text(encoding="utf-8"))
                    from dataclasses import asdict
                    use_cache = cached["config"] == asdict(c)
                if use_cache:
                    field = np.load(cached_npz)
                    result = cached | {k:field[k] for k in ["theta", "primitive", "conserved_average", "history"]}
                    print("  Reusing matching recorded run", flush=True)
                else:
                    result = metrics(ConicalSolver(c).run(), ref)
                save_result(result, ROOT/"results"/name)
                rows.append({k: v for k, v in result.items() if k not in ["theta", "primitive", "conserved_average", "history"]} | {"file": name})
                print(f"  converged={result['converged']} steps={result['iterations']} p/pinf={result['wall_pressure_ratio']:.8f} wall_error={result['wall_pressure_error_percent']:.4f}% shock={result['shock_angle_deg']:.5f} time={perf_counter()-start:.1f}s", flush=True)
                payload = {"research_date": "2026-10-07", "references": refs, "runs": rows,
                    "complete": count == total, "environment": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform(), "precision": "float64", "execution": "sequential CPU; no external CFD package"},
                    "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob("*.py")},
                    "legacy_source_sha256": hashlib.sha256((ROOT.parents[1]/"projects/01-conical-flow/all-codes/Second0.for").read_bytes()).hexdigest(),
                    "interpretation": "Old flux families independently reimplemented on a corrected common conservative FV core; not a successful reproduction of the defective archived executable."}
                (ROOT/"results/comparison.json").write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=int, nargs="+", default=[60, 120, 240])
    parser.add_argument("--max-factor", type=int, default=140)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    run(tuple(args.cells), args.max_factor, args.resume)
