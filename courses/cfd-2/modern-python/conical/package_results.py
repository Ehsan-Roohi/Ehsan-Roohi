"""Check recorded fields and package source, measured results and figure PDFs."""
from pathlib import Path
import hashlib
import json
import zipfile
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT/"results"


def main():
    study = json.loads((OUT/"comparison.json").read_text(encoding="utf-8"))
    stability = json.loads((OUT/"stability.json").read_text(encoding="utf-8"))
    assert study["complete"] and stability["complete"]
    rows = study["runs"]+stability["runs"]
    checks = []
    for row in rows:
        field = np.load(OUT/(row["file"]+".npz"))
        q, u = field["primitive"], field["conserved_average"]
        csv = np.loadtxt(OUT/(row["file"]+".csv"), delimiter=",", skiprows=1)
        assert np.all(np.isfinite(q)) and np.all(np.isfinite(u)) and np.all(np.isfinite(csv))
        assert np.all(q[:, 0] > 0) and np.all(q[:, 3] > 0)
        assert csv.shape == (row["config"]["cells"], 7)
        assert np.allclose(csv[:, 4], csv[:, 1]/csv[:, 2], atol=1e-13, rtol=1e-13)
        assert row["converged"] == (row["residual_rms_max_equation"] < row["config"]["tolerance"])
        checks.append({"file":row["file"], "finite_positive":True, "csv_shape":list(csv.shape),
                       "convergence_label_consistent":True})
    (OUT/"field-integrity.json").write_text(json.dumps({"passed":True,"recorded_runs_checked":len(checks),"checks":checks}, indent=2)+"\n", encoding="utf-8")
    current = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob("*.py")}
    provenance = {"study_execution_source_sha256":study["source_sha256"], "delivered_source_sha256":current,
                  "post_study_changes":"Solver defaults changed to the tested HLLC/CWENO3 combination and 40000 steps; mathematical implementations and explicitly configured common-study runs are unchanged. Report annotations/figures and supplementary diagnostic/packaging utilities were finalized afterward.",
                  "supplementary_budget_note":"Additional iterations are explicit per continuation: one completed 57600-step diagnostic is preserved; other completed stability diagnostics use 10000 additional steps. Intermediate interrupted attempts are not accepted numerical results."}
    (OUT/"source-provenance.json").write_text(json.dumps(provenance, indent=2)+"\n", encoding="utf-8")
    files = list(ROOT.glob("*.py"))+[ROOT/"README.md"]
    for row in rows:
        files.extend(OUT/(row["file"]+suffix) for suffix in [".json", ".npz", ".csv"])
    for name in ["comparison.json", "verification.json", "stability.json", "summary.json", "report.html",
                 "field-integrity.json", "source-provenance.json", "reference-m2.35.csv", "reference-m7.95.csv"]:
        files.append(OUT/name)
    files += list(OUT.glob("*.png"))+list(OUT.glob("*.pdf"))
    files.append(ROOT.parents[1]/"projects/01-conical-flow/all-codes/Second0.for")
    manifest = {str(p.relative_to(ROOT.parents[2])).replace("\\", "/"):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    archive = ROOT.parents[2]/"conical-python-and-results.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in files:
            z.write(p, str(p.relative_to(ROOT.parents[2])))
        z.writestr("package-manifest.json", json.dumps(manifest, indent=2)+"\n")
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
    print(json.dumps({"recorded_fields_checked":len(checks), "package":str(archive), "files":len(files)+1,
                      "bytes":archive.stat().st_size, "sha256":hashlib.sha256(archive.read_bytes()).hexdigest()},indent=2))


if __name__ == "__main__":
    main()
