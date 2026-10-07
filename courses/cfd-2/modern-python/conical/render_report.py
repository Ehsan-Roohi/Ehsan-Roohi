"""Render English evidence tables and scientific figures from measured run files."""
from pathlib import Path
import html
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reference import taylor_maccoll, sample_reference
from state import primitive
from plot_fields import render_fields

ROOT = Path(__file__).resolve().parent
OUT = ROOT/"results"
LABELS = {("roe", "first"): "Roe / first order", ("vanleer", "first"): "Van Leer / first order",
          ("hllc", "first"): "HLLC / first order", ("hllc", "muscl"): "HLLC / MUSCL-MC",
          ("hllc", "cweno3"): "HLLC / CWENO3", ("hllc", "cweno5"): "HLLC / CWENO5"}
COLORS = ["#778899", "#8b5a9b", "#d56c33", "#19998c", "#347bc5", "#c6355c"]
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                    "figure.dpi": 120, "savefig.dpi": 180, "pdf.fonttype": 42, "ps.fonttype": 42})


def main():
    data = json.loads((OUT/"comparison.json").read_text(encoding="utf-8"))
    verification = json.loads((OUT/"verification.json").read_text(encoding="utf-8"))
    stability = json.loads((OUT/"stability.json").read_text(encoding="utf-8")) if (OUT/"stability.json").exists() else None
    if stability and not stability["complete"]:
        raise RuntimeError("Supplementary stability study still running")
    runs = data["runs"]
    finest = max(r["config"]["cells"] for r in runs)
    if not data["complete"]:
        raise RuntimeError("Comparison still running; refusing to label it complete")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    references = {}
    for row, mach in enumerate([2.35, 7.95]):
        ref = taylor_maccoll(mach, max_step=1e-4)
        references[mach] = ref
        outer = 45 if mach == 2.35 else 22
        th = np.linspace(np.radians(10), np.radians(outer), 1800)
        qref = primitive(sample_reference(th, ref))
        pref = qref[:, 3]*1.4*mach*mach
        axes[row, 0].plot(np.degrees(th), pref, "k--", lw=1.8, label="Taylor–Maccoll")
        for color, (key, label) in zip(COLORS, LABELS.items()):
            r = next(r for r in runs if r["config"]["mach"] == mach and r["config"]["cells"] == finest and (r["config"]["flux"], r["config"]["reconstruction"]) == key)
            field = np.load(OUT/(r["file"]+".npz"))
            theta = np.degrees(field["theta"])
            pressure = field["primitive"][:, 3]*1.4*mach*mach
            axes[row, 0].plot(theta, pressure, color=color, lw=1.2, label=label+(" †" if not r["converged"] else ""))
            axes[row, 1].plot(theta, pressure, color=color, lw=1.5)
        axes[row, 1].plot(np.degrees(th), pref, "k--", lw=1.8)
        beta = ref["shock_angle_deg"]
        axes[row, 1].set_xlim(beta-.6, beta+.6)
        pshock = ref["primitive"][-1, 3]*1.4*mach*mach
        axes[row, 1].set_ylim(.98, pshock*1.08)
        axes[row, 0].set_title(f"M∞={mach:g}: full pressure profile, N={finest}")
        axes[row, 1].set_title(f"M∞={mach:g}: resolved shock neighborhood")
        for ax in axes[row]:
            ax.set_xlabel("Polar angle θ [deg]")
            ax.set_ylabel("p / p∞")
            ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8, ncol=2)
    fig.suptitle("Pressure profiles and shock detail · † marks a run that missed the residual gate", fontsize=12)
    fig.savefig(OUT/"pressure-comparison.png")
    fig.savefig(OUT/"pressure-comparison.pdf")
    plt.close(fig)
    render_fields(data, references, LABELS, COLORS, OUT, stability)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    for ax, mach in zip(axes, [2.35, 7.95]):
        for color, (key, label) in zip(COLORS, LABELS.items()):
            r = next(r for r in runs if r["config"]["mach"] == mach and r["config"]["cells"] == finest and (r["config"]["flux"], r["config"]["reconstruction"]) == key)
            field = np.load(OUT/(r["file"]+".npz"))
            history = field["history"]
            ax.semilogy(history[:, 0], history[:, 1], color=color, lw=1.3,
                        label=label+(" (budget exhausted)" if not r["converged"] else ""))
        ax.axhline(2e-8, color="black", ls="--", lw=1, label="Residual tolerance")
        ax.set(title=f"M∞={mach:g}: convergence, N={finest}", xlabel="SSPRK3 iterations", ylabel="Maximum equation RMS residual")
        ax.grid(alpha=.2)
        ax.legend(fontsize=7)
    fig.savefig(OUT/"convergence.png")
    fig.savefig(OUT/"convergence.pdf")
    plt.close(fig)
    if stability and stability["runs"]:
        continued = stability["runs"]
        nrows = (len(continued)+1)//2
        fig, axes = plt.subplots(nrows, 2, figsize=(12, 3.5*nrows), constrained_layout=True, squeeze=False)
        for ax, child in zip(axes.flat, continued):
            parent = next(r for r in runs if r["file"] == child["parent_run"])
            a = np.load(OUT/(parent["file"]+".npz"))["history"]
            b = np.load(OUT/(child["file"]+".npz"))["history"]
            ax.semilogy(a[:, 0], a[:, 1], color="#c96b35", lw=1.3, label="Original CFL 0.35")
            ax.semilogy(parent["iterations"]+b[:, 0], b[:, 1], color="#167b87", lw=1.4, label="Continuation CFL 0.15")
            ax.axvline(parent["iterations"], color="#798a97", ls=":", lw=1)
            ax.axhline(2e-8, color="black", ls="--", lw=1, label="Residual tolerance")
            c = child["config"]
            ax.set(title=f"M∞={c['mach']:g}, N={c['cells']}: {LABELS[(c['flux'],c['reconstruction'])]}",
                   xlabel="Total iterations including original attempt", ylabel="Maximum equation RMS residual")
            ax.grid(alpha=.2)
            ax.legend(fontsize=8)
        for ax in list(axes.flat)[len(continued):]:
            ax.set_visible(False)
        fig.savefig(OUT/"stability-continuation.png")
        fig.savefig(OUT/"stability-continuation.pdf")
        plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    for color, (key, label) in zip(COLORS, LABELS.items()):
        subset = sorted([r for r in runs if r["config"]["mach"] == 7.95 and (r["config"]["flux"], r["config"]["reconstruction"]) == key], key=lambda r:r["config"]["cells"])
        axes[0].loglog([r["config"]["cells"] for r in subset], [r["wall_pressure_error_percent"] for r in subset], "o-", color=color, label=label)
        axes[1].loglog([r["elapsed_seconds"] for r in subset], [r["pressure_ratio_l1"] for r in subset], "o-", color=color)
        unfinished = [r for r in subset if not r["converged"]]
        if unfinished:
            axes[0].scatter([r["config"]["cells"] for r in unfinished], [r["wall_pressure_error_percent"] for r in unfinished], marker="x", s=70, color="black", zorder=4)
            axes[1].scatter([r["elapsed_seconds"] for r in unfinished], [r["pressure_ratio_l1"] for r in unfinished], marker="x", s=70, color="black", zorder=4)
    axes[0].set(xlabel="Angular cells", ylabel="Wall pressure error [%]", title="M∞=7.95: refinement")
    axes[1].set(xlabel="Measured solve time [s]", ylabel="Full-profile L1 error in p/p∞", title="M∞=7.95: error versus cost")
    for kind, color in [("first", COLORS[2]), ("muscl", COLORS[3]), ("cweno3", COLORS[4]), ("cweno5", COLORS[5])]:
        rows = verification["scalar_fv_spatial_operator"][kind]
        order = rows[-1]["observed_order"]
        axes[2].loglog([r["cells"] for r in rows], [r["l2_error"] for r in rows], "o-", color=color, label=f"{kind}: p={order:.2f}")
    axes[2].set(xlabel="Cells", ylabel="Spatial RHS L2 error", title="Smooth scalar FV verification")
    axes[0].legend(fontsize=7)
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2, which="both")
    fig.suptitle("Accuracy and cost · black × marks a run that missed the residual gate", fontsize=12)
    fig.savefig(OUT/"accuracy-cost.png")
    fig.savefig(OUT/"accuracy-cost.pdf")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    x, y = np.meshgrid(np.linspace(.01, 1, 400), np.linspace(0, .41, 220))
    theta = np.arctan2(y, x)
    ref = references[7.95]
    for ax, key in zip(axes, [("roe", "first"), ("hllc", "cweno3")]):
        r = next(r for r in runs if r["config"]["mach"] == 7.95 and r["config"]["cells"] == finest and (r["config"]["flux"], r["config"]["reconstruction"]) == key)
        field = np.load(OUT/(r["file"]+".npz"))
        q = field["primitive"]
        pressure = np.interp(theta, field["theta"], q[:, 3]*1.4*7.95**2, right=1)
        pressure = np.ma.masked_where(theta < np.radians(10), pressure)
        im = ax.pcolormesh(x, y, pressure, shading="auto", cmap="turbo", vmin=1, vmax=ref["wall_pressure_ratio"]*1.015, rasterized=True)
        ax.fill_between([0, 1], [0, 0], [0, np.tan(np.radians(10))], color="#d0d6dc")
        ax.plot([0, 1], [0, np.tan(np.radians(ref["shock_angle_deg"]))], "w--", lw=1, label="Reference shock")
        ax.set(xlabel="Axial coordinate x / L", ylabel="Radius r / L", title=f"{LABELS[key]}, M∞=7.95, N={finest}", xlim=(0, 1), ylim=(0, .41))
        ax.set_aspect("equal")
    fig.colorbar(im, ax=axes, label="p / p∞", shrink=.85)
    fig.suptitle("Meridional reconstruction of the angular conical solution (not a separate 2D solve)", fontsize=11)
    fig.savefig(OUT/"cone-contours.png")
    fig.savefig(OUT/"cone-contours.pdf")
    plt.close(fig)
    def table_rows(rows):
        text = []
        for r in rows:
            c = r["config"]
            label = LABELS[(c["flux"], c["reconstruction"])]
            status = "PASS" if r["converged"] else "NOT CONVERGED"
            text.append(f"<tr><td>{c['mach']:g}</td><td>{c['cells']}</td><td>{label}</td><td>{r['wall_pressure_ratio']:.7f}</td><td>{r['wall_pressure_error_percent']:.4f}%</td><td>{r['shock_angle_error_deg']:.4f}°</td><td>{r['pressure_ratio_l1']:.5g}</td><td>{r['elapsed_seconds']:.2f}</td><td>{status}</td></tr>")
        return "".join(text)
    header = "<tr><th>M∞</th><th>N</th><th>Method</th><th>Wall p/p∞</th><th>Wall error</th><th>Shock-angle error</th><th>Profile L1</th><th>Time [s]</th><th>Residual gate</th></tr>"
    highest = [r for r in runs if r["config"]["cells"] == finest]
    failed = sum(not r["converged"] for r in runs)
    all_attempts = runs+(stability["runs"] if stability else [])
    counts = {"runs": len(runs), "converged": len(runs)-failed, "finest_cells": finest,
              "max_balance_roundoff": max(r["discrete_balance_roundoff"] for r in all_attempts),
              "min_density": min(r["min_density"] for r in all_attempts), "min_pressure": min(r["min_pressure"] for r in all_attempts)}
    counts["lower_cfl_continuations"] = len(stability["runs"]) if stability else 0
    counts["lower_cfl_converged"] = sum(r["converged"] for r in stability["runs"]) if stability else 0
    (OUT/"summary.json").write_text(json.dumps(counts, indent=2)+"\n", encoding="utf-8")
    stability_section = ""
    if stability and stability["runs"]:
        stability_section = f'''<section><h2>Lower-CFL stability checks</h2><p>{counts['lower_cfl_converged']} of {counts['lower_cfl_continuations']} checkpoint continuations reached the same residual tolerance at CFL=0.15. Original CFL=0.35 failures remain in the controlled-study table and convergence plot. This is a separate diagnostic, not a cold-start performance comparison: reported continuation times exclude the original attempt.</p><img src="stability-continuation.png" alt="Original residual histories and lower CFL checkpoint continuations"><div class="table"><table>{header}{table_rows(stability['runs'])}</table></div><p>Six-property figures use the converged lower-CFL continuation where available, identified in the legend. The original pressure-only and cost/error figures retain the common-CFL experiment. See <a href="stability.json">stability.json</a> for parent files, starting residuals, iteration counts and separate costs.</p></section>'''
    page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Conical CFD — Python modernization and measured comparison</title>
<style>body{{font:16px/1.65 system-ui,sans-serif;margin:0;color:#193247;background:#f2f6f9}}header{{background:#102e4c;color:white;padding:40px max(25px,calc((100vw - 1120px)/2));border-bottom:6px solid #f2c550}}h1{{line-height:1.2;font-size:36px}}main{{max-width:1180px;margin:auto;padding:25px}}section{{background:white;padding:28px;margin:22px 0;border:1px solid #dbe5ed;border-radius:10px}}h2{{color:#102e4c;line-height:1.3}}img{{max-width:100%;height:auto}}a{{color:#07689c}}.notice{{background:#fff6d8;border-left:4px solid #edbc2f;padding:16px}}.table{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}th{{background:#17364f;color:white;text-align:left}}td,th{{padding:10px;border-bottom:1px solid #dbe5ed;white-space:nowrap}}tr:nth-child(even){{background:#f4f8fb}}code,pre{{font-family:Consolas,monospace}}pre{{background:#eef3f7;padding:18px;white-space:pre-wrap}}.cards{{display:flex;gap:15px;flex-wrap:wrap}}.cards p{{flex:1;background:#eef5fa;padding:18px;min-width:150px}}.cards b{{display:block;font-size:28px}}small{{color:#587084}}</style></head><body>
<header><small style="color:#d5e5ef">Measured CPU runs · October 7, 2026 · English reproducible evidence</small><h1>Conical flow: from the archived Fortran equations to modular Python</h1><p>Roe and Van Leer baselines, HLLC, limited MUSCL, and third-/fifth-order central WENO finite volumes.</p></header><main>
<section><h2>What was completed</h2><p>The inviscid angular conical Euler equations were independently reimplemented in Python/NumPy, preserving the historical radial/tangential variables and geometric source terms. Flux, reconstruction, source quadrature, boundaries, time integration and diagnostics are separate modules. Original Fortran files were preserved.</p>
<div class="cards"><p><b>{len(runs)}</b>controlled flow runs</p><p><b>{len(runs)-failed}/{len(runs)}</b>reached residual tolerance</p><p><b>3 &amp; 5</b>verified smooth spatial orders</p><p><b>CPU</b>no external CFD solver</p></div>
<p>Cases: Mach 7.95 / 10° cone / 22° outer boundary (historical-style inviscid wall recovery), and Mach 2.35 / 10° cone / 45° outer boundary (independent benchmark). Both use gamma=1.4, rho∞=1, U∞=1 and p∞=1/(gamma M∞²). Angular grids have 60, 120 and 240 physical cells plus ghosts; archived N=240 included boundary entries, so those grids are comparable in scale, not identical.</p>
<div class="notice"><strong>Comparison meaning:</strong> classical flux families are reimplemented on the same corrected conservative core. The archived executable has documented runtime and second-order defects and is not used as a numerical truth baseline. This phase covers its <strong>inviscid cone branch</strong>; the archived viscous-cone/heat-transfer branch has not been migrated or validated.</div></section>
<section><h2>Finest-grid measured results</h2><p>Same physical case, cell count, CFL=0.35, SSPRK3 and residual tolerance 2e−8. Flux-family effects are compared at first order; reconstruction effects hold HLLC fixed. Times are sequential CPU measurements on this host, not portable performance guarantees.</p><div class="table"><table>{header}{table_rows(highest)}</table></div>
<p>Independent Taylor–Maccoll references: M2.35 wall pressure ratio {references[2.35]['wall_pressure_ratio']:.8f}, shock {references[2.35]['shock_angle_deg']:.8f}°; M7.95 wall pressure ratio {references[7.95]['wall_pressure_ratio']:.8f}, shock {references[7.95]['shock_angle_deg']:.8f}°.</p><p>Higher reconstruction order does not guarantee the smallest wall error at every grid spacing. Shock alignment, limiter behavior, boundary closure and metric choice all matter. Assess full-profile error, shock location and cost alongside surface pressure.</p></section>
<section><h2>Pressure profiles and shock resolution</h2><img src="pressure-comparison.png" alt="Measured cone pressure profiles and shock neighborhoods for six numerical methods"><p>Shock location is estimated from the strongest pressure gradient with a local parabolic peak fit. Full-profile L1 includes the shock; the formal high-order claim concerns smooth spatial reconstruction, not pointwise convergence through that discontinuity.</p></section>
<section><h2>Six physical flow properties</h2><p>Pressure and density increase across the compression shock, temperature rises, and the local Mach number decreases. The tangential velocity approaches zero at the cone wall. Temperature is inferred from the ideal-gas relation T/T∞=(p/p∞)/(ρ/ρ∞); local Mach uses both velocity components. These are numerical ideal-gas predictions, not measured experimental data.</p><img src="flow-properties-m2.35.png" alt="Pressure density Mach temperature radial and tangential velocity comparisons at Mach 2.35"><img src="flow-properties-m7.95.png" alt="Six physical flow properties compared at Mach 7.95"><details><summary>Resolve differences close to the compression shock</summary><img src="flow-properties-m2.35-shock.png" alt="Mach 2.35 shock detail in six properties"><img src="flow-properties-m7.95-shock.png" alt="Mach 7.95 shock detail in six properties"></details><p><a href="flow-properties-m2.35.pdf">Mach 2.35 figure PDF</a> · <a href="flow-properties-m7.95.pdf">Mach 7.95 figure PDF</a> · <a href="reference-m2.35.csv">Mach 2.35 reference CSV</a> · <a href="reference-m7.95.csv">Mach 7.95 reference CSV</a>. Every recorded method/grid also has a CSV with the same seven columns alongside its JSON and NPZ.</p></section>
<section><h2>Refinement, accuracy and computational cost</h2><img src="accuracy-cost.png" alt="Grid sensitivity, measured runtime versus error, and independent smooth spatial convergence"><p>The smooth scalar finite-volume spatial operator gives observed order {verification['scalar_fv_spatial_operator']['cweno3'][-1]['observed_order']:.3f} for CWENO3 and {verification['scalar_fv_spatial_operator']['cweno5'][-1]['observed_order']:.3f} for CWENO5. MUSCL-MC gives about 1.5 in this global derivative norm because the limiter reduces accuracy near smooth extrema. SSPRK3 remains third order in time; CWENO5 is not a claim of fifth-order transient time integration.</p><p>The cost/error plot includes every recorded run; consult its convergence status in the table before treating an endpoint as a converged steady solution.</p></section>
<section><h2>Actual convergence histories</h2><img src="convergence.png" alt="All-equation residual histories on the finest grid"><p>The stopping gate is the maximum RMS residual over all four conserved equations, rather than wall-pressure stability alone. A run marked NOT CONVERGED is retained as a measured diagnostic and is not accepted as a converged steady result.</p></section>
{stability_section}
<section><h2>Color contour comparison</h2><img src="cone-contours.png" alt="Pressure contours reconstructed from converged Roe first-order and HLLC CWENO3 angular cone solutions"><p>Both displayed fields reached the residual gate. The 1D angular solution is mapped onto conical rays in a meridional plane. This illustrates shock thickness and pressure distribution; it is not a separate two-dimensional simulation. The dashed ray is the independently calculated Taylor–Maccoll shock.</p></section>
<section><h2>Verification and safeguards</h2><ul><li>Primitive/conserved conversion and all five equal-state flux checks pass near roundoff.</li><li>Polynomial reconstructions preserve their input cell averages.</li><li>Three-point Gaussian source quadrature uses the same cell polynomial as the face states; the spherical geometric terms are retained.</li><li>Taylor–Maccoll shooting uses independent RK4/bisection; ODE step refinement changes wall pressure by less than 1e−7. Shock jump fluxes and zero wall-normal velocity are separately checked.</li><li>All run states remained finite with positive density/pressure; minimum reported density {counts['min_density']:.5g}, pressure {counts['min_pressure']:.5g}.</li><li>Global discrete flux/source balance closes within {counts['max_balance_roundoff']:.2g}. Mass alone is not constant in this angular balance law because of radial/geometric source terms.</li><li>Sampled reconstructed states are scaled toward the cell mean when needed; an inadmissible RK update is rejected and retried with smaller dt. HLLC has a face-local HLLE fallback. These checks are not a general entropy-stability or positivity theorem.</li></ul><p>Roe uses acoustic eigenvalue smoothing. Wall ghosts reflect angular momentum; the wall face has zero mass/energy flux and a pressure force. This boundary closure can limit global order. Uniform Cartesian freestream has a rotating angular basis; its uncorrected discrete residual decreases with refinement. Optional discrete equilibrium subtraction exists but is disabled in every reported comparison.</p></section>
<section><h2>A reference discrepancy found during verification</h2><p>The <a href="https://www.grc.nasa.gov/www/wind/valid/cone10/cone10.html">NASA Mach2.35 cone page</a> prints a theoretical wall pressure ratio of 1.4234 and shock angle 27.1843°. Our independently integrated Taylor–Maccoll solution gives 1.37393637 and 26.73671772°. NASA's own listed Wind-US wall pressure, about 1.3740, agrees closely with our independent pressure value. The printed theoretical table/angle are therefore not used as acceptance truth here. The discrepancy is recorded rather than silently fitting the solver to those numbers.</p><p>Check <a href="verification.json">verification.json</a> for the ODE refinement, shock Rankine–Hugoniot and wall-condition evidence. These establish numerical consistency of the reference; they are not experimental validation of a real hypersonic flow.</p></section>
<section><h2>Legacy defects corrected in the new core</h2><p>In <code>Second0.for</code>, the interactive input unit is reopened as an output file; the second-order tangential left flux uses right-face dissipation coefficients; radial and energy second-order branches overwrite <code>FINVT</code>; energy dissipation also uses radial coefficients. The new implementation assembles all four flux components once per face, applies a limiter, and uses one residual per conserved variable. No original archived source was overwritten.</p><p>The original Roe average and unlimited extrapolation are not copied bit for bit. This is a documented physics/algorithm migration with corrected standard implementations. CWENO belongs to the post-archive literature (2016 preprint / 2018 publication), while HLLC, MUSCL and SSPRK are established methods, not inventions of 2026.</p></section>
<section><h2>Reproduce and inspect</h2><pre>python verify.py
python study.py --cells 60 120 240
python stability_check.py
python render_report.py

python solver.py --mach 7.95 --outer 22 --cells 240 \
  --flux hllc --reconstruction cweno3 --steps 40000 --output results/my-run</pre><p>Use Python with NumPy and Matplotlib. Running these commands from this directory on an ordinary machine is supported; in this restricted host runs were launched from the authorized Courses root with full paths. The default is HLLC/CWENO3, which reached the residual gate in all six case/grid combinations. MUSCL-MC and CWENO5 have recorded shock-case convergence limitations despite their good pressure errors in some runs.</p><p><a href="comparison.json">Complete run configuration and measured metrics</a> · <a href="verification.json">Verification evidence</a> · <a href="summary.json">Compact summary</a> · <a href="../README.md">Code guide and equations</a></p><details><summary>All {len(runs)} measured runs</summary><div class="table"><table>{header}{table_rows(runs)}</table></div></details></section>
<section><h2>Method references and limits</h2><p><a href="https://arxiv.org/abs/1607.07319">Cravero, Puppo, Semplice and Visconti: CWENO for balance laws</a> supplies the mathematical framework for whole-cell nonlinear polynomials and source quadrature. <a href="https://arxiv.org/abs/1503.00736">Cravero and Semplice</a> discuss epsilon scaling with mesh size. Our code uses cell-average Vandermonde reconstruction, Jiang–Shu derivative-integral indicators, positive central weights and epsilon proportional to h²; it is not a copied upstream implementation.</p><p>The study establishes an inviscid steady cone benchmark and smooth operator accuracy. It does not yet establish viscous-cone accuracy, unstructured airfoil performance, full high-order boundary closure, physical validation or an entropy-stability proof.</p></section></main></body></html>'''
    (OUT/"report.html").write_text(page, encoding="utf-8")
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
