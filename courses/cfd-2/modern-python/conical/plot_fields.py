"""Publication-style flow-property figures and reusable dimensional-ratio CSVs."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from reference import sample_reference
from state import primitive


def properties(q, mach, gamma=1.4):
    rho, vr, vt, pressure = q.T
    pr = pressure*gamma*mach**2
    local_mach = np.sqrt(vr**2+vt**2)/np.sqrt(gamma*pressure/rho)
    return np.column_stack([pr, rho, local_mach, pr/rho, vr, vt])


def render_fields(data, references, labels, colors, out, stability=None):
    out = Path(out)
    runs = data["runs"]
    finest = max(r["config"]["cells"] for r in runs)
    names = ["p / p∞", "ρ / ρ∞", "Local Mach number", "T / T∞", "vᵣ / U∞", "vθ / U∞"]
    titles = ["Pressure: compression across the shock", "Density: compressed gas", "Mach: deceleration through the shock",
              "Temperature: compression heating", "Radial velocity: along conical rays", "Tangential velocity: zero at the wall"]
    styles = ["--", "-.", "-", "-", "-", "-"]
    csv_header = "theta_deg,p_over_p_inf,rho_over_rho_inf,mach,T_over_T_inf,vr_over_U_inf,vtheta_over_U_inf"
    continuations = stability["runs"] if stability else []
    for mach, outer in [(2.35, 45), (7.95, 22)]:
        ref = references[mach]
        theta_ref = np.linspace(np.radians(10), np.radians(outer), 2400)
        values_ref = properties(primitive(sample_reference(theta_ref, ref)), mach)
        np.savetxt(out/f"reference-m{mach:g}.csv", np.column_stack([np.degrees(theta_ref), values_ref]),
                   delimiter=",", header=csv_header, comments="")
        profiles = []
        for color, style, (key, label) in zip(colors, styles, labels.items()):
            row = next(r for r in runs if r["config"]["mach"] == mach and r["config"]["cells"] == finest
                       and (r["config"]["flux"], r["config"]["reconstruction"]) == key)
            recovered = next((r for r in continuations if r["parent_run"] == row["file"] and r["converged"]), None)
            if recovered:
                row = recovered
            field = np.load(out/(row["file"]+".npz"))
            theta = np.degrees(field["theta"])
            values = properties(field["primitive"], mach)
            status = " †" if not row["converged"] else ""
            if recovered:
                status += " (CFL 0.15)"
            profiles.append((theta, values, color, style, label+status))
        # CSVs for every grid/method retain sampled cell-center states and all six properties.
        for row in [r for r in runs+continuations if r["config"]["mach"] == mach]:
            field = np.load(out/(row["file"]+".npz"))
            np.savetxt(out/(row["file"]+".csv"), np.column_stack([np.degrees(field["theta"]),
                       properties(field["primitive"], mach)]), delimiter=",", header=csv_header, comments="")
        for zoom in [False, True]:
            fig, axes = plt.subplots(2, 3, figsize=(14, 8))
            for j, ax in enumerate(axes.flat):
                ax.set_facecolor("#f9fbfd")
                for theta, values, color, style, label in profiles:
                    ax.plot(theta, values[:, j], color=color, ls=style, lw=1.65, label=label)
                ax.plot(np.degrees(theta_ref), values_ref[:, j], color="#142332", ls=(0, (4, 3)),
                        lw=2.1, label="Taylor–Maccoll reference")
                ax.axvline(ref["shock_angle_deg"], color="#798a97", ls=":", lw=.9)
                if zoom:
                    half_width = .65 if mach == 2.35 else .3
                    ax.set_xlim(ref["shock_angle_deg"]-half_width, ref["shock_angle_deg"]+half_width)
                    near = abs(np.degrees(theta_ref)-ref["shock_angle_deg"]) < half_width
                    lo, hi = np.min(values_ref[near, j]), np.max(values_ref[near, j])
                    for theta, values, _, _, _ in profiles:
                        selected = abs(theta-ref["shock_angle_deg"]) <= half_width
                        lo = min(lo, np.min(values[selected, j]))
                        hi = max(hi, np.max(values[selected, j]))
                    pad = max(.06*(hi-lo), .005)
                    ax.set_ylim(lo-pad, hi+pad)
                else:
                    ax.set_xlim(10, outer)
                title = titles[j]
                if zoom and j == 4:
                    title = "Radial velocity: continuous across the ideal shock"
                if zoom and j == 5:
                    title = "Tangential velocity: normal to the conical shock"
                ax.set_title(title, fontsize=10, loc="left", fontweight="bold", pad=10)
                ax.set_xlabel("Polar angle θ [deg]")
                ax.set_ylabel(names[j])
                ax.grid(alpha=.18)
                ax.tick_params(labelsize=9)
            handles, legend_labels = axes[0, 0].get_legend_handles_labels()
            fig.legend(handles, legend_labels, loc="lower center", ncol=4, fontsize=9,
                       bbox_to_anchor=(.5, .035), frameon=False)
            fig.suptitle(f"10° cone · M∞ = {mach:g} · {finest} angular cells"
                         +(" · shock detail" if zoom else " · six physical flow properties"),
                         fontsize=17, fontweight="bold", color="#17364f", y=.975)
            foot = "† Not converged: diagnostic field. " if any("†" in p[-1] for p in profiles) else ""
            foot += "CFL 0.15 labels indicate converged checkpoint continuations. " if any("CFL 0.15" in p[-1] for p in profiles) else ""
            fig.text(.5, .012, foot+"Ideal gas, γ=1.4. Cell-center samples; dotted line: reference shock.",
                     ha="center", fontsize=8, color="#536a7c")
            fig.subplots_adjust(left=.065, right=.975, bottom=.17, top=.89, wspace=.28, hspace=.40)
            name = f"flow-properties-m{mach:g}"+("-shock" if zoom else "")
            fig.savefig(out/(name+".png"), dpi=200)
            fig.savefig(out/(name+".pdf"))
            plt.close(fig)


if __name__ == "__main__":
    import json
    from render_report import ROOT, OUT, LABELS, COLORS
    from reference import taylor_maccoll
    data = json.loads((OUT/"comparison.json").read_text(encoding="utf-8"))
    stability = json.loads((OUT/"stability.json").read_text(encoding="utf-8"))
    assert data["complete"] and stability["complete"]
    references = {m:taylor_maccoll(m, max_step=1e-4) for m in [2.35, 7.95]}
    render_fields(data, references, LABELS, COLORS, OUT, stability)
    print("Six-property figures and shock details regenerated without clipping any method's overshoot.")
