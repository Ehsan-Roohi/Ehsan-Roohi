"""Independent property and accuracy checks; writes measured verification evidence."""
from pathlib import Path
import json
import numpy as np
from state import conserved, primitive, physical_flux, source
from flux import FLUXES
from reconstruction import coefficients, evaluate, average_matrix
from reference import taylor_maccoll, shock_state
from solver import Config, ConicalSolver, save_result


def checks():
    rng = np.random.default_rng(20261007)
    q = np.column_stack((rng.uniform(.2, 2, 40), rng.uniform(-2, 2, 40),
                         rng.uniform(-2, 2, 40), rng.uniform(.1, 2, 40)))
    u = conserved(q)
    errors = {"state_roundtrip": float(np.max(abs(primitive(u)-q)))}
    for name, flux in FLUXES.items():
        errors[f"{name}_equal_state_flux"] = float(np.max(abs(flux(u, u)-physical_flux(u))))
    # Cell-average polynomial conservation: quadrature integral must recover its input.
    smooth = np.column_stack((1+.2*np.sin(np.arange(40)*.16), np.cos(np.arange(40)*.16)))
    ext = np.pad(smooth, ((3, 3), (0, 0)), mode="wrap")
    for kind in ["first", "muscl", "cweno3", "cweno5"]:
        c = coefficients(ext, kind, .16)
        average = np.einsum("nik,i->nk", c[1:-1], average_matrix([0], c.shape[1]-1)[0])
        errors[f"{kind}_cell_average_preservation"] = float(np.max(abs(average-smooth)))
    convergence = {}
    for kind in ["first", "muscl", "cweno3", "cweno5"]:
        rows = []
        for n in [20, 40, 80, 160]:
            h = 2*np.pi/n
            x = (np.arange(n)+.5)*h
            mean = (1+.2*np.sinc(h/(2*np.pi))*np.sin(x))[:, None]
            c = coefficients(np.pad(mean, ((3, 3), (0, 0)), mode="wrap"), kind, h)
            f = evaluate(c, .5)[:-1, 0]
            numerical = -(f[1:]-f[:-1])/h
            exact = -.2*(np.sin(x+h/2)-np.sin(x-h/2))/h
            error = float(np.sqrt(np.mean((numerical-exact)**2)))
            order = None if not rows else float(np.log2(rows[-1]["l2_error"]/error))
            rows.append({"cells": n, "l2_error": error, "observed_order": order})
        convergence[kind] = rows
    references = []
    for mach in [2.35, 7.95]:
        coarse = taylor_maccoll(mach, max_step=4e-4)
        fine = taylor_maccoll(mach, max_step=1e-4)
        beta = np.radians(fine["shock_angle_deg"])
        vr, vt, rho, p = shock_state(beta, mach, 1.4)
        upstream = conserved(np.array([[1, np.cos(beta), -np.sin(beta), 1/(1.4*mach*mach)]]))
        downstream = conserved(np.array([[rho, vr, vt, p]]))
        fluxjump = float(np.max(abs(physical_flux(upstream)-physical_flux(downstream))))
        references.append({"mach": mach, "shock_angle_deg": fine["shock_angle_deg"],
            "wall_pressure_ratio": float(fine["wall_pressure_ratio"]), "wall_mach": fine["wall_mach"],
            "ode_step_refinement_pressure_difference": float(abs(fine["wall_pressure_ratio"]-coarse["wall_pressure_ratio"])),
            "shock_rankine_hugoniot_flux_error": fluxjump,
            "wall_normal_velocity": float(fine["primitive"][0, 2])})
    # This is the conical (rotating-basis) freestream, not a constant angular state.
    free_residual = {}
    for kind in ["first", "muscl", "cweno3", "cweno5"]:
        rows = []
        for n in [40, 80, 160]:
            s = ConicalSolver(Config(cells=n, mach=2.35, outer_deg=45, reconstruction=kind, boundary="uniform"))
            r = s.rhs(s.u)
            rows.append({"cells": n, "max_residual": float(np.max(abs(r)))})
        free_residual[kind] = rows
    balanced = ConicalSolver(Config(cells=30, mach=2.35, outer_deg=45, reconstruction="cweno5", boundary="uniform", balance_freestream=True))
    errors["optional_equilibrium_subtraction"] = float(np.max(abs(balanced.rhs(balanced.u))))
    # Regression for decimal-Mach run names: replacing a pathlib suffix loses case identity.
    output_check = Path(__file__).parent/"results"/"io-check-m2.35-n1"
    save_result({"theta": np.array([.2]), "primitive": np.array([[1., 1., 0., .1]]),
                 "conserved_average": np.array([[1., 1., 0., .75]]), "history": [[1, 0, 1]],
                 "config": {"case": "decimal-filename-check"}}, output_check)
    io_preserves_name = (Path(str(output_check)+".json").is_file() and Path(str(output_check)+".npz").is_file())
    passed = (io_preserves_name and max(errors.values()) < 1e-11 and
              convergence["cweno3"][-1]["observed_order"] > 2.7 and
              convergence["cweno5"][-1]["observed_order"] > 4.5 and
              all(row["shock_rankine_hugoniot_flux_error"] < 1e-12 and
                  abs(row["wall_normal_velocity"]) < 1e-10 and
                  row["ode_step_refinement_pressure_difference"] < 1e-7 for row in references))
    return {"passed": bool(passed), "decimal_run_filename_regression": io_preserves_name, "property_errors": errors, "scalar_fv_spatial_operator": convergence,
            "taylor_maccoll_checks": references, "conical_uniform_residual_refinement": free_residual,
            "limitations": ["Spatial order measured on smooth scalar FV operator, not pointwise at cone shock.",
                "SSPRK3 time order is three; the steady CWENO5 solver has nominal fifth-order spatial reconstruction.",
                "Conical geometry has a finite discretization residual without optional equilibrium subtraction."]}


if __name__ == "__main__":
    result = checks()
    out = Path(__file__).parent/"results"/"verification.json"
    out.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    if not result["passed"]:
        raise SystemExit(1)
