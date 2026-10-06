"""Execute the pinned upstream SIMPLE algorithm and record quantitative checks.

Supply a private checkout with --source-root; external solver source is not vendored.
The uniform rectangle has an exact solution. Airfoil runs are diagnostics, not
reference validation: the upstream wall is a slip wall and no mesh study is done.
"""
from pathlib import Path
import argparse
import contextlib
import hashlib
import io
import json
import sys
import time
import traceback
import warnings

COMMIT = '6a1d90684c9462aa48416894f147910ca20e238f'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', required=True)
    parser.add_argument('--mesh', required=True, help='Relative to upstream checkout')
    parser.add_argument('--case', choices=['uniform', 'airfoil'], required=True)
    parser.add_argument('--iterations', type=int, default=1000)
    parser.add_argument('--budget-seconds', type=float, default=180)
    parser.add_argument('--correctors', type=int, default=1)
    parser.add_argument('--output', required=True, help='Output file prefix')
    a = parser.parse_args()
    source = Path(a.source_root).resolve()
    output = Path(a.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(source))
    import numpy as np
    import scipy
    import gmsh
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from src.mesh import Mesh
    from src.boundaries import process_boundaries
    from src.solution import (initialize_soln, Solution, build_momentum_matrix,
                              assemble_matrices, predicted_velocity,
                              pressure_correction, correct_v_and_p)
    from src.helpers import continuity_residual, save_checkpoint, read_XFOIL_data

    mesh_file = source / a.mesh
    mesh = Mesh(mesh_file)
    if a.case == 'uniform':
        u0, v0, p0 = 1.0, 0.0, 101325.0
    else:
        u0 = float(34.3 * np.cos(np.deg2rad(3)))
        v0 = float(34.3 * np.sin(np.deg2rad(3)))
        p0 = 101325.0
    raw_bc = {'velocity-inlet': {'u': u0, 'v': v0},
              'pressure-outlet': {'p': p0}, 'slip-wall': {}}
    BCs = process_boundaries(raw_bc)
    soln = initialize_soln({'u': u0, 'v': v0, 'p': p0}, mesh)
    urf = {'alpha-u': [0, 1000, 0.4, 0.6], 'alpha-v': [0, 1000, 0.4, 0.6],
           'alpha-p': [0, 1000, 0.2, 0.3]}
    record = {'upstream_commit': COMMIT, 'case': a.case, 'mesh': a.mesh,
              'mesh_sha256': hashlib.sha256(mesh_file.read_bytes()).hexdigest(),
              'versions': {'numpy': np.__version__, 'scipy': scipy.__version__,
                           'gmsh': gmsh.__version__, 'matplotlib': matplotlib.__version__},
              'requested_iterations': a.iterations, 'budget_seconds': a.budget_seconds,
              'non_orthogonal_correctors': a.correctors, 'mu': 1.79e-5, 'rho': 1.0,
              'initial_conditions': {'u': u0, 'v': v0, 'p': p0},
              'mesh_cells': mesh.n_cells, 'mesh_nodes': len(mesh.node_coords),
              'minimum_cell_area': float(np.min(mesh.cell_areas)),
              'all_cell_areas_positive': bool(np.all(mesh.cell_areas > 0)),
              'source_hashes': {f.relative_to(source).as_posix(): hashlib.sha256(f.read_bytes()).hexdigest()
                                for f in sorted((source / 'src').rglob('*.py'))},
              'status': 'running', 'history': [], 'warnings': []}
    hist = []
    start = time.monotonic()
    Au = Av = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            for iteration in range(a.iterations):
                with contextlib.redirect_stdout(io.StringIO()):
                    Mu, Mv, bu, bv = build_momentum_matrix(iteration, mesh, soln, 1.79e-5,
                                                         BCs, Au, Av)
                    Au, Av, Hu, Hv = assemble_matrices(Mu, Mv, soln)
                    us, vs = predicted_velocity(Au, Av, Hu, Hv, bu, bv)
                    pc = pressure_correction(mesh, Au, Av, soln.p, us, vs, BCs, a.correctors)
                    u, v, pressure = correct_v_and_p(iteration, mesh, soln.u, soln.v, soln.p,
                                                   us, vs, Au, Av, pc, urf, BCs)
                    change = float(max(np.max(np.abs(u-soln.u)), np.max(np.abs(v-soln.v))))
                    soln = Solution(u, v, pressure)
                    residual, norms = continuity_residual(mesh, soln, Au, Av, BCs)
                finite = bool(np.isfinite(u).all() and np.isfinite(v).all()
                              and np.isfinite(pressure).all() and np.isfinite(residual).all())
                row = {'iteration': iteration+1, 'continuity_max': float(norms[0]),
                       'continuity_L1': float(norms[1]), 'continuity_net': float(norms[2]),
                       'maximum_velocity_change': change, 'all_fields_finite': finite,
                       'maximum_speed': float(np.max(np.hypot(u, v)))}
                record['history'].append(row)
                hist.append([iteration, *norms])
                if (iteration+1) % 10 == 0:
                    print(f'{a.case} {mesh.n_cells} cells: iteration {iteration+1}, '
                          f'L1={norms[1]:.6g}, speed={row["maximum_speed"]:.6g}', flush=True)
                    output.with_suffix('.json').write_text(json.dumps(record, indent=2)+'\n', 'utf8')
                if not finite:
                    record['status'] = 'non_finite_fields'
                    break
                if row['maximum_speed'] > 1e6 * max(1, abs(u0), abs(v0)):
                    record['status'] = 'unbounded_velocity'
                    break
                if time.monotonic()-start >= a.budget_seconds:
                    record['status'] = 'time_budget_reached'
                    break
            else:
                record['status'] = 'requested_iterations_completed'
            record['completed_iterations'] = len(hist)
            record['all_fields_finite'] = bool(all(h['all_fields_finite'] for h in record['history']))
            record['field_ranges'] = {name: [float(np.min(field)), float(np.max(field))]
                                      for name, field in [('u', soln.u), ('v', soln.v), ('p', soln.p)]}
            if a.case == 'uniform':
                errors = {'u_Linf': float(np.max(np.abs(soln.u-u0))),
                          'v_Linf': float(np.max(np.abs(soln.v-v0))),
                          'p_Linf': float(np.max(np.abs(soln.p-p0)))}
                record['exact_solution_errors'] = errors
                record['uniform_solution_preserved'] = bool(errors['u_Linf'] < 1e-8
                                                           and errors['v_Linf'] < 1e-8
                                                           and errors['p_Linf'] < 1e-6)
            np.savez(output.with_suffix('.npz'), u=soln.u, v=soln.v, p=soln.p)
            np.savetxt(output.with_suffix('.csv'), np.array(hist), delimiter=',',
                       header='iteration_zero_based,continuity_max,continuity_L1,continuity_net', comments='')
            if a.case == 'airfoil' and record['all_fields_finite']:
                with contextlib.redirect_stdout(io.StringIO()):
                    save_checkpoint(str(output.parent / (output.name+'_checkpoint')), len(hist),
                                    mesh, soln, residual, np.array(hist), raw_bc,
                                    read_XFOIL_data(source/'data/xfoil/XFOIL_NACA_2412_AoA3.dat'))
                record['upstream_checkpoint_saved'] = True
            if hist:
                fig, ax = plt.subplots(figsize=(7, 4))
                ax.semilogy([h['iteration'] for h in record['history']],
                            [max(h['continuity_L1'], 1e-30) for h in record['history']])
                ax.set(xlabel='Iteration', ylabel='Continuity L1 (sum of cell flux imbalance)',
                       title=f'Upstream SIMPLE: {a.case}, {mesh.n_cells} cells')
                ax.grid(True, alpha=0.3)
                fig.tight_layout()
                fig.savefig(output.with_suffix('.png'), dpi=160)
                plt.close(fig)
        except Exception:
            record['status'] = 'exception'
            record['exception'] = traceback.format_exc().replace(str(source), '<upstream>')
        record['warnings'] = sorted(set(str(w.message) for w in caught))
    record['elapsed_seconds'] = round(time.monotonic()-start, 3)
    output.with_suffix('.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf8')
    print(json.dumps({k: v for k, v in record.items() if k not in ('history', 'source_hashes')}, indent=2))


if __name__ == '__main__':
    main()
