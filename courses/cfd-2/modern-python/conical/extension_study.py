"""Reproducible controlled AUSM/entropy-flux Euler and viscous comparisons."""
from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import hashlib
import json
import sys
import platform
import subprocess
import numpy as np
from solver import Config, ConicalSolver, save_result
from viscous import ViscousConfig, ViscousConicalSolver
from reference import taylor_maccoll, sample_reference
from study import metrics

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results-v2'
METHODS = ['roe', 'hllc', 'ausm', 'ausm-plus', 'ausm-up', 'ausm-up2', 'slau2', 'ec-lf']


def execute(config, steady=False, warm=False):
    if isinstance(config, ViscousConfig) and config.reconstruction == 'minmod' and (OUT/'stop-diagnostics.flag').exists():
        return {'config':asdict(config), 'file':'skipped-'+config.thermal_wall+'-'+config.flux,
                'converged':False, 'skipped':True,
                'reason':'Superseded by corrected quadratic-wall, primitive reconstruction steady comparison'}
    loaded_sources = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}
    viscous = isinstance(config, ViscousConfig)
    mode = config.thermal_wall if viscous else 'euler'
    name = f'{mode}-m{config.mach:g}-n{config.cells}-{config.flux}-{config.reconstruction}'
    if viscous:
        name += f'-g{config.gradient_order}'
    if steady:
        name = 'steady-'+name
    if getattr(config, 'viscous_operator_version', ''):
        name = 'fv2-'+name
    path = OUT/name
    if Path(str(path)+'.json').exists():
        cached = json.loads(Path(str(path)+'.json').read_text(encoding='utf-8'))
        if cached['config'] == asdict(config) and Path(str(path)+'.npz').exists():
            return cached | {'file': name}
    solver = ViscousConicalSolver(config) if viscous else ConicalSolver(config)
    ref = taylor_maccoll(config.mach, max_step=1e-4)
    if viscous:
        # Identical inviscid reference initialization for every viscous flux.
        # Integrate over each cell; this is only an initial guess, not a viscous target.
        x, w = np.polynomial.legendre.leggauss(8)
        points = solver.theta[:, None]+solver.h[:, None]*x/2
        solver.u = np.einsum('nqk,q->nk', sample_reference(points.ravel(), ref).reshape(config.cells, 8, 4), w/2)
    seed_hash = None
    if warm:
        seed_path = OUT/('warm-adiabatic.npz' if config.thermal_wall == 'adiabatic' else 'seed-isothermal.npz')
        data = np.load(seed_path)
        solver.u = np.column_stack([np.interp(solver.theta,data['theta'],data['conserved_average'][:,k]) for k in range(4)])
        seed_hash = hashlib.sha256(seed_path.read_bytes()).hexdigest()
    try:
        if steady:
            from steady import solve
            result = solve(solver, max_iterations=180, verbose=True)
        else:
            result = solver.run()
        if not viscous:
            result = metrics(result, ref)
        save_result(result, path)
        summary = {k:v for k,v in result.items() if k not in ('theta','primitive','conserved_average','history')}
        return summary | {'file': name, 'source_sha256_at_worker_start': loaded_sources,
                          'initialization': 'common wall-matched coarse steady seed, interpolated conservative states' if warm else ('Taylor-Maccoll cell averages' if viscous else 'freestream cell averages'),
                          'initial_seed_sha256':seed_hash}
    except (ArithmeticError, RuntimeError, FloatingPointError) as exc:
        np.savez_compressed(Path(str(path)+'-failure.npz'), conserved_average=solver.u, theta=solver.theta)
        return {'file': name, 'config': asdict(config), 'converged': False, 'failure': str(exc),
                'events': solver.events, 'initialization': 'Taylor-Maccoll cell averages' if viscous else 'freestream cell averages'}


def run(mode, workers, grids):
    OUT.mkdir(exist_ok=True)
    jobs = []
    if mode == 'euler':
        for mach, outer, ns in [(2.35, 45, grids[:1]), (7.95, 22, grids)]:
            for n in ns:
                for method in METHODS:
                    jobs.append(Config(mach=mach, outer_deg=outer, cells=n, flux=method, max_steps=300*n))
    elif mode == 'viscous':
        for wall in ['adiabatic', 'isothermal']:
            for method in METHODS:
                jobs.append(ViscousConfig(cells=grids[0], flux=method, thermal_wall=wall,
                                          reconstruction='minmod', max_steps=45000))
    elif mode == 'refine':
        for wall in ['adiabatic', 'isothermal']:
            for n in grids:
                for method in ['hllc', 'ausm-up2']:
                    jobs.append(ViscousConfig(cells=n, flux=method, thermal_wall=wall,
                                              reconstruction='minmod', max_steps=150000))
    elif mode == 'high-order':
        for method in ['hllc', 'ausm-up2', 'slau2']:
            jobs.append(ViscousConfig(cells=grids[0], flux=method, gradient_order=4, reconstruction='cweno5-char', max_steps=20000))
    elif mode == 'characteristic':
        for method in ['ausm-up2', 'slau2', 'ec-lf']:
            jobs.append(Config(cells=grids[0], flux=method, reconstruction='cweno3-char', max_steps=20000))
        for method in ['hllc', 'ausm-up2']:
            jobs.append(ViscousConfig(cells=grids[0], flux=method, reconstruction='cweno3-char', max_steps=30000))
    elif mode == 'primitive':
        for wall in ['adiabatic', 'isothermal']:
            for method in ['hllc','ausm-up2','slau2']:
                jobs.append(ViscousConfig(cells=grids[0], flux=method, reconstruction='primitive-minmod', thermal_wall=wall, max_steps=35000))
    elif mode in ('steady', 'steady-warm', 'steady-refine', 'steady-high-order'):
        methods = METHODS if mode in ('steady','steady-warm') else ['hllc','ausm-up2']
        for wall in ['adiabatic', 'isothermal']:
            for n in grids:
                for method in methods:
                    jobs.append(ViscousConfig(cells=n, flux=method, thermal_wall=wall,
                                              reconstruction='cweno5-char' if mode == 'steady-high-order' else 'primitive-minmod',
                                              gradient_order=4 if mode == 'steady-high-order' else 2))
    rows = []
    payload = {'mode': mode, 'complete': False, 'planned_runs': len(jobs), 'runs': rows,
               'date': '2026-10-07', 'residual_definition': 'Maximum over equations of RMS(dU/dpseudo_time), absolute, measured every 50 steps; tolerance 2e-8',
               'sources_sha256': {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')},
               'environment': {'python':sys.version, 'numpy':np.__version__, 'platform':platform.platform(), 'workers':workers},
               'model_scope': 'Viscous runs retain a local radial-station angular approximation; not full axisymmetric NS validation.'}
    def child(index, c):
        job = {'viscous': isinstance(c, ViscousConfig), 'config': asdict(c), 'steady': mode.startswith('steady'),
               'warm': mode in ('steady-warm','steady-refine','steady-high-order')}
        output = OUT/f'{mode}-worker-{index}.json'
        command = [sys.executable, '-X', 'utf8', str(Path(__file__).resolve()), mode,
                   '--job', json.dumps(job), '--job-output', str(output)]
        # File-backed child dispatch avoids Windows named-pipe restrictions.
        with (OUT/f'{mode}-worker-{index}.log').open('w', encoding='utf-8') as log:
            subprocess.run(command, stdout=log, stderr=log, check=True)
        return json.loads(output.read_text(encoding='utf-8'))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(child, i, c) for i, c in enumerate(jobs)]
        for f in as_completed(futures):
            row = f.result()
            rows.append(row)
            payload['complete'] = len(rows) == len(jobs)
            (OUT/(mode+'-comparison.json')).write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')
            print(f"{len(rows)}/{len(jobs)} {row['file']} converged={row['converged']} residual={row.get('residual_rms_max_equation', row.get('failure'))} p={row.get('wall_pressure_ratio')} time={row.get('elapsed_seconds')}", flush=True)
    return payload


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['euler','viscous','refine','high-order','characteristic','primitive','steady','steady-warm','steady-refine','steady-high-order'])
    p.add_argument('--workers', type=int, default=3)
    p.add_argument('--grids', type=int, nargs='+', default=[60,120])
    p.add_argument('--job')
    p.add_argument('--job-output')
    a = p.parse_args()
    if a.job:
        job = json.loads(a.job)
        c = (ViscousConfig if job['viscous'] else Config)(**job['config'])
        Path(a.job_output).write_text(json.dumps(execute(c,job.get('steady',False),job.get('warm',False)), indent=2)+'\n', encoding='utf-8')
    else:
        run(a.mode, a.workers, a.grids)
