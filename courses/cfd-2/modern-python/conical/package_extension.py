"""Audit actual fields/statuses and package source, reports and raw results."""
from pathlib import Path
import hashlib
import json
import zipfile
import numpy as np
from solver import save_result

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results-v2'


def repair_restart_histories():
    # Early mesh restart bookkeeping joined two different initial trajectories.
    # Correct only counters/provenance; preserve measured fields and gate status.
    from viscous import ViscousConfig, ViscousConicalSolver
    from steady import remap_checkpoint, norm
    for path in OUT.glob('mesh-continued-fv2-steady-*.json'):
        meta = json.loads(path.read_text(encoding='utf-8'))
        if meta.get('restarted_from_grid_seed'):
            continue
        original_path = OUT/(path.stem.removeprefix('mesh-continued-')+'.json')
        original = json.loads(original_path.read_text(encoding='utf-8'))
        old = np.load(original_path.with_suffix('.npz'))
        data = np.load(path.with_suffix('.npz'))
        c = ViscousConfig(**meta['config'])
        solver = ViscousConicalSolver(c)
        seed = OUT/f'fv2-steady-{c.thermal_wall}-m7.95-n120-ausm-up2-cweno5-char-g4.npz'
        remap_checkpoint(solver,seed)
        r = solver.rhs(solver.u,record=False)
        wall = solver.wall_state()
        tail = data['history'][len(old['history']):].copy()
        tail[:,0] -= original['iterations']
        hist = np.r_[[[0,norm(r),wall[3]*c.gamma*c.mach**2]],tail]
        meta['iterations'] -= original['iterations']
        meta['elapsed_seconds'] -= original['elapsed_seconds']
        meta['explicit_recovery_steps'] -= original.get('explicit_recovery_steps',0)
        meta.update(continued_from=None,original_attempt=original_path.stem,restarted_from_grid_seed=seed.stem,
                    initial_checkpoint_sha256=hashlib.sha256(seed.read_bytes()).hexdigest(),
                    metadata_correction='Separated mesh-restart history from earlier failed trajectory; measured fields and residual unchanged')
        save_result(meta | {k:data[k] for k in ('theta','primitive','conserved_average')} | {'history':hist.tolist()},path.with_suffix(''))


def main():
    repair_restart_histories()
    records = []
    for path in sorted(OUT.glob('*.npz')):
        data = np.load(path)
        if 'primitive' not in data or not path.with_suffix('.json').exists():
            continue
        meta = json.loads(path.with_suffix('.json').read_text(encoding='utf-8'))
        q = data['primitive']
        assert np.all(np.isfinite(q)) and np.all(q[:,0] > 0) and np.all(q[:,3] > 0),path
        if meta.get('converged'):
            assert meta['residual_rms_max_equation'] < meta['config']['tolerance'],path
            if meta.get('steady_solver'):
                assert meta['explicit_20_step_max_state_drift'] < 1e-7,path
        records.append({'file':path.name,'positive_finite':True,'converged':meta.get('converged',False)})
    (OUT/'field-integrity.json').write_text(json.dumps({'checked_runs':len(records),'all_passed':True,'runs':records},indent=2)+'\n',encoding='utf-8')
    from dg import DGConfig, DGConicalSolver
    from state import admissible
    dg_records = []
    for path in sorted((ROOT/'results-dg').glob('cone-*.json')):
        meta = json.loads(path.read_text(encoding='utf-8'))
        data = np.load(path.with_suffix('.npz'))
        modal = np.load(path.with_name(path.stem+'-modes.npz'))['modes']
        solver = DGConicalSolver(DGConfig(**meta['config']))
        assert admissible(data['conserved_average']),path
        assert admissible(solver.values(modal,solver.sample)),path
        if meta['converged']:
            assert meta['residual_rms_max_equation'] < meta['config']['tolerance'],path
            assert meta['explicit_20_step_max_state_drift'] < 1e-7,path
        dg_records.append({'file':path.name,'converged':meta['converged'],'positive_finite':True})
    (ROOT/'results-dg/field-integrity.json').write_text(json.dumps({'all_passed':True,'runs':dg_records},indent=2)+'\n',encoding='utf-8')
    notebook_check = json.loads((ROOT/'colab-validation.json').read_text(encoding='utf-8'))
    assert notebook_check['passed'] and notebook_check['default_run_converged']
    archive = ROOT.parents[2]/'conical-ausm-viscous-results.zip'
    entries = {}
    files = [p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts
             and p.suffix not in ('.log','.flag') and '-worker-' not in p.name]
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(files):
            name = 'CFD-2/modern-python/conical/'+p.relative_to(ROOT).as_posix()
            z.write(p,name)
            entries[name] = {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
        legacy = ROOT.parents[1]/'projects/01-conical-flow/all-codes/Second0.for'
        name = 'CFD-2/projects/01-conical-flow/all-codes/Second0.for'
        z.write(legacy,name)
        entries[name] = {'bytes':legacy.stat().st_size,'sha256':hashlib.sha256(legacy.read_bytes()).hexdigest()}
        manifest = {'date':'2026-10-07','files':entries,'note':'Includes current Python code and retained earlier inviscid results. The original inviscid source release remains in the separate conical-python-and-results.zip. Failed exploratory runs are preserved but excluded from corrected-wall accuracy acceptance.'}
        z.writestr('MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
    with zipfile.ZipFile(archive) as z:
        for name,entry in entries.items():
            assert hashlib.sha256(z.read(name)).hexdigest() == entry['sha256']
    print(json.dumps({'archive':str(archive),'files':len(entries),'bytes':archive.stat().st_size,'checked_runs':len(records)},indent=2))


if __name__ == '__main__':
    main()
