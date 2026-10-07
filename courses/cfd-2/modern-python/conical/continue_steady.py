"""Configuration-matched continuation, preserving the original failed records."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from viscous import ViscousConfig, ViscousConicalSolver
from steady import solve, remap_checkpoint
from solver import save_result

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results-v2'


def main(pattern, iterations, reseed=False):
    rows = []
    for path in sorted(OUT.glob(pattern)):
        original = json.loads(path.read_text(encoding='utf-8'))
        if original['converged']:
            continue
        target = OUT/(('mesh-continued-' if reseed else 'continued-')+path.stem)
        if Path(str(target)+'.json').exists():
            rows.append(json.loads(Path(str(target)+'.json').read_text(encoding='utf-8')))
            continue
        solver = ViscousConicalSolver(ViscousConfig(**original['config']))
        data = np.load(path.with_suffix('.npz'))
        assert np.max(abs(data['theta']-solver.theta)) < 1e-14
        solver.u = data['conserved_average'].copy()
        if reseed:
            wall = original['config']['thermal_wall']
            seed = OUT/f'fv2-steady-{wall}-m7.95-n120-ausm-up2-cweno5-char-g4.npz'
            remap_checkpoint(solver,seed)
        result = solve(solver,iterations)
        input_path = seed if reseed else path.with_suffix('.npz')
        result.update(original_attempt=path.stem,
                      restarted_from_grid_seed=seed.stem if reseed else None,
                      continued_from=None if reseed else path.stem,
                      initial_checkpoint_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
                      source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')})
        if not reseed:
            history = np.array(result['history'])
            history[:,0] += original['iterations']
            result['history'] = np.r_[data['history'],history[1:]].tolist()
            result['iterations'] += original['iterations']
            result['elapsed_seconds'] += original['elapsed_seconds']
            result['explicit_recovery_steps'] += original.get('explicit_recovery_steps',0)
        save_result(result,target)
        rows.append({k:v for k,v in result.items() if k not in ('theta','primitive','conserved_average','history')} | {'file':target.name})
        (OUT/'continuations.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
        print(f"{path.stem}: converged={result['converged']} residual={result['residual_rms_max_equation']:.3e}",flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--pattern',default='fv2-steady-*-primitive-minmod-g2.json')
    p.add_argument('--iterations',type=int,default=240)
    p.add_argument('--reseed',action='store_true')
    a = p.parse_args()
    main(a.pattern,a.iterations,a.reseed)
