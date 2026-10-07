"""Bounded, reproducible DG cone study; retain every nonconverged result."""
from pathlib import Path
import json
import hashlib
from dg import DGConfig, DGConicalSolver, save_dg
from reference import taylor_maccoll
from study import metrics

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results-dg'


def main():
    OUT.mkdir(exist_ok=True)
    ref = taylor_maccoll(2.35,max_step=1e-4)
    rows = []
    for flux in ['hllc','ausm-up2']:
        for degree in [0,1,2]:
            name = f'cone-m2.35-n40-{flux}-p{degree}'
            cfg = DGConfig(flux=flux,degree=degree,cells=40,max_steps=8000)
            solver = DGConicalSolver(cfg)
            result = metrics(solver.run(True),ref)
            save_dg(result,OUT/name)
            rows.append({k:v for k,v in result.items() if k not in ('theta','primitive','conserved_average','history','modal_coefficients')} | {'file':name})
            (OUT/'comparison.json').write_text(json.dumps({'runs':rows,'complete':len(rows)==6,
                'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'dg.py',ROOT/'verify_dg.py',ROOT/'dg_study.py']},
                'reference_wall_pressure_ratio':ref['wall_pressure_ratio'],
                'date':'2026-10-07','initialization':'Freestream Legendre L2 projection, identical for each degree/flux'},indent=2)+'\n',encoding='utf-8')
            print(name,result['converged'],result['residual_rms_max_equation'],flush=True)


if __name__ == '__main__':
    main()
