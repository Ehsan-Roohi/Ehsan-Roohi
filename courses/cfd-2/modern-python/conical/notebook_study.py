"""Fill the three missing flux comparisons and build the complete notebook catalog."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from colab_api import run_case
from flux import FLUXES
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results-notebook'


def jobs():
    result = []
    for flux in ('vanleer','hlle','rusanov'):
        for mach,outer in [(2.35,45),(7.95,22)]:
            result.append(dict(flux=flux,mach=mach,outer_deg=outer,cells=60,max_steps=18000,
                               output_dir=str(OUT/f'm{mach:g}')))
        for wall in ('adiabatic','isothermal'):
            result.append(dict(flux=flux,physics=wall,mach=7.95,outer_deg=22,cells=60,
                               reconstruction='primitive-minmod',max_newton=180,
                               output_dir=str(OUT/wall)))
    for reconstruction in ('cweno3-char','cweno5-char'):
        result.append(dict(flux='hllc',mach=2.35,outer_deg=45,cells=60,max_steps=18000,
                           reconstruction=reconstruction,output_dir=str(OUT/'reconstruction')))
    return result


def record(stem, group):
    data = json.loads(Path(str(stem)+'.json').read_text(encoding='utf-8'))
    # Decimal Mach basenames require appending a suffix, not replacing it.
    npz = Path(str(stem)+'.npz')
    if not npz.exists():
        raise FileNotFoundError(npz)
    q = np.load(npz)['primitive']
    assert np.isfinite(q).all() and (q[:,0]>0).all() and (q[:,3]>0).all()
    return dict(group=group,json_path=Path(str(stem)+'.json').relative_to(ROOT).as_posix(),
                npz_path=npz.relative_to(ROOT).as_posix(),config=data['config'],
                converged=data['converged'],sha256=hashlib.sha256(npz.read_bytes()).hexdigest())


def catalog():
    rows = []
    primary_effective = {}
    for pattern in ['fv2-steady-*.json','continued-fv2-steady-*.json','mesh-continued-fv2-steady-*.json']:
        for p in sorted((ROOT/'results-v2').glob(pattern)):
            c = json.loads(p.read_text())['config']
            if c['cells']==60 and c['reconstruction']=='primitive-minmod' and c['gradient_order']==2:
                primary_effective[(c['thermal_wall'],c['flux'])] = Path(str(p)[:-5])
    for mach in (2.35,7.95):
        for flux in FLUXES:
            stem = (ROOT/'results-v2'/f'euler-m{mach:g}-n60-{flux}-cweno3'
                    if flux not in ('vanleer','hlle','rusanov') else OUT/f'm{mach:g}'/f'inviscid-{flux}-n60-cweno3')
            rows.append(record(stem,'flux'))
    for wall in ('adiabatic','isothermal'):
        for flux in FLUXES:
            stem = (primary_effective[(wall,flux)]
                    if flux not in ('vanleer','hlle','rusanov') else OUT/wall/f'{wall}-{flux}-n60-primitive-minmod')
            rows.append(record(stem,'flux'))
    for reconstruction in ('first','muscl','cweno3','cweno5'):
        rows.append(record(ROOT/'results'/f'm2.35-n60-hllc-{reconstruction}','reconstruction'))
    for reconstruction in ('cweno3-char','cweno5-char'):
        rows.append(record(OUT/'reconstruction'/f'inviscid-hllc-n60-{reconstruction}','reconstruction'))
    for stem in sorted((ROOT/'results-dg').glob('cone-*.json')):
        rows.append(record(Path(str(stem)[:-5]),'dg'))
    # Keep the actual accepted fine-grid continuation where one exists.
    effective = {}
    for pattern in ['fv2-steady-*.json','continued-fv2-steady-*.json','mesh-continued-fv2-steady-*.json']:
        for p in sorted((ROOT/'results-v2').glob(pattern)):
            c = json.loads(p.read_text())['config']
            if c['gradient_order']==4:
                effective[(c['thermal_wall'],c['cells'],c['flux'])] = p
    for p in effective.values():
        rows.append(record(Path(str(p)[:-5]),'viscous-high-order'))
    payload = {'complete':True,'date':'2026-10-07','runs':rows,
               'fluxes':list(FLUXES),'flux_comparison_runs':44,'reconstruction_runs':6,
               'dg_runs':6,'viscous_high_order_runs':len(effective),
               'note':'Executed results loaded by default; optional reruns are separate. Failed gates remain visible.',
               'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}}
    (OUT/'catalog.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in payload.items() if k not in ('runs','source_sha256')},indent=2))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--job',type=int)
    p.add_argument('--catalog-only',action='store_true')
    a = p.parse_args()
    OUT.mkdir(exist_ok=True)
    if a.job is not None:
        run_case(**jobs()[a.job],verbose=True)
    elif a.catalog_only:
        catalog()
    else:
        def child(index):
            with (OUT/f'job-{index}.log').open('w',encoding='utf-8') as log:
                subprocess.run([sys.executable,'-X','utf8',str(Path(__file__).resolve()),'--job',str(index)],
                               stdout=log,stderr=log,check=True)
            print(f'Completed notebook comparison {index+1}/{len(jobs())}',flush=True)
        with ThreadPoolExecutor(max_workers=2) as pool:
            for future in as_completed([pool.submit(child,i) for i in range(len(jobs()))]):
                future.result()
        catalog()


if __name__ == '__main__':
    main()
