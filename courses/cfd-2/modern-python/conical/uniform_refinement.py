"""Fixed uniform-grid refinement using the original planar TubeSolver.

The mesh never moves, splits or merges. All runs use the same CFL, final time
and discretization while doubling N. Metrics use computed FV cell means.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import numpy as np
from shock_suite import ROOT, CASES, FLUXES, TubeConfig, TubeSolver, run_case
from additional_fluxes import wave

GRIDS = [80, 160, 320, 640, 1280, 2560]
METHODS = [*FLUXES, 'jst']


def features(gamma=1.4):
    dl, _, ul, pl = CASES['sod']['left']; dr, _, ur, pr = CASES['sod']['right']
    al, ar = np.sqrt(gamma*pl/dl), np.sqrt(gamma*pr/dr)
    lo, hi = 0., max(pl, pr)
    for _ in range(70):
        p = (lo+hi)/2
        f = wave(p,dl,pl,al,gamma)+wave(p,dr,pr,ar,gamma)+ur-ul
        if f < 0: lo=p
        else: hi=p
    ps = (lo+hi)/2
    us = (ul+ur+wave(ps,dr,pr,ar,gamma)-wave(ps,dl,pl,al,gamma))/2
    ratio=ps/pr; z=(gamma-1)/(gamma+1)
    ds=dr*(ratio+z)/(z*ratio+1)
    speed=ur+ar*np.sqrt((gamma+1)/(2*gamma)*ratio+(gamma-1)/(2*gamma))
    t=CASES['sod']['time']
    return {'contact':{'x':.5+us*t,'high':dl*(ps/pl)**(1/gamma),'low':ds},
            'shock':{'x':.5+speed*t,'high':ds,'low':dr}}


def jump_widths(x, rho):
    """Apparent density 10–90% widths, interpolated between cell centers."""
    answer={}
    for name,f in features().items():
        crossings=[]
        for fraction in [.9,.1]:
            level=f['low']+fraction*(f['high']-f['low'])
            mask=(rho[:-1]>=level)&(rho[1:]<=level)&(abs(.5*(x[:-1]+x[1:])-f['x'])<.065)
            indices=np.flatnonzero(mask)
            if len(indices):
                points=x[indices]+(rho[indices]-level)/(rho[indices]-rho[indices+1])*(x[indices+1]-x[indices])
                crossings.append(float(points[np.argmin(abs(points-f['x']))]))
        width=max(0.,crossings[1]-crossings[0]) if len(crossings)==2 else None
        answer[name]={'exact_location':float(f['x']),'width_10_90':width,
                      'width_in_cells':width*len(x) if width is not None else None}
    return answer


def compute(args):
    flux,n=args; rec='first' if flux=='jst' else 'cweno3'
    name=f'sod-{flux}-{rec}-n{n}'
    target=ROOT/'results-uniform';target.mkdir(exist_ok=True)
    metadata=target/(name+'.json');field=target/(name+'.npz')
    if metadata.exists() and field.exists():
        row=json.loads(metadata.read_text())
        with np.load(field) as data: row.update({k:data[k].copy() for k in data.files})
    else:
        original=ROOT/'results-shocktube'
        if (original/(name+'.json')).exists() and (original/(name+'.npz')).exists():
            row=json.loads((original/(name+'.json')).read_text())
            with np.load(original/(name+'.npz')) as data: row.update({k:data[k].copy() for k in data.files})
            row['provenance']='Reused matching executed original uniform TubeSolver run'
        else:
            row=run_case(flux,rec,n,output_dir=target)
            row['provenance']='Executed fixed uniform TubeSolver run for this refinement study'
        row['jump_widths']=jump_widths(row['x'],row['primitive'][:,0])
        row['mesh']={'kind':'fixed uniform','cells':n,'dx':1/n,'adaptation':False}
        np.savez_compressed(field,**{k:v for k,v in row.items() if isinstance(v,np.ndarray)})
        metadata.write_text(json.dumps({k:v for k,v in row.items() if not isinstance(v,np.ndarray)},indent=2)+'\n')
    q=row['primitive']; x=row['x']
    assert row['config']==vars(TubeConfig(flux=flux,reconstruction=rec,cells=n))
    assert row['integrity_passed'] and abs(row['time']-.2)<1e-12
    assert np.isfinite(q).all() and (q[:,0]>0).all() and (q[:,3]>0).all()
    assert np.allclose(x,(np.arange(n)+.5)/n,atol=1e-14,rtol=0)
    assert row['conservation_budget_defect']<1e-10
    entry={'name':name,'config':row['config'],'integrity_passed':True}
    for suffix in ['json','npz']:
        entry[suffix+'_sha256']=hashlib.sha256((target/(name+'.'+suffix)).read_bytes()).hexdigest()
    return entry


def main():
    rows=[]; jobs=[(f,n) for n in GRIDS for f in METHODS]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for i,future in enumerate(as_completed([pool.submit(compute,a) for a in jobs]),1):
            row=future.result();rows.append(row)
            print(f'{i}/{len(jobs)} {row["name"]}: passed',flush=True)
    catalog={'complete':True,'grids':GRIDS,'methods':METHODS,'runs':sorted(rows,key=lambda r:r['name']),
             'scope':'96 fixed uniform Sod runs: 15 fluxes with CWENO3 and standalone JST; no adaptive mesh'}
    (ROOT/'results-uniform/catalog.json').write_text(json.dumps(catalog,indent=2)+'\n')
    print(json.dumps({'runs':len(rows),'passed':sum(r['integrity_passed'] for r in rows)}),flush=True)


if __name__=='__main__': main()
