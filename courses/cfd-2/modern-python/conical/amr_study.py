"""Same FV geometry-aware operator: uniform 80/160/320/640 versus local AMR."""
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,hashlib
from shock_amr import run_case,FLUXES,ROOT


def job(args):
    flux,n,level=args;name=f'sod-{flux}-cweno3-base{n}-level{level}'
    path=ROOT/'results-amr'/(name+'.json')
    if path.exists(): meta=json.loads(path.read_text())
    else:
        r=run_case(flux=flux,base_cells=n,max_level=level)
        meta={k:v for k,v in r.items() if not hasattr(v,'shape')}
    return name,meta


def main():
    jobs=[(f,n,0) for f in FLUXES for n in [80,160,320,640]]+[(f,80,3) for f in FLUXES]
    rows=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(job,a) for a in jobs]
        for i,future in enumerate(as_completed(futures),1):
            name,r=future.result()
            rows.append({'name':name,'config':r['config'],'completed':r['completed'],'integrity_passed':r['integrity_passed']})
            print(f'{i}/{len(jobs)} {name}: {r["termination"]}',flush=True)
    for row in rows:
        for ext in ['json','npz']:row[ext+'_sha256']=hashlib.sha256((ROOT/'results-amr'/(row['name']+'.'+ext)).read_bytes()).hexdigest()
    rows.sort(key=lambda r:r['name'])
    (ROOT/'results-amr/catalog.json').write_text(json.dumps({'complete':True,'runs':rows,
        'scope':'75 Sod runs: all 15 pointwise fluxes with nonuniform-compatible CWENO3, four uniform controls and density/pressure flagged local AMR'},indent=2)+'\n')
    print(json.dumps({'runs':len(rows),'integrity_passed':sum(r['integrity_passed'] for r in rows)},indent=2))


if __name__=='__main__':main()
