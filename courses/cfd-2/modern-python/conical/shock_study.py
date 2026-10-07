"""Execute every shared flux/reconstruction on Sod, plus grids/DG/stress cases."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import json,hashlib,os
from shock_suite import run_case,FLUXES,RECONSTRUCTIONS,ROOT


def job(args):
    name,rec,n,case,degree=args
    tag=f'{case}-{name}-'+(f'dg{degree}' if degree>=0 else rec)+f'-n{n}'
    target=ROOT/'results-shocktube'
    path=target/(tag+'.json')
    if path.exists(): return tag,json.loads(path.read_text())
    try:
        r=run_case(name,rec,n,case,degree)
        return tag,{k:v for k,v in r.items() if k not in ['x','primitive','conserved_average','exact','history']}
    except Exception as e:
        error={'config':{'flux':name,'reconstruction':rec,'cells':n,'case':case,'degree':degree},
               'completed':False,'integrity_passed':False,'termination':type(e).__name__+': '+str(e)}
        target.mkdir(exist_ok=True)
        path.write_text(json.dumps(error,indent=2)+'\n')
        return tag,error


def main():
    jobs=set()
    for f in FLUXES:
        for rec in RECONSTRUCTIONS: jobs.add((f,rec,80,'sod',-1))
        for n in [40,80,160]: jobs.add((f,'cweno3',n,'sod',-1))
        for case in ['lax','double-rarefaction','stationary-contact']: jobs.add((f,'cweno3',80,case,-1))
        for p in [0,1,2]:
            for n in [40,80,160]: jobs.add((f,'modal-dg',n,'sod',p))
    for rec in RECONSTRUCTIONS:
        for n in [40,80,160]: jobs.add(('hllc',rec,n,'sod',-1))
        for case in ['lax','double-rarefaction']: jobs.add(('hllc',rec,80,case,-1))
    for n in [40,80,160]: jobs.add(('jst','first',n,'sod',-1))
    for case in ['lax','double-rarefaction','stationary-contact']: jobs.add(('jst','first',80,case,-1))
    print(f'Starting {len(jobs)} shock-tube runs, four CPU workers',flush=True)
    rows=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(job,args):args for args in sorted(jobs)}
        for i,f in enumerate(as_completed(futures),1):
            tag,r=f.result();rows.append({'name':tag,'config':r['config'],'completed':r['completed'],
                                         'integrity_passed':r['integrity_passed']})
            if i%15==0 or not r['completed']: print(f'{i}/{len(jobs)} {tag}: {r["termination"]}',flush=True)
    rows.sort(key=lambda r:r['name'])
    for r in rows:
        for suffix in ['json','npz']:
            p=ROOT/'results-shocktube'/(r['name']+'.'+suffix)
            if p.exists(): r[suffix+'_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
    result={'complete':True,'runs':rows,'fluxes':list(FLUXES)+['jst'],'reconstructions':RECONSTRUCTIONS,
            'scope':'Full 15-flux x 20-reconstruction Sod matrix at N80; standalone JST; DG all shared fluxes; targeted grid and stress studies'}
    (ROOT/'results-shocktube/catalog.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'runs':len(rows),'completed':sum(r['completed'] for r in rows),'integrity_passed':sum(r['integrity_passed'] for r in rows)},indent=2),flush=True)

if __name__=='__main__': main()
