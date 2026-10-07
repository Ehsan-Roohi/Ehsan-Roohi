"""Cone extension study: newly added fluxes and FV reconstructions, matched N60."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
from dataclasses import asdict
import numpy as np,json,hashlib
from solver import Config,ConicalSolver,save_result
from viscous import ViscousConfig,ViscousConicalSolver
from advanced_reconstruction import RECONSTRUCTIONS
from study import metrics
from reference import taylor_maccoll
ROOT=Path(__file__).resolve().parent
NEW_FLUXES=['godunov','roe-nc','steger-warming','global-lf']
NEW_RECONS=RECONSTRUCTIONS[6:]


def job(args):
    flux,rec,mach,physics=args
    target=ROOT/'results-additional';target.mkdir(exist_ok=True)
    name=f'{physics}-m{mach:g}-{flux}-{rec}-n60'
    path=target/(name+'.json')
    if path.exists():return name,json.loads(path.read_text())
    if physics=='inviscid':
        c=Config(flux=flux,reconstruction=rec,mach=mach,outer_deg=45 if mach<3 else 22,cells=60,max_steps=6000,cfl=.25)
        solver=ConicalSolver(c)
        catalog=json.loads((ROOT/'results-notebook/catalog.json').read_text())
        seed=next(r for r in catalog['runs'] if r['group']=='flux' and r['config']['flux']=='hllc' and r['config']['mach']==mach and 'thermal_wall' not in r['config'])
        seed_path=ROOT/seed['npz_path'];solver.u=np.load(seed_path)['conserved_average'].copy()
        result=solver.run(False)
        result=metrics(result,taylor_maccoll(mach,max_step=1e-4))
    else:
        c=ViscousConfig(flux=flux,thermal_wall=physics,cells=60)
        solver=ViscousConicalSolver(c)
        seed_path=ROOT/'results-v2'/f'fv2-steady-{physics}-m7.95-n60-ausm-up2-primitive-minmod-g2.npz'
        solver.u=np.load(seed_path)['conserved_average'].copy()
        from steady import solve
        result=solve(solver,max_iterations=180,verbose=False)
    result['initialization']='Configuration-matched published coarse solution; initializer only'
    result['initial_checkpoint_sha256']=hashlib.sha256(seed_path.read_bytes()).hexdigest()
    result['comparison_group']='new-reconstruction' if rec!='cweno3' and flux=='hllc' else 'new-flux'
    if flux=='jst':result['comparison_group']='standalone-jst'
    save_result(result,target/name)
    return name,{k:v for k,v in result.items() if k not in ['theta','primitive','conserved_average','history']}


def main():
    jobs=[]
    for mach in [2.35,7.95]:
        jobs.extend((f,'cweno3',mach,'inviscid') for f in NEW_FLUXES)
        jobs.extend(('hllc',r,mach,'inviscid') for r in NEW_RECONS)
        jobs.append(('jst','first',mach,'inviscid'))
    for wall in ['adiabatic','isothermal']:jobs.extend((f,'primitive-minmod',7.95,wall) for f in NEW_FLUXES)
    print(f'Starting {len(jobs)} new cone cases',flush=True)
    rows=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        future={pool.submit(job,j):j for j in jobs}
        for f in as_completed(future):
            args=future[f]
            try:
                name,r=f.result();rows.append({'name':name,'config':r['config'],'converged':r['converged'],'group':r['comparison_group']})
                print(f'{name}: gate={r["converged"]} R={r["residual_rms_max_equation"]:.3e}',flush=True)
            except Exception as e:
                name=f'{args[3]}-m{args[2]:g}-{args[0]}-{args[1]}-n60'
                r={'name':name,'config':{'flux':args[0],'reconstruction':args[1],'mach':args[2],'physics':args[3],'cells':60},'converged':False,'error':str(e)}
                (ROOT/'results-additional'/(name+'.json')).write_text(json.dumps(r,indent=2)+'\n');rows.append(r)
                print(name+': '+str(e),flush=True)
    (ROOT/'results-additional/catalog.json').write_text(json.dumps({'complete':True,'runs':rows,'new_fluxes':NEW_FLUXES,'new_reconstructions':NEW_RECONS},indent=2)+'\n')
    print(json.dumps({'runs':len(rows),'accepted':sum(r['converged'] for r in rows)},indent=2),flush=True)
if __name__=='__main__':main()
