"""Transient FV/RKDG shock-tube suite using the cone's shared numerical fluxes.

No geometric source and no physical viscosity: this isolates Euler wave capture.
Each run records final time, conservation budget, positivity, error norms and
actual history. Unsteady RHS norms are not interpreted as steady convergence.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
from functools import lru_cache
from time import perf_counter
import json
import numpy as np
from state import conserved,primitive,admissible,physical_flux
from reconstruction import coefficients,evaluate,limit_admissibility
from flux import FLUXES,numerical_flux
from additional_fluxes import riemann_sample,jst_flux
from advanced_reconstruction import RECONSTRUCTIONS
from dg import DGConicalSolver
from characteristics import basis
from numpy.polynomial.legendre import legvander,Legendre

ROOT=Path(__file__).resolve().parent
CASES={
 'sod':{'left':[1.,0.,0.,1.],'right':[.125,0.,0.,.1],'time':.2},
 'lax':{'left':[.445,0.,.698,3.528],'right':[.5,0.,0.,.571],'time':.14},
 'double-rarefaction':{'left':[1.,0.,-2.,.4],'right':[1.,0.,2.,.4],'time':.1},
 'stationary-contact':{'left':[1.,0.,0.,1.],'right':[.125,0.,0.,1.],'time':.2},
}

@dataclass
class TubeConfig:
    case:str='sod'
    cells:int=100
    flux:str='hllc'
    reconstruction:str='cweno3'
    gamma:float=1.4
    cfl:float=.25
    degree:int=-1
    tvb:float=0.
    max_steps:int=100000


@lru_cache(None)
def reference(case,cells,gamma):
    x=(np.arange(cells)+.5)/cells
    time=CASES[case]['time']
    # Integrate conserved exact solution, then convert the exact cell mean.
    def integrate(nq):
        nodes,w=np.polynomial.legendre.leggauss(nq)
        points=x[:,None]+nodes/(2*cells)
        q=riemann_sample(CASES[case]['left'],CASES[case]['right'],(points-.5)/time,gamma)
        return np.einsum('nqv,q->nv',conserved(q,gamma),w/2)
    avg=integrate(64)
    q=primitive(avg,gamma)
    check=primitive(integrate(128),gamma)
    return q,float(np.mean(abs(q-check)))


class TubeSolver:
    def __init__(self,c):
        if c.case not in CASES or c.cells<8: raise ValueError('Invalid shock-tube case/grid')
        if c.flux=='jst' and (c.degree>=0 or c.reconstruction!='first'):
            raise ValueError('JST is a standalone FV stencil; choose first, degree=-1')
        self.config=c;self.h=1/c.cells
        self.x=(np.arange(c.cells)+.5)*self.h
        initial=np.where((self.x<.5)[:,None],CASES[c.case]['left'],CASES[c.case]['right'])
        self.u=conserved(initial,c.gamma)
        self.initial_integral=self.u.sum(axis=0)*self.h
        self.budget=np.zeros(4)
        self.events={'rejected_steps':0,'reconstruction_limited_cells':0}
        self.time=0.
    def ghost(self,u): return np.pad(u,((3,3),(0,0)),mode='edge')
    def rhs(self,u):
        c=self.config;ext=self.ghost(u)
        if c.flux=='jst': f=jst_flux(ext,c.gamma)
        else:
            pol=coefficients(ext,c.reconstruction,self.h,c.gamma)
            pol,count=limit_admissibility(pol,ext[2:-2],c.gamma,[-.5,0,.5])
            self.events['reconstruction_limited_cells']+=count
            f=numerical_flux(c.flux,evaluate(pol,.5)[:-1],evaluate(pol,-.5)[1:],c.gamma,mach_ref=.1)
        return -np.diff(f,axis=0)/self.h,f[0]-f[-1]
    def limit(self,u):
        if not admissible(u,self.config.gamma): raise ArithmeticError('Inadmissible cell mean')
        return u
    def time_step(self,u):
        q=primitive(u,self.config.gamma)
        return self.config.cfl*self.h/np.max(abs(q[:,2])+np.sqrt(self.config.gamma*q[:,3]/q[:,0]))
    def step(self,u,dt):
        r0,b0=self.rhs(u);a=self.limit(u+dt*r0)
        r1,b1=self.rhs(a);b=self.limit(.75*u+.25*(a+dt*r1))
        r2,b2=self.rhs(b);v=self.limit(u/3+2*(b+dt*r2)/3)
        return v,dt*(b0+b1+4*b2)/6
    def averages(self): return self.u
    def run(self):
        start=perf_counter();end=CASES[self.config.case]['time'];history=[];failed=None;steps=0
        for iteration in range(self.config.max_steps):
            if self.time>=end-1e-14: break
            dt=min(self.time_step(self.u),end-self.time)
            for attempt in range(14):
                try:
                    u,increment=self.step(self.u,dt)
                    break
                except (ArithmeticError,FloatingPointError):
                    dt*=.5;self.events['rejected_steps']+=1
            else:
                failed='Positivity recovery exhausted';break
            if dt<1e-14: failed='Time step collapsed';break
            self.u=u;self.budget+=increment;self.time=float(self.time+dt);steps+=1
            if iteration%20==0 or self.time>=end-1e-14:
                q=primitive(self.averages(),self.config.gamma)
                defect=np.max(abs(self.averages().sum(axis=0)*self.h-self.initial_integral-self.budget))
                history.append([self.time,float(q[:,0].min()),float(q[:,3].min()),float(defect)])
        else: failed='Step budget exhausted'
        avg=self.averages();q=primitive(avg,self.config.gamma)
        target,quad=reference(self.config.case,self.config.cells,self.config.gamma)
        reached=bool(self.time>=end-1e-12)
        errors={}
        if reached:
            for name,k in [('density',0),('velocity',2),('pressure',3)]:
                d=abs(q[:,k]-target[:,k]);errors[name]={'L1':float(d.mean()),'L2':float(np.sqrt(np.mean(d*d))),'Linf':float(d.max())}
        defect=float(np.max(abs(avg.sum(axis=0)*self.h-self.initial_integral-self.budget)))
        valid=bool(reached and admissible(avg,self.config.gamma) and defect<1e-10)
        return {'config':asdict(self.config),'completed':reached,'integrity_passed':valid,
                'termination':failed or 'Final physical time reached','time':self.time,
                'steps':steps,'elapsed_seconds':perf_counter()-start,'events':self.events,
                'conservation_budget_defect':defect,'min_density':float(q[:,0].min()),
                'min_pressure':float(q[:,3].min()),'errors':errors,'reference_quadrature_change':quad,
                'scope':'Inviscid transient planar Euler; completion/integrity is not an accuracy certification',
                'x':self.x,'primitive':q,'conserved_average':avg,'exact':target,'history':np.array(history)}


class TubeDG(TubeSolver,DGConicalSolver):
    def __init__(self,c):
        TubeSolver.__init__(self,c)
        self.nodes,self.weights=np.polynomial.legendre.leggauss(max(3,c.degree+2))
        self.v=legvander(self.nodes,c.degree)
        self.dv=np.column_stack([Legendre.basis(k).deriv()(self.nodes) for k in range(c.degree+1)])
        self.ends=legvander(np.array([-1.,1.]),c.degree)
        self.sample=legvander(np.r_[-1.,self.nodes,1.],c.degree)
        self.scale=2*np.arange(c.degree+1)+1
        avg=self.u.copy();self.u=np.zeros((c.cells,c.degree+1,4));self.u[:,0]=avg
        self.events.update(slope_limited_cells=0,positivity_scaled_cells=0)
    project=DGConicalSolver.project
    values=DGConicalSolver.values
    def averages(self): return self.u[:,0]
    def limit(self,modes):
        from dg import minmod
        c=self.config;result=modes.copy();mean=result[:,0]
        if not admissible(mean,c.gamma): raise ArithmeticError('DG mean inadmissible')
        if c.degree:
            prev=np.r_[mean[:1],mean[:-1]];nxt=np.r_[mean[1:],mean[-1:]]
            r,l=basis(mean,c.gamma)
            a=np.einsum('nij,nj->ni',l,result[:,1])
            b=np.einsum('nij,nj->ni',l,mean-prev);d=np.einsum('nij,nj->ni',l,nxt-mean)
            limited=np.where(abs(a)<=c.tvb*self.h**2,a,minmod(a,b,d))
            changed=np.max(abs(limited-a),axis=1)>1e-12
            result[changed,1]=np.einsum('nij,nj->ni',r[changed],limited[changed])
            result[changed,2:]=0
            self.events['slope_limited_cells']+=int(sum(changed))
        samples=self.values(result,self.sample)
        q=primitive(samples,c.gamma)
        bad=np.any((q[...,0]<=1e-12)|(q[...,3]<=1e-12)|(~np.isfinite(q).all(axis=-1)),axis=1)
        if np.any(bad):
            lo,hi=np.zeros(sum(bad)),np.ones(sum(bad));delta=samples[bad]-mean[bad,None]
            for _ in range(36):
                mid=(lo+hi)/2;trial=mean[bad,None]+mid[:,None,None]*delta
                tq=primitive(trial,c.gamma);good=np.all((tq[...,0]>1e-12)&(tq[...,3]>1e-12),axis=1)
                lo,hi=np.where(good,mid,lo),np.where(good,hi,mid)
            result[bad,1:]*=(.999*lo)[:,None,None]
            self.events['positivity_scaled_cells']+=int(sum(bad))
        return result
    def rhs(self,modes):
        c=self.config;values=self.values(modes);ends=self.values(modes,self.ends)
        left=np.r_[ends[:1,0],ends[:,1]];right=np.r_[ends[:,0],ends[-1:,1]]
        f=numerical_flux(c.flux,left,right,c.gamma,mach_ref=.1)
        volume=np.einsum('nqv,qk,q->nkv',physical_flux(values,c.gamma),self.dv,self.weights)
        boundary=f[1:,None]*self.ends[1,:,None]-f[:-1,None]*self.ends[0,:,None]
        return (volume-boundary)*self.scale[None,:,None]/self.h,f[0]-f[-1]
    def time_step(self,u):
        q=primitive(self.values(u,self.sample),self.config.gamma)
        speed=np.max(abs(q[...,2])+np.sqrt(self.config.gamma*q[...,3]/q[...,0]))
        return self.config.cfl*self.h/((2*self.config.degree+1)*speed)


def run_case(flux='hllc',reconstruction='cweno3',cells=100,case='sod',degree=-1,output_dir=None):
    c=TubeConfig(flux=flux,reconstruction=reconstruction,cells=cells,case=case,degree=degree)
    solver=TubeDG(c) if degree>=0 else TubeSolver(c)
    r=solver.run();target=Path(output_dir or ROOT/'results-shocktube');target.mkdir(parents=True,exist_ok=True)
    name=f'{case}-{flux}-'+(f'dg{degree}' if degree>=0 else reconstruction)+f'-n{cells}'
    keys=['x','primitive','conserved_average','exact','history']
    np.savez_compressed(target/(name+'.npz'),**{k:r[k] for k in keys})
    meta={k:v for k,v in r.items() if k not in keys}
    (target/(name+'.json')).write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    return r

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--flux',choices=[*FLUXES,'jst'],default='hllc')
    p.add_argument('--reconstruction',choices=RECONSTRUCTIONS,default='cweno3')
    p.add_argument('--cells',type=int,default=100);p.add_argument('--case',choices=CASES,default='sod')
    p.add_argument('--degree',type=int,default=-1,choices=[-1,0,1,2]);a=p.parse_args()
    r=run_case(**vars(a));print(json.dumps({k:v for k,v in r.items() if not isinstance(v,np.ndarray)},indent=2))
