"""Conservative leaf-cell h refinement for the shared planar Euler fluxes.

One active nonoverlapping partition and a global SSPRK3 step: each interface
uses one shared flux, so no coarse/fine subcycling or reflux correction is needed.
Density and pressure jumps flag cells; the exact answer never flags the mesh.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
from time import perf_counter
import json
import numpy as np
from fv_mesh import Mesh
from reconstruction import evaluate, limit_admissibility
from state import primitive, conserved, admissible
from flux import numerical_flux, FLUXES
from shock_suite import CASES, ROOT
from additional_fluxes import riemann_sample, wave

SUPPORTED = ['first','muscl','minmod','cweno3','cweno5','cweno3-char','cweno5-char']


class LeafMesh(Mesh):
    """Reuse FV moment fits and CWENO polynomials on explicitly supplied faces."""
    def __init__(self, faces):
        self.faces = np.asarray(faces, float).copy()
        w = np.diff(self.faces)
        if len(w)<8 or np.any(w<=0): raise ValueError('Invalid leaf mesh')
        self.n = len(w)
        self.ext_faces = np.r_[self.faces[0]-np.cumsum(w[:3])[::-1], self.faces,
                               self.faces[-1]+np.cumsum(w[-3:][::-1])]
        self.width = np.diff(self.ext_faces)
        self.center = .5*(self.ext_faces[:-1]+self.ext_faces[1:])
        self.fits = {}
        for offsets in [(-1,0,1),(-1,0),(0,1),(-2,-1,0,1,2),(-2,-1,0),(0,1,2)]:
            indices = np.arange(2,self.n+4)[:,None]+np.array(offsets)
            lo = (self.ext_faces[indices]-self.center[2:-2,None])/self.width[2:-2,None]
            hi = (self.ext_faces[indices+1]-self.center[2:-2,None])/self.width[2:-2,None]
            matrix = np.stack([(hi**(k+1)-lo**(k+1))/((k+1)*(hi-lo)) for k in range(len(offsets))],-1)
            self.fits[offsets] = (indices,np.linalg.inv(matrix))


@dataclass
class AMRConfig:
    case: str = 'sod'
    flux: str = 'hllc'
    reconstruction: str = 'cweno3'
    base_cells: int = 80
    max_level: int = 3
    gamma: float = 1.4
    cfl: float = .25
    refine_threshold: float = .02
    coarsen_threshold: float = .004
    buffer_cells: int = 4
    regrid_every: int = 8
    max_steps: int = 100000


def exact_means(faces, case, gamma=1.4, nq=64):
    centers = .5*(faces[:-1]+faces[1:]); widths = np.diff(faces)
    nodes, weights = np.polynomial.legendre.leggauss(nq)
    points = centers[:,None]+widths[:,None]*nodes/2
    q = riemann_sample(CASES[case]['left'],CASES[case]['right'],(points-.5)/CASES[case]['time'],gamma)
    return np.einsum('nqv,q->nv',conserved(q,gamma),weights/2)


def sod_features(gamma=1.4):
    """Reference locations/plateaus for measurement only, never refinement."""
    dl,_,ul,pl = CASES['sod']['left']; dr,_,ur,pr = CASES['sod']['right']
    al,ar = np.sqrt(gamma*pl/dl),np.sqrt(gamma*pr/dr)
    lo,hi = 0.,max(pl,pr)
    for _ in range(70):
        p = (lo+hi)/2
        f = wave(p,dl,pl,al,gamma)+wave(p,dr,pr,ar,gamma)+ur-ul
        if f<0: lo=p
        else: hi=p
    ps = (lo+hi)/2; us = (ul+ur+wave(ps,dr,pr,ar,gamma)-wave(ps,dl,pl,al,gamma))/2
    ratio=ps/pr; z=(gamma-1)/(gamma+1)
    ds = dr*(ratio+z)/(z*ratio+1)
    speed=ur+ar*np.sqrt((gamma+1)/(2*gamma)*ratio+(gamma-1)/(2*gamma))
    t=CASES['sod']['time']
    return {'contact':{'x':.5+us*t,'high':dl*(ps/pl)**(1/gamma),'low':ds},
            'shock':{'x':.5+speed*t,'high':ds,'low':dr}}


def jump_widths(x, rho, gamma=1.4):
    answer = {}
    for name, f in sod_features(gamma).items():
        crossings=[]
        for fraction in [.9,.1]:
            level=f['low']+fraction*(f['high']-f['low'])
            mask=(rho[:-1]>=level)&(rho[1:]<=level)&(abs(.5*(x[:-1]+x[1:])-f['x'])<.065)
            indices=np.flatnonzero(mask)
            if len(indices):
                points=x[indices]+(rho[indices]-level)/(rho[indices]-rho[indices+1])*(x[indices+1]-x[indices])
                crossings.append(float(points[np.argmin(abs(points-f['x']))]))
        answer[name] = {'exact_location':float(f['x']),
                        'width_10_90':max(0.,crossings[1]-crossings[0]) if len(crossings)==2 else None}
    return answer


class AdaptiveTube:
    def __init__(self, c):
        if c.flux not in FLUXES: raise ValueError('AMR supports the shared pointwise fluxes; JST has no nonuniform stencil')
        if c.reconstruction not in SUPPORTED: raise ValueError('Nonuniform reconstruction must be one of '+str(SUPPORTED))
        if c.case not in CASES or c.base_cells<8 or not 0<=c.max_level<=6: raise ValueError('Invalid AMR configuration')
        if not 0<c.coarsen_threshold<c.refine_threshold or c.regrid_every<1 or c.buffer_cells<1: raise ValueError('Invalid refinement controls')
        self.config=c; self.finest=c.base_cells*2**c.max_level
        self.start=np.arange(c.base_cells,dtype=int)*2**c.max_level
        self.span=np.full(c.base_cells,2**c.max_level,dtype=int)
        self.update_mesh()
        fraction=np.clip((.5-self.faces[:-1])/self.dx,0,1)
        self.u=fraction[:,None]*conserved(CASES[c.case]['left'],c.gamma)+(1-fraction[:,None])*conserved(CASES[c.case]['right'],c.gamma)
        self.initial_integral=self.integral();self.budget=np.zeros(4);self.time=0.
        self.events={'rejected_steps':0,'limited_cells':0,'regrids':0,'split_cells':0,'merged_pairs':0,
                     'rhs_cell_updates':0,'max_remap_defect':0.,'cumulative_abs_remap_defect':0.}
        self.mesh_history=[];self.snapshots=[]
        for _ in range(c.max_level): self.adapt()
        self.snapshot()
    def update_mesh(self):
        self.faces=np.r_[self.start/self.finest,1.]
        self.mesh=LeafMesh(self.faces);self.dx=np.diff(self.faces);self.x=.5*(self.faces[:-1]+self.faces[1:])
    def integral(self): return np.einsum('n,nv->v',self.dx,self.u)
    def sensor(self):
        q=primitive(self.u,self.config.gamma)[:,[0,3]]
        jumps=np.max(abs(q[1:]-q[:-1])/(abs(q[1:])+abs(q[:-1])+1e-14),axis=1)
        s=np.maximum(np.r_[0.,jumps],np.r_[jumps,0.])
        padded=np.pad(s,self.config.buffer_cells,mode='edge')
        return np.max(np.lib.stride_tricks.sliding_window_view(padded,2*self.config.buffer_cells+1),axis=1)
    def prolong(self, split):
        ext=np.pad(self.u,((3,3),(0,0)),mode='edge')
        pol=self.mesh.reconstruct(ext,'muscl',self.config.gamma)
        pol,_=limit_admissibility(pol,ext[2:-2],self.config.gamma,[-.5,0,.5])
        left=evaluate(pol,-.25)[1:-1];right=2*self.u-left
        ql,qr=primitive(left,self.config.gamma),primitive(right,self.config.gamma)
        bad=(~np.isfinite(ql).all(axis=1))|(~np.isfinite(qr).all(axis=1))|(ql[:,0]<=1e-12)|(qr[:,0]<=1e-12)|(ql[:,3]<=1e-12)|(qr[:,3]<=1e-12)
        left[bad]=self.u[bad];right[bad]=self.u[bad]
        starts=[];spans=[];states=[]
        for i in range(len(self.u)):
            if split[i]:
                h=self.span[i]//2
                starts.extend([self.start[i],self.start[i]+h]);spans.extend([h,h]);states.extend([left[i],right[i]])
            else: starts.append(self.start[i]);spans.append(self.span[i]);states.append(self.u[i])
        self.events['split_cells']+=int(sum(split))
        self.start=np.array(starts);self.span=np.array(spans);self.u=np.array(states);self.update_mesh()
    def adapt(self):
        if not self.config.max_level: return
        before=self.integral(); sensor=self.sensor(); old_start=self.start.copy();old_span=self.span.copy()
        # Coarsen only exact sibling pairs, outside the buffered wave region.
        starts=[];spans=[];states=[];i=0;merged=0
        while i<len(self.u):
            h=self.span[i]
            merge=(i+1<len(self.u) and h<2**self.config.max_level and self.span[i+1]==h
                   and self.start[i]%(2*h)==0 and max(sensor[i],sensor[i+1])<self.config.coarsen_threshold)
            if merge:
                starts.append(self.start[i]);spans.append(2*h);states.append((self.u[i]+self.u[i+1])/2);i+=2;merged+=1
            else: starts.append(self.start[i]);spans.append(h);states.append(self.u[i]);i+=1
        if merged:
            self.start=np.array(starts);self.span=np.array(spans);self.u=np.array(states);self.update_mesh()
            self.events['merged_pairs']+=merged
        split=(self.sensor()>self.config.refine_threshold)&(self.span>1)
        if np.any(split): self.prolong(split)
        # Enforce adjacent width ratio <=2, including cells created by merging.
        while True:
            split=np.zeros(len(self.span),bool)
            split[:-1]|=self.span[:-1]>2*self.span[1:];split[1:]|=self.span[1:]>2*self.span[:-1]
            if not np.any(split): break
            self.prolong(split)
        defect=float(np.max(abs(self.integral()-before)))
        self.events['max_remap_defect']=max(self.events['max_remap_defect'],defect)
        self.events['cumulative_abs_remap_defect']+=defect
        if not admissible(self.u,self.config.gamma): raise ArithmeticError('Remap means inadmissible')
        if not(np.array_equal(old_start,self.start) and np.array_equal(old_span,self.span)): self.events['regrids']+=1
        self.mesh_history.append([self.time,len(self.u),float(self.dx.min()),defect])
    def rhs(self,u):
        ext=np.pad(u,((3,3),(0,0)),mode='edge')
        pol=self.mesh.reconstruct(ext,self.config.reconstruction,self.config.gamma)
        pol,count=limit_admissibility(pol,ext[2:-2],self.config.gamma,[-.5,-.25,0,.25,.5])
        self.events['limited_cells']+=count;self.events['rhs_cell_updates']+=len(u)
        f=numerical_flux(self.config.flux,evaluate(pol,.5)[:-1],evaluate(pol,-.5)[1:],self.config.gamma,mach_ref=.1)
        return -np.diff(f,axis=0)/self.dx[:,None],f[0]-f[-1]
    def valid(self,u):
        if not admissible(u,self.config.gamma): raise ArithmeticError('Inadmissible RK stage')
        return u
    def step(self,u,dt):
        r0,b0=self.rhs(u);a=self.valid(u+dt*r0)
        r1,b1=self.rhs(a);b=self.valid(.75*u+.25*(a+dt*r1))
        r2,b2=self.rhs(b);v=self.valid(u/3+2*(b+dt*r2)/3)
        return v,dt*(b0+b1+4*b2)/6
    def snapshot(self):
        self.snapshots.append((self.time,self.faces.copy(),primitive(self.u,self.config.gamma)[:,0].copy()))
    def run(self):
        start=perf_counter();end=CASES[self.config.case]['time'];history=[];failed=None;steps=0;next_picture=end/4
        peak=len(self.u);cell_sum=0
        for it in range(self.config.max_steps):
            if self.time>=end-1e-14: break
            q=primitive(self.u,self.config.gamma)
            dt=min(self.config.cfl*float(np.min(self.dx/(abs(q[:,2])+np.sqrt(self.config.gamma*q[:,3]/q[:,0])))),end-self.time)
            for attempt in range(14):
                try: new,increment=self.step(self.u,dt);break
                except (ArithmeticError,FloatingPointError): dt*=.5;self.events['rejected_steps']+=1
            else: failed='Positivity recovery exhausted';break
            if dt<1e-14: failed='Time step collapsed';break
            cell_sum+=len(self.u);self.u=new;self.budget+=increment;self.time+=dt;steps+=1
            if steps%self.config.regrid_every==0 and self.time<end-1e-14: self.adapt()
            peak=max(peak,len(self.u))
            if self.time>=next_picture-1e-14:
                self.snapshot();next_picture+=end/4
            if steps%20==0 or self.time>=end-1e-14:
                q=primitive(self.u,self.config.gamma)
                history.append([self.time,float(q[:,0].min()),float(q[:,3].min()),float(np.max(abs(self.integral()-self.initial_integral-self.budget))),len(self.u)])
        else: failed='Step budget exhausted'
        q=primitive(self.u,self.config.gamma)
        exact=primitive(exact_means(self.faces,self.config.case,self.config.gamma),self.config.gamma)
        check=primitive(exact_means(self.faces,self.config.case,self.config.gamma,128),self.config.gamma)
        reached=bool(self.time>=end-1e-12);errors={}
        if reached:
            for name,k in [('density',0),('velocity',2),('pressure',3)]:
                d=abs(q[:,k]-exact[:,k]);errors[name]={'L1':float(self.dx@d),'L2':float(np.sqrt(self.dx@(d*d))),'Linf':float(d.max())}
        defect=float(np.max(abs(self.integral()-self.initial_integral-self.budget)))
        r={'config':asdict(self.config),'completed':reached,'integrity_passed':bool(reached and admissible(self.u,self.config.gamma) and defect<1e-10),
           'termination':failed or 'Final physical time reached','time':self.time,'steps':steps,'elapsed_seconds':perf_counter()-start,
           'active_cells':len(self.u),'peak_cells':peak,'mean_active_cells':cell_sum/max(steps,1),'min_cell_width':float(self.dx.min()),
           'conservation_budget_defect':defect,'events':self.events,'errors':errors,
           'reference_quadrature_change':float(self.dx@np.mean(abs(exact-check),axis=1)),
           'jump_widths':jump_widths(self.x,q[:,0],self.config.gamma) if self.config.case=='sod' and reached else {},
           'scope':'Nonuniform FV CWENO/MUSCL; global stepping, 2:1 leaf mesh, conservative regridding; no DG or nonuniform WENO/TENO/JST',
           'faces':self.faces,'x':self.x,'primitive':q,'conserved_average':self.u,'exact':exact,'history':np.array(history),
           'mesh_history':np.array(self.mesh_history).reshape(-1,4)}
        for i,(t,faces,rho) in enumerate(self.snapshots): r[f'snapshot_{i}_faces']=faces;r[f'snapshot_{i}_density']=rho;r[f'snapshot_{i}_time']=np.array(t)
        return r


def run_case(flux='hllc',reconstruction='cweno3',base_cells=80,max_level=3,case='sod',output_dir=None):
    solver=AdaptiveTube(AMRConfig(flux=flux,reconstruction=reconstruction,base_cells=base_cells,max_level=max_level,case=case))
    r=solver.run();target=Path(output_dir or ROOT/'results-amr');target.mkdir(parents=True,exist_ok=True)
    name=f'{case}-{flux}-{reconstruction}-base{base_cells}-level{max_level}'
    arrays={k:v for k,v in r.items() if isinstance(v,np.ndarray)};meta={k:v for k,v in r.items() if k not in arrays}
    np.savez_compressed(target/(name+'.npz'),**arrays)
    (target/(name+'.json')).write_text(json.dumps(meta,indent=2)+'\n')
    return r


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--flux',choices=list(FLUXES),default='hllc')
    p.add_argument('--reconstruction',choices=SUPPORTED,default='cweno3');p.add_argument('--base-cells',type=int,default=80)
    p.add_argument('--max-level',type=int,default=3);p.add_argument('--case',choices=list(CASES),default='sod')
    args=p.parse_args();r=run_case(**vars(args));print(json.dumps({k:v for k,v in r.items() if not isinstance(v,np.ndarray)},indent=2))
