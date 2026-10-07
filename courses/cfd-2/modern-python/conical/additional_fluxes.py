"""Independent FV Euler fluxes. Tangential velocity is advected across contact.

Godunov samples the exact ideal-gas Riemann solution at x/t=0, including vacuum.
Global LF is the semidiscrete flux with one maximum spectral radius per batch.
Roe without entropy smoothing is retained as an educational comparison.
"""
import numpy as np
from state import primitive, conserved, physical_flux


def wave(p, rho, pk, a, gamma):
    shock = (p-pk)*np.sqrt(2/((gamma+1)*rho)/(p+(gamma-1)/(gamma+1)*pk))
    rare = 2*a/(gamma-1)*((np.maximum(p,0)/pk)**((gamma-1)/(2*gamma))-1)
    return np.where(p > pk, shock, rare)


def riemann_sample(left, right, xi=0., gamma=1.4):
    """Primitive exact solution for independently broadcast Riemann problems."""
    shape=np.broadcast_shapes(np.shape(left)[:-1],np.shape(right)[:-1],np.shape(xi))
    l=np.broadcast_to(np.asarray(left,float),shape+(4,))
    r=np.broadcast_to(np.asarray(right,float),shape+(4,))
    dl, tl, ul, pl = np.moveaxis(l,-1,0)
    dr, tr, ur, pr = np.moveaxis(r,-1,0)
    al, ar = np.sqrt(gamma*pl/dl),np.sqrt(gamma*pr/dr)
    xi = np.broadcast_to(xi,dl.shape)
    vacuum = ur-ul >= 2*(al+ar)/(gamma-1)
    lo = np.zeros_like(dl)
    hi = np.maximum(pl,pr)
    for _ in range(24):
        f = wave(hi,dl,pl,al,gamma)+wave(hi,dr,pr,ar,gamma)+ur-ul
        hi = np.where(f<0,2*hi,hi)
    for _ in range(52):
        p = (lo+hi)/2
        f = wave(p,dl,pl,al,gamma)+wave(p,dr,pr,ar,gamma)+ur-ul
        lo,hi = np.where(f<0,p,lo),np.where(f>=0,p,hi)
    ps = (lo+hi)/2
    us = (ul+ur+wave(ps,dr,pr,ar,gamma)-wave(ps,dl,pl,al,gamma))/2
    side = xi <= us
    # Reflect the right-side solution into left coordinates; then one fan formula.
    d,t,u,p,a = [np.where(side,x,y) for x,y in
                [(dl,dr),(tl,tr),(ul,-ur),(pl,pr),(al,ar)]]
    s = np.where(side,xi,-xi)
    vstar = np.where(side,us,-us)
    ratio = ps/p
    shock_speed = u-a*np.sqrt((gamma+1)/(2*gamma)*ratio+(gamma-1)/(2*gamma))
    ds_shock = d*(ratio+(gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*ratio+1)
    ds_rare = d*ratio**(1/gamma)
    tail = vstar-a*ratio**((gamma-1)/(2*gamma))
    head = u-a
    vf = 2/(gamma+1)*(a+(gamma-1)*u/2+s)
    af = np.maximum(2/(gamma+1)*(a+(gamma-1)*(u-s)/2),0)
    df = d*(af/a)**(2/(gamma-1))
    pf = p*(af/a)**(2*gamma/(gamma-1))
    initial = np.where(ps>p,s<=shock_speed,s<=head)
    fan = (ps<=p)&(s>head)&(s<tail)
    density = np.where(initial,d,np.where(fan,df,np.where(ps>p,ds_shock,ds_rare)))
    pressure = np.where(initial,p,np.where(fan,pf,ps))
    velocity = np.where(initial,u,np.where(fan,vf,vstar))*np.where(side,1,-1)
    # Vacuum consists of the two independent rarefaction fans and an empty gap.
    vl_end,vr_end = ul+2*al/(gamma-1),ur-2*ar/(gamma-1)
    vleft = 2/(gamma+1)*(al+(gamma-1)*ul/2+xi)
    cleft = np.maximum(2/(gamma+1)*(al+(gamma-1)*(ul-xi)/2),0)
    vright = 2/(gamma+1)*(-ar+(gamma-1)*ur/2+xi)
    cright = np.maximum(2/(gamma+1)*(ar+(gamma-1)*(xi-ur)/2),0)
    vd = np.where(xi<ul-al,dl,np.where(xi<vl_end,dl*(cleft/al)**(2/(gamma-1)),
         np.where(xi<=vr_end,0,np.where(xi<ur+ar,dr*(cright/ar)**(2/(gamma-1)),dr))))
    vp = np.where(xi<ul-al,pl,np.where(xi<vl_end,pl*(cleft/al)**(2*gamma/(gamma-1)),
         np.where(xi<=vr_end,0,np.where(xi<ur+ar,pr*(cright/ar)**(2*gamma/(gamma-1)),pr))))
    vv = np.where(xi<ul-al,ul,np.where(xi<vl_end,vleft,np.where(xi<=vr_end,xi,np.where(xi<ur+ar,vright,ur))))
    density,pressure,velocity = [np.where(vacuum,v,x) for v,x in [(vd,density),(vp,pressure),(vv,velocity)]]
    return np.stack((density,t,velocity,pressure),-1)


def godunov(left,right,gamma=1.4):
    q = riemann_sample(primitive(left,gamma),primitive(right,gamma),gamma=gamma)
    rho,t,v,p = np.moveaxis(q,-1,0)
    energy = p/(gamma-1)+rho*(t*t+v*v)/2
    return np.stack((rho*v,rho*t*v,rho*v*v+p,(energy+p)*v),-1)


def roe_uncorrected(left,right,gamma=1.4):
    from flux import roe
    return roe(left,right,gamma,entropy_fix=False)


def global_lf(left,right,gamma=1.4):
    l,r = primitive(left,gamma),primitive(right,gamma)
    speed = max(float(np.max(abs(l[...,2])+np.sqrt(gamma*l[...,3]/l[...,0]))),
                float(np.max(abs(r[...,2])+np.sqrt(gamma*r[...,3]/r[...,0]))))
    return (physical_flux(left,gamma)+physical_flux(right,gamma)-speed*(right-left))/2


def steger_warming(left,right,gamma=1.4):
    def split(u,sign):
        q = primitive(u,gamma)
        rho,t,v,p = np.moveaxis(q,-1,0)
        a = np.sqrt(gamma*p/rho)
        h = (u[...,3]+p)/rho
        one = np.ones_like(rho)
        acoustic_minus = np.stack((one,t,v-a,h-v*a),-1)
        acoustic_plus = np.stack((one,t,v+a,h+v*a),-1)
        contact = np.stack((one,t,v,(t*t+v*v)/2),-1)
        def lam(x): return (x+sign*abs(x))/2
        return rho[...,None]/(2*gamma)*(lam(v-a)[...,None]*acoustic_minus+
            lam(v+a)[...,None]*acoustic_plus+2*(gamma-1)*lam(v)[...,None]*contact)
    return split(left,1)+split(right,-1)


def jst_flux(ext,gamma=1.4,k2=.5,k4=.02):
    """FV central flux with pressure-sensed second/fourth artificial differences.

    ext has three ghost cells on each side. This is a stencil scheme, not a
    pointwise Riemann solver or a flux to compose with high-order reconstruction.
    """
    q = primitive(ext,gamma)
    p = q[:,3]
    sensor = abs(p[2:]-2*p[1:-1]+p[:-2])/(abs(p[2:])+2*abs(p[1:-1])+abs(p[:-2])+1e-30)
    nu = np.maximum(sensor[1:-2],sensor[2:-1])
    l,r = ext[2:-3],ext[3:-2]
    speed = abs(q[:,2])+np.sqrt(gamma*p/q[:,0])
    scale = np.maximum(speed[2:-3],speed[3:-2])
    eps2 = k2*nu
    eps4 = np.maximum(0,k4-eps2)
    d3 = ext[4:-1]-3*r+3*l-ext[1:-4]
    return physical_flux((l+r)/2,gamma)-scale[:,None]*(eps2[:,None]*(r-l)-eps4[:,None]*d3)
