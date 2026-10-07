"""Independent Sod shock-tube benchmark for the shared Euler face fluxes.

Exact ideal-gas Riemann solution; conservative FV CWENO3 and SSPRK3.
The planar benchmark tests shock/contact/rarefaction handling separately from
conical geometry. No external CFD solver is invoked.
"""
from pathlib import Path
import json
import numpy as np
from state import conserved, primitive, admissible
from flux import numerical_flux
from reconstruction import coefficients, evaluate, limit_admissibility

ROOT = Path(__file__).resolve().parent
METHODS = ['roe', 'hllc', 'ausm', 'ausm-plus', 'ausm-up', 'ausm-up2', 'slau2', 'ec-lf']
GAMMA = 1.4


def exact(x, time=.2):
    g = GAMMA
    dl, dr, pl, pr = 1., .125, 1., .1
    cl, cr = np.sqrt(g*pl/dl), np.sqrt(g*pr/dr)
    def wave(p, d, pk, ck):
        if p > pk:
            return (p-pk)*np.sqrt(2/((g+1)*d)/(p+(g-1)/(g+1)*pk))
        return 2*ck/(g-1)*((p/pk)**((g-1)/(2*g))-1)
    lo, hi = pr, pl
    for _ in range(80):
        mid = (lo+hi)/2
        if wave(mid,dl,pl,cl)+wave(mid,dr,pr,cr) > 0:
            hi = mid
        else:
            lo = mid
    ps = (lo+hi)/2
    us = .5*(wave(ps,dr,pr,cr)-wave(ps,dl,pl,cl))
    sl = us-cl*(ps/pl)**((g-1)/(2*g))
    sr = cr*np.sqrt((g+1)/(2*g)*ps/pr+(g-1)/(2*g))
    xi = (x-.5)/time
    q = np.zeros((len(x),4))
    q[:, 0], q[:, 3] = dr, pr
    mask = xi < -cl
    q[mask, 0], q[mask, 3] = dl, pl
    mask = (xi >= -cl) & (xi < sl)
    vel = 2/(g+1)*(cl+xi[mask])
    cs = 2/(g+1)*(cl-(g-1)*xi[mask]/2)
    q[mask,0],q[mask,2],q[mask,3] = dl*(cs/cl)**(2/(g-1)),vel,pl*(cs/cl)**(2*g/(g-1))
    mask = (xi >= sl) & (xi < us)
    q[mask,0],q[mask,2],q[mask,3] = dl*(ps/pl)**(1/g),us,ps
    mask = (xi >= us) & (xi < sr)
    ratio = ps/pr
    q[mask,0],q[mask,2],q[mask,3] = dr*(ratio+(g-1)/(g+1))/((g-1)/(g+1)*ratio+1),us,ps
    return q, ps, us


def run(n=200):
    x = (np.arange(n)+.5)/n
    q0 = np.zeros((n,4))
    q0[:,0] = np.where(x < .5,1,.125)
    q0[:,3] = np.where(x < .5,1,.1)
    exact_q, ps, us = exact(x)
    out = ROOT/'results-v2'
    out.mkdir(exist_ok=True)
    rows = []
    for name in METHODS:
        u = conserved(q0)
        def rhs(state):
            ext = np.pad(state,((3,3),(0,0)),mode='edge')
            pol = coefficients(ext,'cweno3',1/n)
            pol,_ = limit_admissibility(pol,ext[2:-2],GAMMA,[-.5,0,.5])
            f = numerical_flux(name,evaluate(pol,.5)[:-1],evaluate(pol,-.5)[1:],mach_ref=.1)
            return -np.diff(f,axis=0)*n
        time, steps, rejected = 0.,0,0
        while time < .2-1e-14:
            q = primitive(u)
            dt = min(.25/n/np.max(abs(q[:,2])+np.sqrt(GAMMA*q[:,3]/q[:,0])),.2-time)
            for attempt in range(13):
                a = u+dt*rhs(u)
                if admissible(a):
                    b = .75*u+.25*(a+dt*rhs(a))
                    if admissible(b):
                        c = u/3+2/3*(b+dt*rhs(b))
                        if admissible(c):
                            break
                rejected += 1
                dt *= .5
            else:
                raise RuntimeError(f'{name}: rejected Sod step')
            u = c
            time += dt
            steps += 1
        q = primitive(u)
        rows.append({'flux':name,'cells':n,'time':time,'steps':steps,'rejected':rejected,
                     'density_l1':float(np.mean(abs(q[:,0]-exact_q[:,0]))),
                     'pressure_l1':float(np.mean(abs(q[:,3]-exact_q[:,3]))),
                     'min_density':float(q[:,0].min()),'min_pressure':float(q[:,3].min())})
        np.savez_compressed(out/f'sod-{name}.npz',x=x,primitive=q,exact=exact_q)
    assert all(r['density_l1'] < .02 and r['pressure_l1'] < .02 for r in rows)
    result = {'passed':True,'reference':'Exact ideal-gas Riemann solution by pressure bisection',
              'star_pressure':ps,'star_velocity':us,'runs':rows}
    (out/'sod-verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    run()
