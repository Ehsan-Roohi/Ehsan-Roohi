"""Independent modal-DG checks: FV reduction, conservation, smooth order, FR link."""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import legvander, Legendre
from dg import DGConfig, DGConicalSolver
from solver import Config, ConicalSolver
from state import admissible

ROOT = Path(__file__).resolve().parent


def advection_rhs(modes, h, nodes, weights, v, dv):
    p = modes.shape[1]-1
    ends = legvander(np.array([-1., 1.]), p)
    values = modes @ v.T
    face_right = modes @ ends[1]
    face_left = np.roll(face_right, 1)
    vol = (values*weights) @ dv
    return (vol-face_right[:, None]*ends[1]+face_left[:, None]*ends[0])*(2*np.arange(p+1)+1)/h


def smooth_order(degree):
    rows = []
    nodes, weights = np.polynomial.legendre.leggauss(degree+3)
    v = legvander(nodes, degree)
    dv = np.column_stack([Legendre.basis(k).deriv()(nodes) for k in range(degree+1)])
    for n in [12, 24, 48]:
        h = 1/n
        centers = (np.arange(n)+.5)*h
        points = centers[:, None]+h*nodes/2
        u = ((np.sin(2*np.pi*points)*weights) @ v)*(2*np.arange(degree+1)+1)/2
        t = .15
        steps = int(np.ceil(t/(.15*h/(2*degree+1))))
        dt = t/steps
        for _ in range(steps):
            a = u+dt*advection_rhs(u,h,nodes,weights,v,dv)
            b = .75*u+.25*(a+dt*advection_rhs(a,h,nodes,weights,v,dv))
            u = u/3+2*(b+dt*advection_rhs(b,h,nodes,weights,v,dv))/3
        error = float(np.sqrt(np.sum((u@v.T-np.sin(2*np.pi*(points-t)))**2*weights)/(2*n)))
        rows.append({'cells':n, 'l2_error':error, 'order':None if not rows else float(np.log2(rows[-1]['l2_error']/error))})
    return rows


def fr_equivalence(degree):
    nodes, weights = np.polynomial.legendre.leggauss(degree+1)
    v = legvander(nodes,degree)
    dv = np.column_stack([Legendre.basis(k).deriv()(nodes) for k in range(degree+1)])
    modes = np.random.default_rng(19).normal(size=(8,degree+1))
    ends = legvander(np.array([-1.,1.]),degree)
    f_left, f_right = (modes@ends.T).T
    numerical_left = np.roll(f_right,1)
    gleft = ((-1)**degree/2*(Legendre.basis(degree)-Legendre.basis(degree+1))).deriv()(nodes)
    gright = (.5*(Legendre.basis(degree)+Legendre.basis(degree+1))).deriv()(nodes)
    # Right upwind correction is zero for positive unit advection velocity.
    fr = -2*(modes@dv.T+(numerical_left-f_left)[:,None]*gleft+0*gright)
    dg = advection_rhs(modes,1,nodes,weights,v,dv)@v.T
    return float(np.max(abs(fr-dg)))


def verify():
    cfg = DGConfig(degree=0,cells=30,mach=2.35,outer_deg=45)
    dg = DGConicalSolver(cfg)
    fv = ConicalSolver(Config(cells=30,mach=2.35,outer_deg=45,reconstruction='first'))
    fv.u = dg.u[:,0].copy()
    a, b = dg.rhs(dg.u)[:,0], fv.rhs(fv.u,record=False)
    # FV outer ghosts are cell averages; DG exterior is a point boundary value.
    # The interior and wall operators must match; the last cell is excluded.
    reduction = float(np.max(abs(a[:-1]-b[:-1])))
    degree = DGConicalSolver(DGConfig(degree=2,cells=30))
    degree.rhs(degree.u)
    conservation = degree.balance_error
    original = degree.u.copy()
    original[:,1:] *= 50
    limited = degree.limit(original)
    limiter_mean_error = float(np.max(abs(limited[:,0]-original[:,0])))
    positive = bool(admissible(degree.values(limited,degree.sample)))
    order = {str(p):smooth_order(p) for p in [1,2]}
    fr = {str(p):fr_equivalence(p) for p in [0,1,2]}
    result = {'fv_degree_zero_interior_wall_error':reduction,
              'integrated_mean_balance_error':conservation,
              'limiter_mean_preservation_error':limiter_mean_error,
              'limited_samples_positive':positive, 'periodic_advection_order':order,
              'linear_fr_radau_dg_equivalence_error':fr,
              'scope':'Linear smooth time-dependent tests establish DG order; no high-order shock or viscous DG claim.'}
    result['passed'] = reduction < 1e-12 and conservation < 1e-12 and limiter_mean_error < 1e-14 and positive and max(fr.values()) < 1e-11 and order['1'][-1]['order'] > 1.8 and order['2'][-1]['order'] > 2.7
    return result


if __name__ == '__main__':
    result = verify()
    out = ROOT/'results-dg'
    out.mkdir(exist_ok=True)
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    if not result['passed']:
        raise SystemExit('DG verification failed')
