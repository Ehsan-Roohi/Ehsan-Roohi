"""Independent modal RKDG companion for the same angular Euler cone equations.

Degree 0/1/2 Legendre polynomials, overintegrated weak form, shared Euler
fluxes, characteristic TVB slope limiting and sampled positivity scaling.
This is a one-dimensional angular DG implementation, not unstructured 2-D DG,
not a viscous DG method, and not an implementation of a separate FR scheme.
"""
from dataclasses import dataclass, asdict
from time import perf_counter
import argparse
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import legvander, Legendre
from state import primitive, physical_flux, source, freestream, admissible
from flux import FLUXES, numerical_flux
from characteristics import basis
from solver import Config, save_result


@dataclass
class DGConfig(Config):
    degree: int = 1
    cfl: float = .25
    cells: int = 40
    mach: float = 2.35
    outer_deg: float = 45
    max_steps: int = 20000
    reconstruction: str = 'modal-dg'
    tvb: float = 10


def minmod(a, b, c):
    same = (np.sign(a) == np.sign(b)) & (np.sign(a) == np.sign(c))
    return np.where(same, np.sign(a)*np.minimum(abs(a), np.minimum(abs(b), abs(c))), 0)


class DGConicalSolver:
    def __init__(self, config):
        c = self.config = config
        if c.degree not in (0, 1, 2) or c.cells < 4:
            raise ValueError('DG supports degree 0, 1 or 2 and at least four cells')
        if c.balance_freestream:
            raise ValueError('DG equilibrium subtraction is not implemented')
        self.h = np.radians(c.outer_deg-c.cone_deg)/c.cells
        self.theta = np.radians(c.cone_deg)+(np.arange(c.cells)+.5)*self.h
        self.nodes, self.weights = np.polynomial.legendre.leggauss(max(3, c.degree+2))
        self.v = legvander(self.nodes, c.degree)
        self.dv = np.column_stack([Legendre.basis(k).deriv()(self.nodes) for k in range(c.degree+1)])
        self.ends = legvander(np.array([-1., 1.]), c.degree)
        self.sample = legvander(np.r_[-1., self.nodes, 1.], c.degree)
        self.scale = 2*np.arange(c.degree+1)+1
        self.u = self.project(freestream(self.theta[:, None]+self.h*self.nodes/2, c.mach, c.gamma))
        self.events = {'slope_limited_cells': 0, 'positivity_scaled_cells': 0, 'rejected_steps': 0}

    def project(self, values):
        return np.einsum('nqv,qk,q,k->nkv', values, self.v, self.weights, self.scale/2)

    def values(self, modes, matrix=None):
        return np.einsum('nkv,qk->nqv', modes, self.v if matrix is None else matrix)

    def limit(self, modes):
        c = self.config
        result = modes.copy()
        mean = result[:, 0]
        if not admissible(mean, c.gamma):
            raise ArithmeticError('DG cell mean is inadmissible')
        if c.degree:
            prev = np.r_[mean[:1], mean[:-1]].copy()
            prev[0, 2] *= -1
            nxt = np.r_[mean[1:], freestream(np.array([np.radians(c.outer_deg)+self.h/2]), c.mach, c.gamma)]
            r, l = basis(mean, c.gamma)
            a = np.einsum('nij,nj->ni', l, result[:, 1])
            b = np.einsum('nij,nj->ni', l, mean-prev)
            d = np.einsum('nij,nj->ni', l, nxt-mean)
            limited = np.where(abs(a) <= c.tvb*self.h**2, a, minmod(a, b, d))
            changed = np.max(abs(limited-a), axis=1) > 1e-12
            result[changed, 1] = np.einsum('nij,nj->ni', r[changed], limited[changed])
            result[changed, 2:] = 0
            self.events['slope_limited_cells'] += int(np.sum(changed))
        # Scaling preserves every conserved cell mean. It is a sampled safeguard,
        # not a theorem of positivity for arbitrary time steps or all interior x.
        samples = self.values(result, self.sample)
        bad = (~np.all(np.isfinite(samples), axis=(1, 2))) | (np.min(samples[..., 0], axis=1) <= 1e-12) | (np.min(primitive(samples, c.gamma)[..., 3], axis=1) <= 1e-12)
        for i in np.flatnonzero(bad):
            lo, hi = 0., 1.
            for _ in range(45):
                mid = (lo+hi)/2
                if admissible(mean[i]+mid*(samples[i]-mean[i]), c.gamma):
                    lo = mid
                else:
                    hi = mid
            result[i, 1:] *= lo*(1-1e-12)
            self.events['positivity_scaled_cells'] += 1
        return result

    def rhs(self, modes):
        c = self.config
        values = self.values(modes)
        ends = self.values(modes, self.ends)
        left = np.r_[ends[:1, 0], ends[:, 1]]
        left[0, 2] *= -1
        right = np.r_[ends[:, 0], freestream(np.array([np.radians(c.outer_deg)]), c.mach, c.gamma)]
        face = numerical_flux(c.flux, left, right, c.gamma, mach_ref=c.mach)
        if c.boundary == 'wall':
            face[0] = [0, 0, primitive(ends[0, 0], c.gamma)[3], 0]
        f = physical_flux(values, c.gamma)
        s = source(values, self.theta[:, None]+self.h*self.nodes/2, c.gamma)
        volume = np.einsum('nqv,qk,q->nkv', f, self.dv, self.weights)
        boundary = face[1:, None]*self.ends[1, :, None]-face[:-1, None]*self.ends[0, :, None]
        src = np.einsum('nqv,qk,q->nkv', s, self.v, self.weights)*self.scale[None, :, None]/2
        residual = (volume-boundary)*self.scale[None, :, None]/self.h+src
        self.balance_error = float(np.max(abs(np.sum(residual[:, 0]*self.h, axis=0)+face[-1]-face[0]-np.sum(src[:, 0]*self.h, axis=0))))
        return residual

    def time_step(self, modes):
        q = primitive(self.values(modes, self.sample), self.config.gamma)
        speed = np.max(abs(q[..., 2])+np.sqrt(self.config.gamma*q[..., 3]/q[..., 0]))
        return self.config.cfl*self.h/((2*self.config.degree+1)*speed)

    def step(self, modes, dt):
        a = self.limit(modes+dt*self.rhs(modes))
        b = self.limit(.75*modes+.25*(a+dt*self.rhs(a)))
        return self.limit(modes/3+2*(b+dt*self.rhs(b))/3)

    def run(self, verbose=False):
        c = self.config
        start = perf_counter()
        history = []
        self.u = self.limit(self.u)
        converged = False
        drift = None
        for iteration in range(1, c.max_steps+1):
            old = self.u
            dt = self.time_step(old)
            for attempt in range(13):
                try:
                    self.u = self.step(old, dt)
                    break
                except (ArithmeticError, FloatingPointError):
                    self.events['rejected_steps'] += 1
                    dt *= .5
            else:
                raise RuntimeError('DG positivity recovery exhausted')
            if iteration % 50 == 0 or iteration == 1 or iteration == c.max_steps:
                # Includes all modes, not only cell means; shock limiting may
                # leave a nonzero semidiscrete residual at a fixed time map.
                residual = float(np.max(np.sqrt(np.mean(self.rhs(self.u)**2, axis=0))))
                change = float(np.max(abs(self.u-old)))
                wall = primitive(self.values(self.u[:1], self.ends)[0, 0], c.gamma)
                history.append([iteration, residual, wall[3]*c.gamma*c.mach**2, change])
                if verbose and iteration % 1000 == 0:
                    print(f'DG p={c.degree} {c.flux} N={c.cells} step={iteration} R={residual:.3e}', flush=True)
                if residual < c.tolerance:
                    check = self.u.copy()
                    for _ in range(20):
                        check = self.step(check, self.time_step(check))
                    drift = float(np.max(abs(check-self.u)))
                    converged = drift < 1e-7
                    if converged:
                        break
        q = primitive(self.values(self.u, legvander(np.array([0.]), c.degree))[:, 0], c.gamma)
        wall = primitive(self.values(self.u[:1], self.ends)[0, 0], c.gamma)
        return {'config': asdict(c), 'converged': converged, 'iterations': iteration,
                'residual_rms_max_equation': residual, 'elapsed_seconds': perf_counter()-start,
                'wall_pressure_ratio': float(wall[3]*c.gamma*c.mach**2),
                'wall_mach': float(np.hypot(wall[1], wall[2])/np.sqrt(c.gamma*wall[3]/wall[0])),
                'min_density': float(q[:, 0].min()), 'min_pressure': float(q[:, 3].min()),
                'discrete_balance_roundoff': self.balance_error, 'events': self.events.copy(),
                'explicit_20_step_max_state_drift': drift,
                'scope': 'Inviscid angular modal RKDG companion; not viscous or 2-D unstructured DG',
                'history': history, 'theta': self.theta, 'primitive': q,
                'conserved_average': self.u[:, 0].copy(), 'modal_coefficients': self.u.copy()}


def save_dg(result, path):
    save_result(result | {'modal_coefficients': 'stored separately in modal NPZ'}, path)
    np.savez_compressed(str(path)+'-modes.npz', modes=result['modal_coefficients'])


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--degree', type=int, choices=[0, 1, 2], default=1)
    p.add_argument('--flux', choices=FLUXES, default='hllc')
    p.add_argument('--cells', type=int, default=40)
    p.add_argument('--steps', type=int, default=20000)
    p.add_argument('--mach', type=float, default=2.35)
    p.add_argument('--outer', type=float, default=45)
    p.add_argument('--output', default='results-dg/single')
    a = p.parse_args()
    result = DGConicalSolver(DGConfig(degree=a.degree, flux=a.flux, cells=a.cells,
                                    max_steps=a.steps, mach=a.mach, outer_deg=a.outer)).run(True)
    save_dg(result, a.output)
    print(json.dumps({k:v for k,v in result.items() if k not in ('history','theta','primitive','conserved_average','modal_coefficients')}, indent=2))
