"""Local radial-station viscous angular cone model, not a full 2-D NS solver.

The archived assumption d/dr=0 for primitive variables is retained. Stresses
still scale as mu/r; their radial geometric derivatives are included. Re is
local to the chosen radial station. The model cannot predict downstream
boundary-layer growth or replace a full axisymmetric cone validation.
"""
from dataclasses import dataclass
import argparse
import numpy as np
from solver import Config, ConicalSolver, save_result
from state import primitive, source, freestream_average, conserved
from reconstruction import evaluate, limit_admissibility
from flux import FLUXES, numerical_flux
from fv_mesh import Mesh


@dataclass
class ViscousConfig(Config):
    reconstruction: str = 'primitive-minmod'
    viscous_operator_version: str = 'fv-wall-quadratic-v2'
    reynolds: float = 420000
    prandtl: float = .72
    temperature_inf_K: float = 775.56/(1+.2*7.95**2)
    viscosity: str = 'sutherland'
    thermal_wall: str = 'adiabatic'
    wall_temperature_ratio: float = 1
    stretch: float = 3
    gradient_order: int = 2
    cfl: float = .2
    max_steps: int = 80000


def viscosity_ratio(temperature, config):
    if config.viscosity == 'constant':
        return np.ones_like(temperature)
    return temperature**1.5*(config.temperature_inf_K+110.4)/(temperature*config.temperature_inf_K+110.4)


def transport(q, derivative, theta, c):
    """q=[rho,vr,vt,T/Tinf]; derivative is with respect to theta, r=1."""
    rho, vr, vt, temp = np.moveaxis(q, -1, 0)
    dvr, dvt, dtemp = derivative[..., 1], derivative[..., 2], derivative[..., 3]
    mu = viscosity_ratio(temp, c)/c.reynolds
    cot = 1/np.tan(theta)
    rr = -2/3*mu*(dvt+2*vr+vt*cot)
    tt = 2/3*mu*(2*dvt+vr-vt*cot)
    pp = 2/3*mu*(-dvt+vr+2*vt*cot)
    rt = mu*(dvr-vt)
    heat = -mu*dtemp/(c.prandtl*(c.gamma-1)*c.mach**2)
    f = np.stack((np.zeros_like(rr), rt, tt, rt*vr+tt*vt-heat), -1)
    s = np.stack((np.zeros_like(rr), rr-tt-pp+rt*cot,
                  2*rt+(tt-pp)*cot, rr*vr+rt*vt+cot*(rt*vr+tt*vt-heat)), -1)
    return f, s, np.stack((rr, tt, pp, rt, heat), -1)


def thermal_primitive(u, c):
    q = primitive(u, c.gamma)
    return np.concatenate((q[..., :3], (c.gamma*c.mach**2*q[..., 3]/q[..., 0])[..., None]), -1)


class ViscousConicalSolver(ConicalSolver):
    def __init__(self, config):
        c = config
        if c.reynolds <= 0 or c.prandtl <= 0 or c.wall_temperature_ratio <= 0:
            raise ValueError('Re, Pr and wall temperature must be positive')
        if c.balance_freestream:
            raise ValueError('Equilibrium subtraction is not implemented for viscous runs')
        self.config = c
        self.mesh = Mesh(np.radians(c.cone_deg), np.radians(c.outer_deg), c.cells, c.stretch)
        self.h = self.mesh.width[3:-3]
        self.theta = self.mesh.center[3:-3]
        self.ext_theta = self.mesh.center
        self.upstream_ext = freestream_average(self.ext_theta, self.mesh.width, c.mach, c.gamma)
        self.u = self.upstream_ext[3:-3].copy()
        self.nodes, self.weights = np.polynomial.legendre.leggauss(3)
        self.nodes *= .5
        self.weights *= .5
        self.samples = np.r_[-.5, self.nodes, .5]
        self.events = {'reconstruction_limited_cells': 0, 'rejected_steps': 0, 'rhs_evaluations': 0}
        self.correction = np.zeros_like(self.u)
        ratio = self.h[1]/self.h[0]
        self.wall_moments = np.array([[.5, 1/3], [1+ratio/2, ((1+ratio)**3-1)/(3*ratio)]])
        self.wall_dirichlet_fit = np.linalg.inv(self.wall_moments)

    def ghost(self, u, uniform_boundary=False):
        ext = self.upstream_ext.copy()
        ext[3:-3] = u
        if not uniform_boundary:
            q = primitive(u[:3][::-1], self.config.gamma)
            q[:, 1:3] *= -1
            # Convective wall ghosts preserve pressure, density and temperature.
            # Thermal BC is imposed separately in the physical viscous flux and
            # first-cell source polynomial; it cannot affect wall convective heat.
            ext[:3] = conserved(q, self.config.gamma)
        return ext

    def polynomials(self, u, uniform_boundary=False, record=True):
        ext = self.ghost(u, uniform_boundary)
        pol = self.mesh.reconstruct(ext, self.config.reconstruction, self.config.gamma)
        pol, count = limit_admissibility(pol, ext[2:-2], self.config.gamma, self.samples)
        if record:
            self.events['reconstruction_limited_cells'] += count
        return pol

    def viscous_rhs(self, u, pol=None):
        c = self.config
        if pol is None:
            pol = self.polynomials(u, record=False)
        q = thermal_primitive(self.ghost(u), c)
        if c.gradient_order == 4:
            # Deconvolve conservative averages with cell polynomials before
            # averaging primitives. Boundary closure remains second order.
            q[2:-2] = np.einsum('nqk,q->nk', thermal_primitive(evaluate(pol, self.nodes), c), self.weights)
        face_q, face_d = self.mesh.face_value_gradient(q, c.gradient_order)
        first = thermal_primitive(u[:2], c)
        face_q[0, 1:3] = 0
        wall_coeff = np.zeros((3, 4))
        wall_coeff[0] = first[0]
        velocity_fit = self.wall_dirichlet_fit@q[3:5, 1:3]
        wall_coeff[0, 1:3] = .5*velocity_fit[0]+.25*velocity_fit[1]
        wall_coeff[1, 1:3] = velocity_fit[0]+velocity_fit[1]
        wall_coeff[2, 1:3] = velocity_fit[1]
        face_d[0, 1:3] = velocity_fit[0]/self.h[0]
        if c.thermal_wall == 'adiabatic':
            curvature = (q[4, 3]-q[3, 3])/(self.wall_moments[1, 1]-self.wall_moments[0, 1])
            wall_temperature = q[3, 3]-curvature/3
            wall_coeff[:, 3] = [wall_temperature+curvature/4, curvature, curvature]
            face_q[0, 3] = wall_temperature
            face_d[0, 3] = 0
        else:
            b, curvature = self.wall_dirichlet_fit@(q[3:5, 3]-c.wall_temperature_ratio)
            wall_coeff[:, 3] = [c.wall_temperature_ratio+b/2+curvature/4, b+curvature, curvature]
            face_q[0, 3] = c.wall_temperature_ratio
            face_d[0, 3] = b/self.h[0]
        if np.any(face_q[:, 3] <= 0):
            raise ArithmeticError('Non-positive face temperature in diffusion fit')
        fv, _, wall_data = transport(face_q, face_d, self.mesh.faces, c)
        qp = self.mesh.reconstruct(q, 'cweno3' if c.gradient_order == 2 else 'cweno5')[1:-1]
        qp[0] = 0
        qp[0, :3] = wall_coeff
        vals = evaluate(qp, self.nodes)
        deriv = evaluate(qp[:, 1:]*np.arange(1, qp.shape[1])[None, :, None], self.nodes)/self.h[:, None, None]
        if np.any(vals[..., 3] <= 0):
            raise ArithmeticError('Non-positive source temperature')
        _, src, _ = transport(vals, deriv, self.theta[:, None]+self.h[:, None]*self.nodes, c)
        sv = np.einsum('nqk,q->nk', src, self.weights)
        self.last_transport = wall_data[0]
        self.wall_thermal_value = face_q[0, 3]
        return (fv[1:]-fv[:-1])/self.h[:, None]+sv, fv, sv

    def rhs(self, u, uniform_boundary=False, record=True):
        c = self.config
        pol = self.polynomials(u, record=record)
        f = numerical_flux(c.flux, evaluate(pol, .5)[:-1], evaluate(pol, -.5)[1:], c.gamma, c.mach)
        wall_pressure = primitive(evaluate(pol[1:2], -.5), c.gamma)[0, 3]
        f[0] = [0, 0, wall_pressure, 0]
        values = evaluate(pol[1:-1], self.nodes)
        s = np.einsum('nqk,q->nk', source(values, self.theta[:, None]+self.h[:, None]*self.nodes, c.gamma), self.weights)
        visc, fv, sv = self.viscous_rhs(u, pol)
        residual = -(f[1:]-f[:-1])/self.h[:, None]+s+visc
        self.balance_error = np.max(abs(np.sum(residual*self.h[:, None], axis=0)+(f[-1]-fv[-1])-(f[0]-fv[0])-np.sum((s+sv)*self.h[:, None], axis=0)))
        if record:
            self.events['rhs_evaluations'] += 1
        return residual

    def time_step(self, u):
        c = self.config
        q = primitive(u, c.gamma)
        adv = c.cfl*np.min(self.h/(abs(q[:, 2])+np.sqrt(c.gamma*q[:, 3]/q[:, 0])))
        temperature = c.gamma*c.mach**2*q[:, 3]/q[:, 0]
        diffusivity = max(4/3, c.gamma/c.prandtl)*viscosity_ratio(temperature, c)/(c.reynolds*q[:, 0])
        diff = .12*np.min(self.h**2/diffusivity)
        return min(adv, diff)

    def wall_state(self):
        q = super().wall_state()
        q[1:3] = 0
        tw = getattr(self, 'wall_thermal_value', self.config.gamma*self.config.mach**2*q[3]/q[0])
        q[0] = self.config.gamma*self.config.mach**2*q[3]/tw
        return q

    def run(self, verbose=False):
        result = super().run(verbose)
        self.rhs(self.u, record=False)
        wall = self.wall_state()
        result.update(model='fixed-radial-station angular viscous approximation',
                      skin_friction_coefficient=float(2*self.last_transport[3]),
                      wall_heat_flux_into_fluid=float(self.last_transport[4]),
                      wall_temperature_ratio=float(self.config.gamma*self.config.mach**2*wall[3]/wall[0]),
                      first_cell_width_deg=float(np.degrees(self.h[0])),
                      wall_no_slip=True, wall_closure_order=2)
        return result


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--cells', type=int, default=120)
    p.add_argument('--mach', type=float, default=7.95)
    p.add_argument('--outer', type=float, default=22)
    p.add_argument('--flux', choices=FLUXES, default='ausm-up2')
    p.add_argument('--reconstruction', choices=['first', 'minmod', 'primitive-minmod', 'muscl', 'cweno3', 'cweno5', 'cweno3-char', 'cweno5-char'], default='primitive-minmod')
    p.add_argument('--wall', choices=['adiabatic', 'isothermal'], default='adiabatic')
    p.add_argument('--wall-temperature', type=float, default=1)
    p.add_argument('--reynolds', type=float, default=420000)
    p.add_argument('--prandtl', type=float, default=.72)
    p.add_argument('--temperature-inf', type=float, default=775.56/(1+.2*7.95**2))
    p.add_argument('--gradient-order', type=int, choices=[2, 4], default=2)
    p.add_argument('--stretch', type=float, default=3)
    p.add_argument('--steps', type=int, default=80000)
    p.add_argument('--output', default='viscous-results/single')
    a = p.parse_args()
    c = ViscousConfig(cells=a.cells, mach=a.mach, outer_deg=a.outer, flux=a.flux, reconstruction=a.reconstruction,
                      thermal_wall=a.wall, wall_temperature_ratio=a.wall_temperature, reynolds=a.reynolds,
                      prandtl=a.prandtl, temperature_inf_K=a.temperature_inf,
                      gradient_order=a.gradient_order, stretch=a.stretch, max_steps=a.steps)
    r = ViscousConicalSolver(c).run(True)
    save_result(r, a.output)
    print({k:v for k,v in r.items() if k not in ('theta','primitive','conserved_average','history')})
