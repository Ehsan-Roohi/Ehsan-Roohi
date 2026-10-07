"""Conservative conical Euler state: [rho, rho*vr, rho*vtheta, rho*E]."""
import numpy as np


def conserved(q, gamma=1.4):
    q = np.asarray(q)
    rho, vr, vt, p = np.moveaxis(q, -1, 0)
    return np.stack((rho, rho*vr, rho*vt,
                     p/(gamma-1) + .5*rho*(vr*vr+vt*vt)), axis=-1)


def primitive(u, gamma=1.4):
    rho = u[..., 0]
    vr, vt = u[..., 1]/rho, u[..., 2]/rho
    p = (gamma-1)*(u[..., 3]-.5*rho*(vr*vr+vt*vt))
    return np.stack((rho, vr, vt, p), axis=-1)


def admissible(u, gamma=1.4, floor=1e-12):
    return np.all(np.isfinite(u)) and np.all(u[..., 0] > floor) and np.all(primitive(u, gamma)[..., 3] > floor)


def physical_flux(u, gamma=1.4):
    q = primitive(u, gamma)
    vt, p = q[..., 2], q[..., 3]
    f = u*vt[..., None]
    f[..., 2] += p
    f[..., 3] += p*vt
    return f


def source(u, theta, gamma=1.4):
    """Unweighted angular balance law, matching the legacy inviscid sources."""
    rho, vr, vt, p = np.moveaxis(primitive(u, gamma), -1, 0)
    cot = 1/np.tan(theta)
    enthalpy_density = u[..., 3]+p
    return np.stack((-rho*(2*vr+vt*cot),
                     -rho*(2*vr*vr-vt*vt+vr*vt*cot),
                     -rho*(3*vr*vt+vt*vt*cot),
                     -enthalpy_density*(2*vr+vt*cot)), axis=-1)


def freestream(theta, mach, gamma=1.4):
    return conserved(np.stack((np.ones_like(theta), np.cos(theta),
                               -np.sin(theta), np.full_like(theta, 1/(gamma*mach*mach))), axis=-1), gamma)


def freestream_average(theta, h, mach, gamma=1.4):
    """Exact angular cell averages of a Cartesian-uniform upstream state."""
    u = freestream(theta, mach, gamma)
    u[..., 1:3] *= np.sinc(np.asarray(h)/(2*np.pi))[..., None]
    return u
