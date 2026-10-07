"""Chandrashekar (2013) entropy-conservative flux with optional LF dissipation.

The face-level entropy identity is verified separately. Reconstruction, walls,
sources and time integration are not covered by a global entropy-stability claim.
"""
import numpy as np
from state import primitive


def logmean(a, b):
    mean = .5*(a+b)
    z = (b-a)/(a+b)
    small = np.abs(z) < 1e-4
    ratio = np.divide(z, np.arctanh(z), out=np.ones_like(z), where=~small)
    return mean*np.where(small, 1-z*z/3-4*z**4/45, ratio)


def entropy_variables(u, gamma=1.4):
    q = primitive(u, gamma)
    rho, vr, vt, p = np.moveaxis(q, -1, 0)
    beta = rho/(2*p)
    s = np.log(p)-gamma*np.log(rho)
    return np.stack(((gamma-s)/(gamma-1)-beta*(vr*vr+vt*vt), 2*beta*vr, 2*beta*vt, -2*beta), axis=-1)


def entropy_conservative(left, right, gamma=1.4):
    l, r = primitive(left, gamma), primitive(right, gamma)
    mean = .5*(l+r)
    beta_l, beta_r = l[..., 0]/(2*l[..., 3]), r[..., 0]/(2*r[..., 3])
    rho_log, beta_log = logmean(l[..., 0], r[..., 0]), logmean(beta_l, beta_r)
    pressure = mean[..., 0]/(beta_l+beta_r)
    mass = rho_log*mean[..., 2]
    radial = mean[..., 1]*mass
    normal = mean[..., 2]*mass+pressure
    velocity_square_mean = .5*(np.sum(l[..., 1:3]**2, axis=-1)+np.sum(r[..., 1:3]**2, axis=-1))
    energy = (1/(2*(gamma-1)*beta_log)-.5*velocity_square_mean)*mass+mean[..., 1]*radial+mean[..., 2]*normal
    return np.stack((mass, radial, normal, energy), axis=-1)


def ec_lf(left, right, gamma=1.4):
    l, r = primitive(left, gamma), primitive(right, gamma)
    speed = np.maximum(np.abs(l[..., 2])+np.sqrt(gamma*l[..., 3]/l[..., 0]),
                       np.abs(r[..., 2])+np.sqrt(gamma*r[..., 3]/r[..., 0]))
    return entropy_conservative(left, right, gamma)-.5*speed[..., None]*(right-left)
