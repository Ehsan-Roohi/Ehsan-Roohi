"""Local Euler right/left eigenvectors for angular normal velocity.

Characteristic CWENO is an optional reconstruction; flux formulas are unchanged.
"""
import numpy as np
from state import primitive


def basis(u, gamma=1.4):
    q = primitive(u, gamma)
    rho, vr, vt, p = np.moveaxis(q, -1, 0)
    c = np.sqrt(gamma*p/rho)
    h = (u[..., 3]+p)/rho
    one, zero = np.ones_like(rho), np.zeros_like(rho)
    rm = np.stack((one, vr, vt-c, h-vt*c), -1)
    rc = np.stack((one, vr, vt, .5*(vr**2+vt**2)), -1)
    rs = np.stack((zero, one, zero, vr), -1)
    rp = np.stack((one, vr, vt+c, h+vt*c), -1)
    right = np.stack((rm, rc, rs, rp), -1)
    return right, np.linalg.inv(right)
