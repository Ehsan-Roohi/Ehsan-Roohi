"""Independent AUSM-family mathematical implementations for angular Euler faces.

References: Liou & Steffen (1993), Liou (1996, 2006), Kitamura & Shima
(2013). These are flux algorithms, not turbulence or viscosity models.
The total-enthalpy critical sound speed is used by the AUSM+ variants.
SLAU2 uses arithmetic sound speed; no optional low-dissipation sensor.
"""
import numpy as np
from state import primitive


def mach_split(m, sign, beta=0):
    sonic = .5*(m+sign*np.abs(m))
    sub = sign*(.25*(m+sign)**2+beta*(m*m-1)**2)
    return np.where(np.abs(m) >= 1, sonic, sub)


def pressure_split(m, sign, alpha=0):
    sonic = .5*(1+sign*np.sign(m))
    sub = .25*(m+sign)**2*(2-sign*m)+sign*alpha*m*(m*m-1)**2
    return np.where(np.abs(m) >= 1, sonic, sub)


def ingredients(left, right, gamma):
    l, r = primitive(left, gamma), primitive(right, gamma)
    hl = (left[..., 3]+l[..., 3])/l[..., 0]
    hr = (right[..., 3]+r[..., 3])/r[..., 0]
    cl = np.sqrt(gamma*l[..., 3]/l[..., 0])
    cr = np.sqrt(gamma*r[..., 3]/r[..., 0])
    return l, r, hl, hr, cl, cr


def assemble(mass, pressure, l, r, hl, hr):
    positive = mass >= 0
    transport = np.stack((np.ones_like(mass), np.where(positive, l[..., 1], r[..., 1]),
                          np.where(positive, l[..., 2], r[..., 2]), np.where(positive, hl, hr)), axis=-1)
    flux = mass[..., None]*transport
    flux[..., 2] += pressure
    return flux


def critical_sound(l, r, hl, hr, gamma):
    al = np.sqrt(2*(gamma-1)/(gamma+1)*hl)
    ar = np.sqrt(2*(gamma-1)/(gamma+1)*hr)
    return np.minimum(al*al/np.maximum(al, l[..., 2]), ar*ar/np.maximum(ar, -r[..., 2]))


def ausm(left, right, gamma=1.4):
    """Original local-sound-speed AUSM with the cubic pressure split."""
    l, r, hl, hr, cl, cr = ingredients(left, right, gamma)
    ml, mr = l[..., 2]/cl, r[..., 2]/cr
    mf = mach_split(ml, 1)+mach_split(mr, -1)
    mass = mf*np.where(mf >= 0, l[..., 0]*cl, r[..., 0]*cr)
    pf = pressure_split(ml, 1)*l[..., 3]+pressure_split(mr, -1)*r[..., 3]
    return assemble(mass, pf, l, r, hl, hr)


def plus_family(left, right, gamma, variant, mach_ref):
    l, r, hl, hr, _, _ = ingredients(left, right, gamma)
    a = critical_sound(l, r, hl, hr, gamma)
    ml, mr = l[..., 2]/a, r[..., 2]/a
    mean_m2 = .5*(ml*ml+mr*mr)
    rho = .5*(l[..., 0]+r[..., 0])
    fa = np.ones_like(a)
    alpha = 3/16
    if variant != "plus":
        if mach_ref <= 0:
            raise ValueError("AUSM+-up requires a positive reference Mach number")
        m0 = np.sqrt(np.minimum(1, np.maximum(mean_m2, mach_ref**2)))
        fa = m0*(2-m0)
        alpha = (3/16)*(5*fa*fa-4)
    pl, pr = pressure_split(ml, 1, alpha), pressure_split(mr, -1, alpha)
    mf = mach_split(ml, 1, 1/8)+mach_split(mr, -1, 1/8)
    if variant != "plus":
        mf -= (.25/fa)*np.maximum(1-mean_m2, 0)*(r[..., 3]-l[..., 3])/(rho*a*a)
    mass = a*mf*np.where(mf >= 0, l[..., 0], r[..., 0])
    pf = pl*l[..., 3]+pr*r[..., 3]
    if variant == "up":
        pf -= .75*fa*pl*pr*(2*rho)*a*(r[..., 2]-l[..., 2])
    elif variant == "up2":
        speed = np.sqrt(.5*(np.sum(l[..., 1:3]**2, axis=-1)+np.sum(r[..., 1:3]**2, axis=-1)))
        pf = .5*(l[..., 3]+r[..., 3])+.5*(pl-pr)*(l[..., 3]-r[..., 3])+speed*(pl+pr-1)*rho*a
    return assemble(mass, pf, l, r, hl, hr)


def ausm_plus(left, right, gamma=1.4):
    return plus_family(left, right, gamma, "plus", 1)


def ausm_up(left, right, gamma=1.4, mach_ref=.1):
    return plus_family(left, right, gamma, "up", mach_ref)


def ausm_up2(left, right, gamma=1.4, mach_ref=.1):
    return plus_family(left, right, gamma, "up2", mach_ref)


def slau2(left, right, gamma=1.4):
    l, r, hl, hr, cl, cr = ingredients(left, right, gamma)
    a = .5*(cl+cr)
    ml, mr = l[..., 2]/a, r[..., 2]/a
    speed = np.sqrt(.5*(np.sum(l[..., 1:3]**2, axis=-1)+np.sum(r[..., 1:3]**2, axis=-1)))
    chi = (1-np.minimum(1, speed/a))**2
    g = -np.clip(ml, -1, 0)*np.clip(mr, 0, 1)
    density_velocity = (l[..., 0]*np.abs(l[..., 2])+r[..., 0]*np.abs(r[..., 2]))/(l[..., 0]+r[..., 0])
    abs_l = (1-g)*density_velocity+g*np.abs(l[..., 2])
    abs_r = (1-g)*density_velocity+g*np.abs(r[..., 2])
    mass = .5*(l[..., 0]*(l[..., 2]+abs_l)+r[..., 0]*(r[..., 2]-abs_r)-chi*(r[..., 3]-l[..., 3])/a)
    pl, pr = pressure_split(ml, 1), pressure_split(mr, -1)
    pf = .5*(l[..., 3]+r[..., 3])+.5*(pl-pr)*(l[..., 3]-r[..., 3])+speed*(pl+pr-1)*.5*(l[..., 0]+r[..., 0])*a
    return assemble(mass, pf, l, r, hl, hr)
