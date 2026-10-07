"""Independent vectorized implementations of classical and robust Euler fluxes."""
import numpy as np
from state import primitive, physical_flux


def waves(u, gamma):
    q = primitive(u, gamma)
    return q, np.sqrt(gamma*q[..., 3]/q[..., 0])


def hlle(left, right, gamma=1.4):
    ql, cl = waves(left, gamma)
    qr, cr = waves(right, gamma)
    sl = np.minimum(np.minimum(ql[..., 2]-cl, qr[..., 2]-cr), 0)
    sr = np.maximum(np.maximum(ql[..., 2]+cl, qr[..., 2]+cr), 0)
    fl, fr = physical_flux(left, gamma), physical_flux(right, gamma)
    return (sr[..., None]*fl-sl[..., None]*fr+
            (sl*sr)[..., None]*(right-left))/(sr-sl)[..., None]


def rusanov(left, right, gamma=1.4):
    ql, cl = waves(left, gamma)
    qr, cr = waves(right, gamma)
    a = np.maximum(abs(ql[..., 2])+cl, abs(qr[..., 2])+cr)
    return .5*(physical_flux(left, gamma)+physical_flux(right, gamma)-a[..., None]*(right-left))


def roe(left, right, gamma=1.4, entropy_fix=True):
    ql, cl = waves(left, gamma)
    qr, cr = waves(right, gamma)
    rl, rr = np.sqrt(ql[..., 0]), np.sqrt(qr[..., 0])
    denom = rl+rr
    vr = (rl*ql[..., 1]+rr*qr[..., 1])/denom
    vt = (rl*ql[..., 2]+rr*qr[..., 2])/denom
    hl = (left[..., 3]+ql[..., 3])/ql[..., 0]
    hr = (right[..., 3]+qr[..., 3])/qr[..., 0]
    h = (rl*hl+rr*hr)/denom
    c2 = (gamma-1)*(h-.5*(vr*vr+vt*vt))
    c = np.sqrt(np.maximum(c2, 1e-14))
    r = rl*rr
    dr, dvr, dvt, dp = np.moveaxis(qr-ql, -1, 0)
    am = (dp-r*c*dvt)/(2*c2)
    ap = (dp+r*c*dvt)/(2*c2)
    ac = dr-dp/c2
    ash = r*dvr
    ones, zeros = np.ones_like(vr), np.zeros_like(vr)
    rm = np.stack((ones, vr, vt-c, h-vt*c), axis=-1)
    rp = np.stack((ones, vr, vt+c, h+vt*c), axis=-1)
    rc = np.stack((ones, vr, vt, .5*(vr*vr+vt*vt)), axis=-1)
    rsh = np.stack((zeros, ones, zeros, vr), axis=-1)
    # Harten-style smoothing of acoustic eigenvalues; not an entropy-stability proof.
    delta = .1*c
    def acoustic(lam):
        a = abs(lam)
        return np.where(a < delta, .5*(a*a/delta+delta), a) if entropy_fix else a
    diss = ((acoustic(vt-c)*am)[..., None]*rm+
            (acoustic(vt+c)*ap)[..., None]*rp+
            (abs(vt)*ac)[..., None]*rc+(abs(vt)*ash)[..., None]*rsh)
    return .5*(physical_flux(left, gamma)+physical_flux(right, gamma)-diss)


def vanleer_split(u, sign, gamma=1.4):
    q, c = waves(u, gamma)
    rho, vr, vt, p = np.moveaxis(q, -1, 0)
    m = vt/c
    mass = sign*rho*c*.25*(m+sign)**2
    normal = ((gamma-1)*vt+sign*2*c)/gamma
    energy = ((gamma-1)*vt+sign*2*c)**2/(2*(gamma*gamma-1))+.5*vr*vr
    sub = np.stack((mass, mass*vr, mass*normal, mass*energy), axis=-1)
    full = physical_flux(u, gamma)
    return np.where((sign*m >= 1)[..., None], full,
                    np.where((sign*m <= -1)[..., None], 0, sub))


def vanleer(left, right, gamma=1.4):
    return vanleer_split(left, 1, gamma)+vanleer_split(right, -1, gamma)


def hllc(left, right, gamma=1.4):
    ql, cl = waves(left, gamma)
    qr, cr = waves(right, gamma)
    rl, ul, vl, pl = np.moveaxis(ql, -1, 0)
    rr, ur, vr, pr = np.moveaxis(qr, -1, 0)
    sl = np.minimum(vl-cl, vr-cr)
    sr = np.maximum(vl+cl, vr+cr)
    sm = (pr-pl+rl*vl*(sl-vl)-rr*vr*(sr-vr))/(rl*(sl-vl)-rr*(sr-vr))
    def star(u, rho, tang, vel, p, s):
        rstar = rho*(s-vel)/(s-sm)
        e = u[..., 3]/rho+(sm-vel)*(sm+p/(rho*(s-vel)))
        return np.stack((rstar, rstar*tang, rstar*sm, rstar*e), axis=-1)
    ls = star(left, rl, ul, vl, pl, sl)
    rs = star(right, rr, ur, vr, pr, sr)
    fl, fr = physical_flux(left, gamma), physical_flux(right, gamma)
    f = np.where((sl >= 0)[..., None], fl,
        np.where((sm >= 0)[..., None], fl+sl[..., None]*(ls-left),
        np.where((sr > 0)[..., None], fr+sr[..., None]*(rs-right), fr)))
    # Face-local fallback is conservative: one replacement flux for both cells.
    valid = ((ls[..., 0] > 0) & (rs[..., 0] > 0) &
             (primitive(ls, gamma)[..., 3] > 0) & (primitive(rs, gamma)[..., 3] > 0) & np.all(np.isfinite(f), axis=-1))
    return np.where(valid[..., None], f, hlle(left, right, gamma))


FLUXES = {"roe": roe, "vanleer": vanleer, "hllc": hllc, "hlle": hlle, "rusanov": rusanov}

from ausm import ausm, ausm_plus, ausm_up, ausm_up2, slau2
from entropy_flux import ec_lf
FLUXES.update({"ausm":ausm, "ausm-plus":ausm_plus, "ausm-up":ausm_up,
               "ausm-up2":ausm_up2, "slau2":slau2, "ec-lf":ec_lf})
from additional_fluxes import godunov, roe_uncorrected, steger_warming, global_lf
FLUXES.update({'godunov':godunov, 'roe-nc':roe_uncorrected,
               'steger-warming':steger_warming, 'global-lf':global_lf})


def numerical_flux(name, left, right, gamma=1.4, mach_ref=.1):
    if name in ["ausm-up", "ausm-up2"]:
        return FLUXES[name](left, right, gamma, mach_ref=mach_ref)
    return FLUXES[name](left, right, gamma)
