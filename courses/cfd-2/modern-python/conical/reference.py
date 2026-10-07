"""Independent Taylor-Maccoll shooting reference; RK4 and bisection, no CFD package."""
from math import sin, cos, tan, sqrt, asin, radians, degrees, ceil
import numpy as np


def shock_state(beta, mach, gamma):
    mn2 = (mach*sin(beta))**2
    density = (gamma+1)*mn2/(2+(gamma-1)*mn2)
    p = (1+2*gamma/(gamma+1)*(mn2-1))/(gamma*mach*mach)
    return cos(beta), -sin(beta)/density, density, p


def integrate(beta, cone, mach, gamma, max_step=2e-4):
    vr, vt, _, _ = shock_state(beta, mach, gamma)
    n = int(ceil((beta-cone)/max_step))
    step = (cone-beta)/n
    theta = beta
    theta_values, velocity = [theta], [(vr, vt)]
    total = 1+2/((gamma-1)*mach*mach)
    def rhs(t, a, b):
        c2 = .5*(gamma-1)*(total-a*a-b*b)
        den = c2-b*b
        if c2 <= 0 or abs(den) < 1e-10:
            raise ArithmeticError("Taylor-Maccoll sonic/singular branch")
        return b, (b*b*a-c2*(2*a+b/tan(t)))/den
    for _ in range(n):
        k1 = rhs(theta, vr, vt)
        k2 = rhs(theta+step/2, vr+step*k1[0]/2, vt+step*k1[1]/2)
        k3 = rhs(theta+step/2, vr+step*k2[0]/2, vt+step*k2[1]/2)
        k4 = rhs(theta+step, vr+step*k3[0], vt+step*k3[1])
        vr += step*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        vt += step*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        theta += step
        theta_values.append(theta)
        velocity.append((vr, vt))
    return np.array(theta_values), np.array(velocity)


def taylor_maccoll(mach, cone_deg=10, gamma=1.4, max_step=2e-4):
    cone = radians(cone_deg)
    begin = max(asin(1/mach), cone)+1e-5
    last = None
    bracket = None
    for beta in np.linspace(begin, radians(75), 140):
        try:
            end = integrate(beta, cone, mach, gamma, max_step)[1][-1, 1]
        except (ArithmeticError, OverflowError):
            last = None
            continue
        if last is not None and last[1]*end < 0:
            bracket = (last[0], beta, last[1])
            break
        last = (beta, end)
    if bracket is None:
        raise ValueError("No attached weak-shock bracket found")
    lo, hi, flo = bracket
    for _ in range(40):
        mid = .5*(lo+hi)
        fm = integrate(mid, cone, mach, gamma, max_step)[1][-1, 1]
        if fm*flo > 0:
            lo, flo = mid, fm
        else:
            hi = mid
    beta = .5*(lo+hi)
    theta, v = integrate(beta, cone, mach, gamma, max_step)
    vr2, vt2, rho2, p2 = shock_state(beta, mach, gamma)
    total = 1+2/((gamma-1)*mach*mach)
    a2 = .5*(gamma-1)*(total-np.sum(v*v, axis=1))
    a2shock = .5*(gamma-1)*(total-vr2*vr2-vt2*vt2)
    ratio = a2/a2shock
    rho = rho2*ratio**(1/(gamma-1))
    p = p2*ratio**(gamma/(gamma-1))
    q = np.column_stack((rho, v, p))
    return {"mach": mach, "gamma": gamma, "cone_deg": cone_deg,
            "shock_angle_deg": degrees(beta), "theta": theta[::-1], "primitive": q[::-1],
            "wall_pressure_ratio": p[-1]*gamma*mach*mach,
            "wall_mach": sqrt(np.sum(v[-1]*v[-1])/a2[-1]), "wall_density_ratio": rho[-1]}


def sample_reference(theta, reference):
    from state import freestream, conserved
    gamma, mach = reference["gamma"], reference["mach"]
    u = freestream(theta, mach, gamma)
    inside = theta <= np.radians(reference["shock_angle_deg"])
    q = np.column_stack([np.interp(theta[inside], reference["theta"], reference["primitive"][:, k]) for k in range(4)])
    u[inside] = conserved(q, gamma)
    return u


if __name__ == "__main__":
    for m in [2.35, 7.95]:
        r = taylor_maccoll(m)
        print({k: v for k, v in r.items() if k not in ["theta", "primitive"]})
