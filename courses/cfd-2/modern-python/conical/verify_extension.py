"""Independent flux, finite-volume geometry and constitutive verification.

These checks do not establish full axisymmetric cone validation. Results and
observed rates are saved separately from the earlier inviscid verification.
"""
from pathlib import Path
from dataclasses import replace
import json
import numpy as np
from state import conserved, primitive, physical_flux
from flux import FLUXES, numerical_flux
from entropy_flux import entropy_variables, entropy_conservative, ec_lf
from fv_mesh import Mesh
from reconstruction import evaluate
from viscous import ViscousConfig, viscosity_ratio, transport

ROOT = Path(__file__).resolve().parent


def average(fun, mesh):
    x, w = np.polynomial.legendre.leggauss(12)
    theta = mesh.center[:, None]+mesh.width[:, None]*x/2
    return np.einsum('nqk,q->nk', fun(theta), w/2)


def verify():
    rng = np.random.default_rng(725)
    q = np.c_[rng.uniform(.2, 3, 200), rng.uniform(-2, 2, (200, 2)), rng.uniform(.05, 3, 200)]
    u = conserved(q)
    v = conserved(q[rng.permutation(len(q))])
    errors = {}
    flip = np.array([1, 1, -1, 1])
    for name in FLUXES:
        f = numerical_flux(name, u, v)
        errors[name] = dict(equal_state=float(np.max(abs(numerical_flux(name, u, u)-physical_flux(u)))),
                            orientation=float(np.max(abs(numerical_flux(name, v*flip, u*flip)+f*flip))))
    dv = entropy_variables(v)-entropy_variables(u)
    potential = v[:, 2]-u[:, 2]
    identity = np.sum(dv*entropy_conservative(u, v), axis=1)-potential
    inequality = np.sum(dv*ec_lf(u, v), axis=1)-potential
    entropy = {'identity_max_error': float(np.max(abs(identity))),
               'largest_dissipative_face_production': float(np.max(inequality))}
    c = ViscousConfig()
    theta = np.linspace(.2, 1, 100)
    rigid = np.stack((np.ones_like(theta), np.cos(theta), -np.sin(theta), np.ones_like(theta)), -1)
    derivative = np.stack((np.zeros_like(theta), -np.sin(theta), -np.cos(theta), np.zeros_like(theta)), -1)
    _, _, stresses = transport(rigid, derivative, theta, c)
    random_d = rng.uniform(-1, 1, (100, 4))
    _, _, traces = transport(rigid, random_d, theta, c)
    constitutive = {'uniform_cartesian_velocity_stress': float(np.max(abs(stresses))),
                    'stress_trace': float(np.max(abs(traces[:, :3].sum(axis=1)))),
                    'sutherland_normalization': float(abs(viscosity_ratio(np.array([1.]), c)[0]-1))}
    f1,s1,t1 = transport(rigid,random_d,theta,c)
    f2,s2,t2 = transport(rigid,random_d,theta,replace(c,reynolds=2*c.reynolds))
    _,_,tp = transport(rigid,random_d,theta,replace(c,prandtl=2*c.prandtl))
    constitutive['reynolds_scaling'] = float(max(np.max(abs(f2-.5*f1)),np.max(abs(s2-.5*s1))))
    constitutive['prandtl_heat_scaling'] = float(np.max(abs(tp[:,4]-.5*t1[:,4])))
    rates = {}
    for order in (2, 4):
        rows = []
        for n in (24, 48, 96, 192):
            mesh = Mesh(.2, 1, n, 2)
            mean = average(lambda x: np.sin(2*x)[..., None], mesh)
            value, grad = mesh.face_value_gradient(mean, order)
            # Exact cell average of second derivative = endpoint derivative difference / h.
            numerical = (grad[1:, 0]-grad[:-1, 0])/mesh.width[3:-3]
            exact = 2*(np.cos(2*mesh.faces[1:])-np.cos(2*mesh.faces[:-1]))/mesh.width[3:-3]
            error = float(np.sqrt(np.mean((numerical[2:-2]-exact[2:-2])**2)))
            rows.append({'cells': n, 'l2_error': error,
                         'order': None if not rows else float(np.log2(rows[-1]['l2_error']/error))})
        rates[str(order)] = rows
    reconstruction = {}
    for kind in ('cweno3', 'cweno5'):
        rows = []
        for n in (24, 48, 96, 192):
            mesh = Mesh(.2, 1, n, 2)
            mean = average(lambda x: (1+.2*np.sin(2*x))[..., None], mesh)
            pol = mesh.reconstruct(mean, kind)
            x, w = np.polynomial.legendre.leggauss(6)
            recovered = np.einsum('nqk,q->nk', evaluate(pol, x/2), w/2)
            conservation = float(np.max(abs(recovered-mean[2:-2])))
            face = evaluate(pol, .5)[:-1, 0]
            rhs = (face[1:]-face[:-1])/mesh.width[3:-3]
            exact = .2*(np.sin(2*mesh.faces[1:])-np.sin(2*mesh.faces[:-1]))/mesh.width[3:-3]
            error = float(np.sqrt(np.mean((rhs-exact)**2)))
            rows.append({'cells': n, 'l2_error': error, 'average_error': conservation,
                         'order': None if not rows else float(np.log2(rows[-1]['l2_error']/error))})
        reconstruction[kind] = rows
    manufactured = {}
    def profile(x):
        return np.stack((1+.1*np.cos(x), .6+.05*np.cos(2*x), .02*np.sin(3*x), 1+.1*np.sin(x)), -1)
    def profile_derivative(x):
        return np.stack((-.1*np.sin(x), -.1*np.sin(2*x), .06*np.cos(3*x), .1*np.cos(x)), -1)
    for order in (2, 4):
        rows = []
        for n in (24, 48, 96, 192):
            mesh = Mesh(.2, 1, n, 2)
            mean = average(profile, mesh)
            face, derivative = mesh.face_value_gradient(mean, order)
            fv, _, _ = transport(face, derivative, mesh.faces, c)
            qp = mesh.reconstruct(mean, 'cweno3' if order == 2 else 'cweno5')[1:-1]
            x, w = np.polynomial.legendre.leggauss(8)
            vals = evaluate(qp, x/2)
            deriv = evaluate(qp[:, 1:]*np.arange(1, qp.shape[1])[None, :, None], x/2)/mesh.width[3:-3, None, None]
            points = mesh.center[3:-3, None]+mesh.width[3:-3, None]*x/2
            _, sv, _ = transport(vals, deriv, points, c)
            numerical = np.diff(fv, axis=0)/mesh.width[3:-3, None]+np.einsum('nqk,q->nk', sv, w/2)
            exact_f, _, _ = transport(profile(mesh.faces), profile_derivative(mesh.faces), mesh.faces, c)
            _, exact_s, _ = transport(profile(points), profile_derivative(points), points, c)
            exact = np.diff(exact_f, axis=0)/mesh.width[3:-3, None]+np.einsum('nqk,q->nk', exact_s, w/2)
            error = float(np.sqrt(np.mean((numerical[2:-2]-exact[2:-2])**2)))
            rows.append({'cells': n, 'l2_error': error, 'order':None if not rows else float(np.log2(rows[-1]['l2_error']/error))})
        manufactured[str(order)] = rows
    # Constant-property planar Couette flow: viscous heating balanced by heat conduction.
    # Exact primitive cell averages, cubic FV face operator, independent analytic temperature.
    cc = ViscousConfig(viscosity='constant', mach=2, reynolds=1000)
    mesh = Mesh(0, 1, 64, 1)
    amplitude = .5*cc.prandtl*(cc.gamma-1)*cc.mach**2
    qm = average(lambda y: np.stack((np.ones_like(y), y, np.zeros_like(y), 1+amplitude*y*(1-y)), -1), mesh)
    face, derivative = mesh.face_value_gradient(qm, 4)
    tau = derivative[:, 1]/cc.reynolds
    heat = -derivative[:, 3]/(cc.reynolds*cc.prandtl*(cc.gamma-1)*cc.mach**2)
    energy_flux = tau*face[:, 1]-heat
    couette = {'shear_error': float(np.max(abs(tau-1/cc.reynolds))),
               'heat_into_lower_wall': float(heat[0]), 'heat_into_upper_wall': float(-heat[-1]),
               'energy_residual_max': float(np.max(abs(np.diff(energy_flux)/mesh.width[3:-3])))}
    characteristic = {}
    for kind in ('cweno3-char', 'cweno5-char', 'primitive-minmod'):
        mesh = Mesh(.2, 1, 80, 2)
        means = average(lambda x: conserved(np.stack((1+.1*np.sin(x), .6+.02*np.cos(x),
                            .1+.02*np.sin(x), .3+.05*np.cos(x)), -1)), mesh)
        pol = mesh.reconstruct(means, kind)
        x, w = np.polynomial.legendre.leggauss(8)
        recovered = np.einsum('nqk,q->nk', evaluate(pol, x/2), w/2)
        characteristic[kind] = {'average_error':float(np.max(abs(recovered-means[2:-2])))}
    char_rates = {}
    def smooth_euler(x):
        return conserved(np.stack((1+.1*np.sin(2*x),.6+.02*np.cos(2*x),
                                    .1+.02*np.sin(2*x),.3+.05*np.cos(2*x)),-1))
    for kind in ('cweno3-char','cweno5-char'):
        rows = []
        for n in (24,48,96,192):
            mesh = Mesh(.2,1,n,2)
            means = average(smooth_euler,mesh)
            pol = mesh.reconstruct(means,kind)
            flux = numerical_flux('hllc',evaluate(pol,.5)[:-1],evaluate(pol,-.5)[1:])
            numerical = np.diff(flux,axis=0)/mesh.width[3:-3,None]
            exact = np.diff(physical_flux(smooth_euler(mesh.faces)),axis=0)/mesh.width[3:-3,None]
            error = float(np.sqrt(np.mean((numerical-exact)**2)))
            rows.append({'cells':n,'l2_error':error,'order':None if not rows else float(np.log2(rows[-1]['l2_error']/error))})
        char_rates[kind] = rows
    assert char_rates['cweno3-char'][-1]['order'] > 2.7 and char_rates['cweno5-char'][-1]['order'] > 4.3
    assert all(row['average_error'] < 1e-12 for row in characteristic.values())
    # Independently audit colored Jacobian support against one-column differences.
    from viscous import ViscousConicalSolver, thermal_primitive
    from steady import colored_jacobian
    js = ViscousConicalSolver(ViscousConfig(cells=16,reconstruction='primitive-minmod'))
    xx = js.theta
    js.u = conserved(np.stack((1+.1*xx**2,.6+.02*xx**2,-.1-.05*xx**2,.04+.01*xx**2),-1))
    rr = js.rhs(js.u,record=False)
    colored = colored_jacobian(js,js.u,rr)
    brute = np.zeros_like(colored)
    eps = 2e-7*np.maximum(abs(js.u),.01)
    for cell in range(len(js.u)):
        for component in range(4):
            perturbed = js.u.copy()
            perturbed[cell,component] += eps[cell,component]
            brute[:,4*cell+component] = ((js.rhs(perturbed,record=False)-rr)/eps[cell,component]).ravel()
    jacobian_error = float(np.max(abs(colored-brute)))
    assert jacobian_error < 1e-8
    # Exact quadratic cell-average Dirichlet closure, not a nodal derivative.
    vals = 2+3*js.wall_moments[:,0]+.7*js.wall_moments[:,1]
    recovered = js.wall_dirichlet_fit@(vals-2)
    wall_fit_error = float(np.max(abs(recovered-np.array([3,.7]))))
    assert wall_fit_error < 1e-12
    assert all(max(e.values()) < 1e-10 for e in errors.values()), errors
    assert entropy['identity_max_error'] < 1e-9 and entropy['largest_dissipative_face_production'] <= 1e-10
    assert max(constitutive.values()) < 1e-12
    assert rates['2'][-1]['order'] > 1.8 and rates['4'][-1]['order'] > 3.7
    assert manufactured['2'][-1]['order'] > 1.8 and manufactured['4'][-1]['order'] > 3.7
    assert all(rows[-1]['order'] > int(name[-1])-.3 for name, rows in reconstruction.items())
    assert couette['energy_residual_max'] < 1e-9
    result = {'passed': True, 'flux_properties': errors, 'face_entropy': entropy,
              'constitutive': constitutive, 'fv_diffusion_rates': rates,
              'smooth_viscous_angular_operator': manufactured,
              'nonuniform_cweno_rates': reconstruction, 'couette_heating': couette,
              'additional_reconstruction_conservation': characteristic,
              'characteristic_euler_operator_rates':char_rates,
              'colored_vs_brute_jacobian_max_error':jacobian_error,
              'quadratic_fv_wall_fit_error':wall_fit_error,
              'scope': 'Smooth interior diffusion order; wall closure remains second order. No global high-order entropy theorem.'}
    out = ROOT/'results-v2'
    out.mkdir(exist_ok=True)
    (out/'verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    verify()
