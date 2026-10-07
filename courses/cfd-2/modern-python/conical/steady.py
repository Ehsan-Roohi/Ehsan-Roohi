"""Damped pseudo-transient Newton solve of the SAME conservative FV residual.

NumPy linear algebra only. Colored finite-difference Jacobian, positivity-aware
line search, independent residual acceptance, and explicit SSPRK drift check.
The iteration axis is Newton iterations, not SSPRK steps or physical time.
"""
from pathlib import Path
from dataclasses import asdict
from time import perf_counter
import argparse
import json
import numpy as np
from viscous import ViscousConicalSolver, ViscousConfig
from state import admissible, primitive, conserved, freestream
from solver import save_result
from reference import taylor_maccoll, sample_reference

ROOT = Path(__file__).resolve().parent


def colored_jacobian(solver, u, residual, radius=5):
    n = len(u)
    jac = np.zeros((4*n,4*n))
    perturbation = max(5e-10,min(2e-7,1e-3*norm(residual)))
    eps = perturbation*np.maximum(abs(u), .01)
    stride = 2*radius+1
    for component in range(4):
        for color in range(min(stride,n)):
            cells = np.arange(color,n,stride)
            perturbed = u.copy()
            perturbed[cells, component] += eps[cells,component]
            difference = solver.rhs(perturbed,record=False)-residual
            for cell in cells:
                rows = np.arange(max(0,cell-radius),min(n,cell+radius+1))
                jac[(4*rows[:,None]+np.arange(4)).ravel(),4*cell+component] = (difference[rows]/eps[cell,component]).ravel()
    return jac


def norm(r):
    return float(np.max(np.sqrt(np.mean(r*r,axis=0))))


def solve(solver, max_iterations=180, verbose=True):
    c = solver.config
    start = perf_counter()
    history = []
    tau = .01
    merit = np.inf
    jacobian_direction_error = None
    line_reductions = 0
    explicit_recovery_steps = 0
    recovery_count = 0
    for iteration in range(max_iterations+1):
        u = solver.u
        r = solver.rhs(u,record=False)
        residual = norm(r)
        wall = solver.wall_state()
        history.append([iteration,residual,float(wall[3]*c.gamma*c.mach**2)])
        if verbose and iteration%10 == 0:
            print(f'{c.flux}/{c.reconstruction} N={c.cells} Newton={iteration} residual={residual:.3e}',flush=True)
        if residual < c.tolerance or iteration == max_iterations:
            break
        jac = colored_jacobian(solver,u,r)
        if jacobian_direction_error is None:
            direction = np.random.default_rng(71).normal(size=u.shape)*np.maximum(abs(u),.01)
            numerical = (solver.rhs(u+1e-7*direction,record=False)-r)/1e-7
            predicted = (jac@direction.ravel()).reshape(u.shape)
            jacobian_direction_error = float(np.linalg.norm(predicted-numerical)/max(np.linalg.norm(numerical),1e-15))
        merit = np.linalg.norm(r)
        accepted = False
        for retry in range(8):
            try:
                delta = np.linalg.solve(np.eye(u.size)/tau-jac,r.ravel()).reshape(u.shape)
            except np.linalg.LinAlgError:
                tau *= .2
                continue
            alpha = 1.
            for _ in range(20):
                candidate = u+alpha*delta
                if admissible(candidate,c.gamma):
                    try:
                        new_r = solver.rhs(candidate,record=False)
                        if np.linalg.norm(new_r) <= merit*(1-1e-4*alpha):
                            accepted = True
                            solver.u = candidate
                            break
                    except ArithmeticError:
                        pass
                alpha *= .5
                line_reductions += 1
            if accepted:
                tau = min(1e8,tau*(1.6 if alpha >= .5 else .8))
                break
            tau *= .2
        if not accepted:
            if recovery_count >= 6:
                break
            # The nonnormal transient can initially increase the residual norm.
            # A short explicit evolution supplies a safer new Newton iterate.
            for _ in range(200):
                dt = solver.time_step(solver.u)
                for attempt in range(13):
                    try:
                        solver.u = solver.step(solver.u,dt)
                        break
                    except ArithmeticError:
                        dt *= .5
                else:
                    raise RuntimeError('Explicit Newton recovery remained inadmissible')
            explicit_recovery_steps += 200
            recovery_count += 1
            tau = .005
    solver.rhs(solver.u,record=False)
    q = solver.center_primitive()
    wall = solver.wall_state()
    result = {'config':asdict(c),'converged':residual < c.tolerance,'iterations':iteration,
              'residual_rms_max_equation':residual,'elapsed_seconds':perf_counter()-start,
              'steady_solver':'damped pseudo-transient Newton','history_axis':'Newton iterations',
              'max_newton_iterations':max_iterations,
              'termination':'residual tolerance' if residual < c.tolerance else ('iteration budget' if iteration == max_iterations else 'line search stalled'),
              'jacobian_direction_relative_error':jacobian_direction_error,'line_search_reductions':line_reductions,
              'explicit_recovery_steps':explicit_recovery_steps,
              'wall_pressure_ratio':float(wall[3]*c.gamma*c.mach**2),'wall_mach':0.,
              'wall_temperature_ratio':float(c.gamma*c.mach**2*wall[3]/wall[0]),
              'skin_friction_coefficient':float(2*solver.last_transport[3]),
              'wall_heat_flux_into_fluid':float(solver.last_transport[4]),
              'min_density':float(q[:,0].min()),'min_pressure':float(q[:,3].min()),
              'discrete_balance_roundoff':float(solver.balance_error),
              'model':'fixed-radial-station angular viscous approximation',
              'wall_closure_order':2,'first_cell_width_deg':float(np.degrees(solver.h[0])),
              'events':solver.events.copy(),'history':history,'theta':solver.theta,
              'primitive':q,'conserved_average':solver.u.copy()}
    # Independently ask the explicit time integrator whether the solved state drifts.
    if result['converged']:
        original = solver.u.copy()
        dt = solver.time_step(original)
        check = original.copy()
        result['steady_residual_converged'] = True
        try:
            for _ in range(20):
                check = solver.step(check,dt)
            result['explicit_20_step_max_state_drift'] = float(np.max(abs(check-original)))
            result['explicit_20_step_residual'] = norm(solver.rhs(check,record=False))
            result['explicit_drift_check_passed'] = result['explicit_20_step_max_state_drift'] < 1e-7
        except ArithmeticError as error:
            result['explicit_drift_check_passed'] = False
            result['explicit_drift_check_error'] = str(error)
        result['converged'] = result['explicit_drift_check_passed']
        if not result['converged']:
            result['termination'] = 'steady residual reached; independent explicit drift check failed'
    return result


def initialize(solver, checkpoint=None):
    if checkpoint:
        data = np.load(checkpoint)
        if data['conserved_average'].shape != solver.u.shape:
            raise ValueError('Checkpoint grid does not match')
        solver.u = data['conserved_average'].copy()
    else:
        ref = taylor_maccoll(solver.config.mach,max_step=1e-4)
        x,w = np.polynomial.legendre.leggauss(8)
        points = solver.theta[:,None]+solver.h[:,None]*x/2
        solver.u = np.einsum('nqk,q->nk',sample_reference(points.ravel(),ref).reshape(len(solver.u),8,4),w/2)


def remap_checkpoint(solver, checkpoint):
    """Wall-aware cell averaging of a convex conservative interpolant; initializer only."""
    path = Path(checkpoint)
    meta = json.loads(path.with_suffix('.json').read_text(encoding='utf-8'))
    data = np.load(path)
    coarse = ViscousConicalSolver(ViscousConfig(**meta['config']))
    coarse.u = data['conserved_average'].copy()
    coarse.rhs(coarse.u,record=False)
    wall = conserved(coarse.wall_state()[None,:],coarse.config.gamma)
    x = np.r_[coarse.mesh.faces[0],data['theta'],coarse.mesh.faces[-1]]
    values = np.r_[wall,data['conserved_average'],freestream(np.array([coarse.mesh.faces[-1]]),coarse.config.mach,coarse.config.gamma)]
    nodes,weights = np.polynomial.legendre.leggauss(8)
    points = solver.theta[:,None]+solver.h[:,None]*nodes/2
    mapped = np.stack([np.interp(points,x,values[:,k]) for k in range(4)],-1)
    solver.u = np.einsum('nqk,q->nk',mapped,weights/2)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--flux',default='ausm-up2')
    p.add_argument('--cells',type=int,default=60)
    p.add_argument('--reconstruction',default='primitive-minmod')
    p.add_argument('--wall',default='adiabatic',choices=['adiabatic','isothermal'])
    p.add_argument('--gradient-order',type=int,default=2,choices=[2,4])
    p.add_argument('--checkpoint')
    p.add_argument('--iterations',type=int,default=180)
    p.add_argument('--relax-steps',type=int,default=0)
    p.add_argument('--output',default=str(ROOT/'results-v2/newton-pilot'))
    a = p.parse_args()
    solver = ViscousConicalSolver(ViscousConfig(cells=a.cells,flux=a.flux,reconstruction=a.reconstruction,
                                              thermal_wall=a.wall,gradient_order=a.gradient_order))
    initialize(solver,a.checkpoint)
    if a.relax_steps:
        solver.config.max_steps = a.relax_steps
        relaxation = solver.run(True)
        save_result(relaxation,str(a.output)+'-relaxation')
    result = solve(solver,a.iterations)
    save_result(result,a.output)
    print(json.dumps({k:v for k,v in result.items() if k not in ('theta','primitive','conserved_average','history')},indent=2))
