"""Pseudo-time FV solver for the archived inviscid angular conical equations."""
from dataclasses import dataclass, asdict
from time import perf_counter
import argparse
import json
from pathlib import Path
import numpy as np
from state import primitive, admissible, source, freestream_average
from flux import FLUXES, numerical_flux
from reconstruction import coefficients, evaluate, limit_admissibility


@dataclass
class Config:
    mach: float = 7.95
    cone_deg: float = 10
    outer_deg: float = 22
    cells: int = 240
    gamma: float = 1.4
    flux: str = "hllc"
    reconstruction: str = "cweno3"
    cfl: float = .35
    tolerance: float = 2e-8
    max_steps: int = 40000
    boundary: str = "wall"
    balance_freestream: bool = False


class ConicalSolver:
    def __init__(self, config):
        self.config = config
        c = config
        self.h = np.radians(c.outer_deg-c.cone_deg)/c.cells
        self.theta = np.radians(c.cone_deg)+(np.arange(c.cells)+.5)*self.h
        self.ext_theta = np.radians(c.cone_deg)+(np.arange(c.cells+6)-2.5)*self.h
        self.upstream_ext = freestream_average(self.ext_theta, self.h, c.mach, c.gamma)
        self.u = self.upstream_ext[3:-3].copy()
        self.nodes, self.weights = np.polynomial.legendre.leggauss(3)
        self.nodes *= .5
        self.weights *= .5
        self.samples = np.r_[-.5, self.nodes, .5]
        self.events = {"reconstruction_limited_cells": 0, "rejected_steps": 0, "rhs_evaluations": 0}
        self.correction = np.zeros_like(self.u)
        if c.balance_freestream:
            # Discrete equilibrium subtraction is optional and explicitly recorded.
            self.correction = self.rhs(self.u, uniform_boundary=True, record=False)

    def ghost(self, u, uniform_boundary=False):
        ext = self.upstream_ext.copy()
        ext[3:-3] = u
        if self.config.boundary == "wall" and not uniform_boundary:
            ext[:3] = u[:3][::-1]
            ext[:3, 2] *= -1
        return ext

    def polynomials(self, u, uniform_boundary=False, record=True):
        ext = self.ghost(u, uniform_boundary)
        c = coefficients(ext, self.config.reconstruction, self.h, self.config.gamma)
        c, limited = limit_admissibility(c, ext[2:-2], self.config.gamma, self.samples)
        if record:
            self.events["reconstruction_limited_cells"] += limited
        return c

    def rhs(self, u, uniform_boundary=False, record=True):
        c = self.config
        pol = self.polynomials(u, uniform_boundary, record)
        left = evaluate(pol, .5)[:-1]
        right = evaluate(pol, -.5)[1:]
        if c.flux=='jst':
            from additional_fluxes import jst_flux
            f=jst_flux(self.ghost(u,uniform_boundary),c.gamma)
        else:
            f = numerical_flux(c.flux, left, right, c.gamma, mach_ref=c.mach)
        if c.boundary == "wall" and not uniform_boundary:
            pwall = primitive(evaluate(pol[1:2], -.5), c.gamma)[0, 3]
            f[0] = [0, 0, pwall, 0]
        values = evaluate(pol[1:-1], self.nodes)
        theta = self.theta[:, None]+self.h*self.nodes[None, :]
        s = np.einsum("nqk,q->nk", source(values, theta, c.gamma), self.weights)
        residual = -(f[1:]-f[:-1])/self.h+s-self.correction
        if record:
            self.events["rhs_evaluations"] += 1
        self.balance_error = np.max(abs(np.sum(residual*self.h, axis=0)+f[-1]-f[0]-np.sum((s-self.correction)*self.h, axis=0)))
        return residual

    def center_primitive(self):
        return primitive(evaluate(self.polynomials(self.u, record=False)[1:-1], 0), self.config.gamma)

    def wall_state(self):
        return primitive(evaluate(self.polynomials(self.u, record=False)[1:2], -.5), self.config.gamma)[0]

    def step(self, u, dt):
        c = self.config
        a = u+dt*self.rhs(u)
        if not admissible(a, c.gamma):
            raise ArithmeticError("inadmissible RK stage 1")
        b = .75*u+.25*(a+dt*self.rhs(a))
        if not admissible(b, c.gamma):
            raise ArithmeticError("inadmissible RK stage 2")
        result = u/3+2*(b+dt*self.rhs(b))/3
        if not admissible(result, c.gamma):
            raise ArithmeticError("inadmissible RK stage 3")
        return result

    def time_step(self, u):
        c = self.config
        q = primitive(u, c.gamma)
        return c.cfl*self.h/np.max(abs(q[:, 2])+np.sqrt(c.gamma*q[:, 3]/q[:, 0]))

    def run(self, verbose=False):
        start = perf_counter()
        c = self.config
        history = []
        converged = False
        residual = float("inf")
        for iteration in range(1, c.max_steps+1):
            dt = self.time_step(self.u)
            old = self.u
            for attempt in range(13):
                try:
                    self.u = self.step(old, dt)
                    break
                except ArithmeticError:
                    self.events["rejected_steps"] += 1
                    dt *= .5
            else:
                raise RuntimeError("Step remained inadmissible after 12 reductions")
            if iteration % 50 == 0 or iteration == 1 or iteration == c.max_steps:
                residual = float(np.max(np.sqrt(np.mean(self.rhs(self.u, record=False)**2, axis=0))))
                wall = self.wall_state()
                history.append([iteration, residual, wall[3]*c.gamma*c.mach*c.mach])
                if verbose and iteration % 1000 == 0:
                    print(f"{c.flux}/{c.reconstruction} N={c.cells} step={iteration} residual={residual:.3e} pwall/pinf={history[-1][2]:.8f}", flush=True)
                if residual < c.tolerance:
                    converged = True
                    break
        elapsed = perf_counter()-start
        q = self.center_primitive()
        wall = self.wall_state()
        machwall = abs(wall[1])/np.sqrt(c.gamma*wall[3]/wall[0])
        return {"config": asdict(c), "converged": converged, "iterations": iteration,
                "residual_rms_max_equation": residual, "elapsed_seconds": elapsed,
                "wall_pressure_ratio": float(wall[3]*c.gamma*c.mach*c.mach), "wall_mach": float(machwall),
                "min_density": float(np.min(q[:, 0])), "min_pressure": float(np.min(q[:, 3])),
                "discrete_balance_roundoff": float(self.balance_error), "events": self.events.copy(),
                "history": history, "theta": self.theta, "primitive": q, "conserved_average": self.u.copy()}


def save_result(result, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Append, rather than replace a suffix: a basename may contain a decimal Mach.
    np.savez_compressed(Path(str(path)+".npz"), theta=result["theta"], primitive=result["primitive"],
                        conserved_average=result["conserved_average"], history=np.array(result["history"]))
    summary = {k: v for k, v in result.items() if k not in ["theta", "primitive", "conserved_average", "history"]}
    Path(str(path)+".json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mach", type=float, default=7.95)
    parser.add_argument("--cells", type=int, default=240)
    parser.add_argument("--outer", type=float, default=22)
    parser.add_argument("--flux", choices=[*FLUXES,'jst'], default="hllc")
    from advanced_reconstruction import RECONSTRUCTIONS
    parser.add_argument("--reconstruction", choices=RECONSTRUCTIONS, default="cweno3")
    parser.add_argument("--steps", type=int, default=40000)
    parser.add_argument("--output", default="results/single")
    args = parser.parse_args()
    config = Config(mach=args.mach, cells=args.cells, outer_deg=args.outer, flux=args.flux,
                    reconstruction=args.reconstruction, max_steps=args.steps)
    result = ConicalSolver(config).run(verbose=True)
    save_result(result, args.output)
    print(json.dumps({k: v for k, v in result.items() if k not in ["theta", "primitive", "conserved_average", "history"]}, indent=2))
