"""Student-facing API shared by the self-contained Colab notebook and local runs."""
from dataclasses import asdict
from pathlib import Path
import json
import numpy as np
from flux import FLUXES
from solver import Config, ConicalSolver, save_result
from viscous import ViscousConfig, ViscousConicalSolver
from dg import DGConfig, DGConicalSolver, save_dg
from reference import taylor_maccoll
from study import metrics

ROOT = Path(__file__).resolve().parent


def run_case(flux='hllc', physics='inviscid', cells=60, mach=2.35,
             outer_deg=45, reconstruction='cweno3', max_steps=12000,
             degree=1, gradient_order=2, reynolds=420000, prandtl=.72,
             wall_temperature_ratio=1., max_newton=180, use_seed=True,
             output_dir=None, verbose=False):
    if flux not in FLUXES and flux!='jst':
        raise ValueError(f'Unknown flux: {flux}')
    if flux=='jst' and physics!='inviscid':
        raise ValueError('JST is implemented only as a standalone uniform-grid inviscid FV scheme')
    target = Path(output_dir or ROOT/'student-results')
    target.mkdir(parents=True, exist_ok=True)
    seed_name = None
    if physics == 'dg-inviscid':
        solver = DGConicalSolver(DGConfig(flux=flux,cells=cells,mach=mach,outer_deg=outer_deg,
                                         degree=degree,max_steps=max_steps))
        result = solver.run(verbose)
    elif physics == 'inviscid':
        if flux=='jst' and reconstruction!='first':
            raise ValueError('JST is a standalone central stencil scheme; select first')
        if reconstruction == 'primitive-minmod':
            raise ValueError('Select first, muscl or CWENO for the uniform inviscid FV solver')
        solver = ConicalSolver(Config(flux=flux,cells=cells,mach=mach,outer_deg=outer_deg,
                                      reconstruction=reconstruction,max_steps=max_steps))
        result = solver.run(verbose)
    elif physics in ('adiabatic','isothermal'):
        if reconstruction not in ('primitive-minmod','first','minmod','cweno3','cweno5','cweno3-char','cweno5-char'):
            raise ValueError('Select primitive-minmod or a supported CWENO viscous reconstruction')
        cfg = ViscousConfig(flux=flux,cells=cells,mach=mach,outer_deg=outer_deg,
                            thermal_wall=physics,reconstruction=reconstruction,
                            max_steps=max_steps,gradient_order=gradient_order,
                            reynolds=reynolds,prandtl=prandtl,wall_temperature_ratio=wall_temperature_ratio)
        solver = ViscousConicalSolver(cfg)
        seed = ROOT/'results-v2'/f'fv2-steady-{physics}-m7.95-n60-ausm-up2-primitive-minmod-g2.npz'
        if use_seed and seed.exists():
            meta = json.loads(seed.with_suffix('.json').read_text(encoding='utf-8'))
            keys = ('mach','cone_deg','outer_deg','cells','gamma','reynolds','prandtl',
                    'temperature_inf_K','viscosity','thermal_wall','wall_temperature_ratio','stretch')
            if all(asdict(cfg)[k] == meta['config'][k] for k in keys):
                solver.u = np.load(seed)['conserved_average'].copy()
                seed_name = seed.name
        if seed_name is None:
            from steady import initialize
            # initialize() integrates the reference over the nonuniform FV cells.
            initialize(solver)
        from steady import solve
        result = solve(solver,max_iterations=max_newton,verbose=verbose)
        result['initialization'] = ('Published, wall/physics/grid-matched coarse steady seed: '+seed_name
                                    if seed_name else 'Taylor-Maccoll cell averages, initial guess only')
    else:
        raise ValueError('Choose inviscid, adiabatic, isothermal or dg-inviscid')
    if physics in ('inviscid','dg-inviscid'):
        reference = taylor_maccoll(mach,cone_deg=10,gamma=1.4,max_step=1e-4)
        result = metrics(result,reference)
    tag = f'{physics}-{flux}-n{cells}-'+(f'p{degree}' if physics == 'dg-inviscid' else reconstruction)
    if physics == 'dg-inviscid':
        save_dg(result,target/tag)
    else:
        save_result(result,target/tag)
    result['file'] = tag
    print(f'{tag}: converged={result["converged"]}, residual={result["residual_rms_max_equation"]:.3e}')
    if not result['converged']:
        print('Exploratory result: the acceptance gate was not met. Do not treat the plotted curve as validated accuracy.')
    if seed_name:
        print('A converged coarse seed was used. Iteration counts are not fair flux-speed rankings.')
    return result


def plot_comparison(results, output_dir=None):
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18})
    target = Path(output_dir or ROOT/'student-results')
    target.mkdir(parents=True,exist_ok=True)
    fig, axes = plt.subplots(2,4,figsize=(16,8.8),layout='constrained')
    colors = plt.cm.tab20(np.linspace(0,1,max(1,len(results))))
    for result,color in zip(results,colors):
        c,q = result['config'],result['primitive']
        theta = np.degrees(result['theta'])
        mach = np.hypot(q[:,1],q[:,2])/np.sqrt(c['gamma']*q[:,3]/q[:,0])
        values = [q[:,3]*c['gamma']*c['mach']**2,q[:,0],mach,
                  c['gamma']*c['mach']**2*q[:,3]/q[:,0],q[:,1],q[:,2]]
        label = c['flux']+(' / DG p='+str(c['degree']) if 'degree' in c else '')
        label += ' [accepted]' if result['converged'] else ' [gate not met]'
        style = '-' if result['converged'] else '--'
        for ax,value in zip(axes.flat,values):
            ax.plot(theta,value,style,color=color,lw=1.5,label=label)
        hist = np.array(result['history'])
        axes[1,2].semilogy(hist[:,0],np.maximum(hist[:,1],1e-18),color=color,lw=1.4,label=label)
    labels = [r'$p/p_\infty$',r'$\rho/\rho_\infty$','Local Mach number',r'$T/T_\infty$',
              r'$v_r/U_\infty$',r'$v_\theta/U_\infty$']
    for ax,label in zip(axes.flat,labels):
        ax.set(xlabel=r'Polar angle $\theta$ [deg]',ylabel=label)
    axes[1,2].axhline(2e-8,color='black',ls=':',lw=1)
    axes[1,2].set(xlabel='SSPRK steps or Newton attempts (see run metadata)',ylabel='Maximum RMS residual')
    axes[1,3].axis('off')
    handles,labels = axes[0,0].get_legend_handles_labels()
    axes[1,3].legend(handles,labels,loc='center left',frameon=False,fontsize=8)
    fig.suptitle('Conical flow: computed properties and convergence',fontsize=20)
    fig.savefig(target/'flow-comparison.png',dpi=190)
    fig.savefig(target/'flow-comparison.pdf')
    plt.show()
    return fig
