"""English scientific figures and offline report from actual recorded fields."""
from pathlib import Path
import html
import json
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reference import taylor_maccoll, sample_reference
from state import primitive

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results-v2'
METHODS = ['roe','hllc','ausm','ausm-plus','ausm-up','ausm-up2','slau2','ec-lf']
NAMES = dict(zip(METHODS,['Roe','HLLC','AUSM','AUSM+','AUSM+-up','AUSM+-up2','SLAU2','EC + LF']))
COLORS = dict(zip(METHODS,['#475569','#1664a0','#aa641d','#8d5999','#ba4058','#007e83','#55932c','#ef8700']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
    'axes.spines.right':False,'axes.grid':True,'grid.alpha':.15,'figure.facecolor':'white',
    'axes.labelcolor':'#263746','text.color':'#263746','axes.titleweight':'bold',
    'savefig.dpi':210,'pdf.fonttype':42})


def read_rows(pattern):
    rows = []
    for p in sorted(OUT.glob(pattern)):
        row = json.loads(p.read_text(encoding='utf-8'))
        if 'config' in row and 'residual_rms_max_equation' in row and p.with_suffix('.npz').exists():
            row['file'] = p.stem
            rows.append(row)
    return rows


def field(row):
    return np.load(OUT/(row['file']+'.npz'))


def properties(q,mach):
    p = q[:,3]*1.4*mach**2
    return [p,q[:,0],np.linalg.norm(q[:,1:3],axis=1)/np.sqrt(1.4*q[:,3]/q[:,0]),
            p/q[:,0],q[:,1],q[:,2]]


LABELS = [r'$p/p_\infty$',r'$\rho/\rho_\infty$',r'$M$',r'$T/T_\infty$',r'$v_r/U_\infty$',r'$v_\theta/U_\infty$']


def save(fig,name,title,subtitle):
    fig.suptitle(title,x=.075,ha='left',fontsize=18,fontweight='bold',y=.995)
    tall = fig.get_size_inches()[1] >= 7
    fig.text(.075,.925 if tall else .89,subtitle,fontsize=10,color='#5a6c7a')
    fig.tight_layout(rect=[.02,.16,.98,.87 if tall else .82],h_pad=2,w_pad=2)
    fig.savefig(OUT/(name+'.png'))
    fig.savefig(OUT/(name+'.pdf'))
    plt.close(fig)


def rowkey(row):
    c = row['config']
    return (c['thermal_wall'],c['cells'],c['flux'],c['reconstruction'],c['gradient_order'])


def main():
    euler = read_rows('euler-m*-cweno3.json')
    char = read_rows('euler-m*-cweno3-char.json')
    visc = read_rows('fv2-steady-*.json')
    continued = read_rows('continued-fv2-steady-*.json')+read_rows('mesh-continued-fv2-steady-*.json')
    effective = {rowkey(r):r for r in visc}
    for r in continued:
        effective[rowkey(r)] = r
    visc_effective = list(effective.values())
    verification = json.loads((OUT/'verification.json').read_text(encoding='utf-8'))
    for mach,n in [(2.35,60),(7.95,120)]:
        fig,axes = plt.subplots(2,3,figsize=(13,8))
        ref = taylor_maccoll(mach,max_step=1e-4)
        end = 45 if mach < 3 else 22
        theta = np.radians(np.linspace(10,end,1200))
        qref = primitive(sample_reference(theta,ref))
        for ax,value,label in zip(axes.ravel(),properties(qref,mach),LABELS):
            ax.plot(np.degrees(theta),value,'k--',lw=1.6,label='Taylor–Maccoll')
            ax.set_ylabel(label)
            ax.set_xlabel(r'$\theta$ (deg)')
        selected = [r for r in euler if r['config']['mach']==mach and r['config']['cells']==n]
        for r in sorted(selected,key=lambda a:METHODS.index(a['config']['flux'])):
            data = field(r)
            method = r['config']['flux']
            for ax,value in zip(axes.ravel(),properties(data['primitive'],mach)):
                ax.plot(np.degrees(data['theta']),value,color=COLORS[method],lw=1.25,
                        ls='-' if r['converged'] else ':',alpha=1 if r['converged'] else .6,
                        label=NAMES[method]+(' †' if not r['converged'] else ''))
        handles,labels = axes[0,0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='lower center',ncol=5,frameon=False,fontsize=9)
        save(fig,f'euler-properties-m{mach:g}',f'Inviscid cone · Mach {mach:g}',
             f'{n} uniform FV cells · common CWENO3 · dotted † curves exhausted their convergence budget')
    for wall in ['adiabatic','isothermal']:
        selected = [r for r in visc_effective if r['config']['cells']==60 and r['config']['thermal_wall']==wall
                    and r['config']['reconstruction']=='primitive-minmod']
        fig,axes = plt.subplots(2,3,figsize=(13,8))
        for r in sorted(selected,key=lambda a:METHODS.index(a['config']['flux'])):
            data = field(r)
            method = r['config']['flux']
            for ax,value,label in zip(axes.ravel(),properties(data['primitive'],7.95),LABELS):
                ax.plot(np.degrees(data['theta']),value,color=COLORS[method],lw=1.3,
                        ls='-' if r['converged'] else ':',label=NAMES[method]+(' †' if not r['converged'] else ''))
                ax.set_ylabel(label)
                ax.set_xlabel(r'$\theta$ (deg)')
        fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False)
        save(fig,f'viscous-properties-{wall}',f'Viscous cone station · {wall} wall',
             'Mach 7.95 · Re = 420,000 · 60 wall-clustered cells · primitive minmod · local angular model')
    fig,axes = plt.subplots(2,2,figsize=(12,8))
    for column,wall in enumerate(['adiabatic','isothermal']):
        for r in visc_effective:
            c = r['config']
            if c['thermal_wall']!=wall or c['cells']!=60 or c['reconstruction']!='primitive-minmod':
                continue
            data = field(r)
            q = data['primitive']
            x = np.r_[0,np.degrees(data['theta'])-10]
            style = '-' if r['converged'] else ':'
            axes[0,column].plot(x,np.r_[0,q[:,1]],color=COLORS[c['flux']],ls=style,label=NAMES[c['flux']])
            axes[1,column].plot(x,np.r_[r['wall_temperature_ratio'],properties(q,7.95)[3]],color=COLORS[c['flux']],ls=style)
        axes[0,column].set_title(wall.capitalize())
        for ax,label in zip(axes[:,column],[LABELS[4],LABELS[3]]):
            ax.set_xlim(0,1)
            ax.set_xlabel(r'Angular distance from wall $\theta-10^\circ$ (deg)')
            ax.set_ylabel(label)
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False)
    save(fig,'near-wall-profiles','No-slip velocity and wall-temperature effects',
         'True wall values included · cold wall T/T∞ = 1 · angular distance at a fixed radial station')
    fig,axes = plt.subplots(1,3,figsize=(13,5.5))
    for wall,marker in [('adiabatic','o'),('isothermal','s')]:
        for method in ['hllc','ausm-up2']:
            for recon,style in [('primitive-minmod','-'),('cweno5-char','--')]:
                rows = sorted([r for r in visc_effective if r['config']['thermal_wall']==wall and
                              r['config']['flux']==method and r['config']['reconstruction']==recon],key=lambda r:r['config']['cells'])
                for ax,key in zip(axes,['wall_pressure_ratio','skin_friction_coefficient','wall_heat_flux_into_fluid']):
                    x = [r['config']['cells'] for r in rows]
                    y = [r[key] if r['converged'] else np.nan for r in rows]
                    ax.plot(x,y,color=COLORS[method],ls=style,alpha=.7,
                            label=f'{NAMES[method]}, {wall}, '+('limited' if style=='-' else 'CW5-char/G4'))
                    for r,a in zip(rows,x):
                        b = r[key]
                        ax.plot(a,b,marker=marker,ms=5,color=COLORS[method],mfc=COLORS[method] if r['converged'] else 'white')
                    ax.set_xlabel('FV cells')
                    ax.set_xscale('log',base=2)
                    ax.set_xticks([60,120,240],labels=['60','120','240'])
        axes[0].set_ylabel(r'$p_w/p_\infty$')
        axes[1].set_ylabel(r'$C_f=2\tau_{r\theta,w}$')
        axes[2].set_ylabel(r'$q_{w,into\ fluid}/(\rho_\infty U_\infty^3)$')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,fontsize=8,frameon=False)
    save(fig,'grid-and-wall-transport','Grid sensitivity of pressure, shear and heat flux',
         'Filled: residual and drift gates passed · open: not accepted · o adiabatic / s isothermal')
    fig,axes = plt.subplots(1,2,figsize=(12,6))
    for ax,wall in zip(axes,['adiabatic','isothermal']):
        for r in visc_effective:
            c = r['config']
            if c['cells']!=60 or c['thermal_wall']!=wall or c['reconstruction']!='primitive-minmod':
                continue
            h = field(r)['history']
            method = c['flux']
            ax.semilogy(h[:,0],h[:,1],color=COLORS[method],lw=1.4,label=NAMES[method]+(' †' if not r['converged'] else ''))
            if len(h)==1:
                ax.plot(h[0,0],h[0,1],'o',color=COLORS[method])
        ax.axhline(2e-8,color='#384957',ls='--',lw=1)
        ax.set_xlabel('Newton attempts (explicit recoveries logged separately)')
        ax.set_ylabel('Maximum equation RMS residual')
        ax.set_title(wall.capitalize())
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False)
    save(fig,'steady-convergence','Convergence to the same strict residual gate',
         'Shared wall-matched AUSM+-up2 seed · seed is already converged for AUSM+-up2 · no speed ranking')
    fig,axes = plt.subplots(1,2,figsize=(12,6))
    for ax,mach in zip(axes,[2.35,7.95]):
        for r in euler:
            if r['config']['mach']!=mach or r['config']['cells']!=60:
                continue
            h = field(r)['history']
            method = r['config']['flux']
            ax.semilogy(h[:,0],h[:,1],color=COLORS[method],label=NAMES[method],ls='-' if r['converged'] else ':')
        if mach > 3:
            for r in char:
                h = field(r)['history']
                ax.semilogy(h[:,0],h[:,1],color=COLORS[r['config']['flux']],ls='--',lw=1.5,label=NAMES[r['config']['flux']]+' char')
        ax.axhline(2e-8,color='#384957',ls='--',lw=1)
        ax.set_xlabel('SSPRK3 steps')
        ax.set_ylabel('Maximum equation RMS residual')
        ax.set_title(f'Mach {mach:g}')
    fig.legend(*axes[1].get_legend_handles_labels(),loc='lower center',ncol=4,fontsize=8,frameon=False)
    save(fig,'euler-convergence','Flux choice and characteristic reconstruction',
         'Same CWENO3, 60 cells · dotted: not converged · dashed colored curves: characteristic CWENO3')
    fig,axes = plt.subplots(1,3,figsize=(13,5.5))
    for method in METHODS:
        data = np.load(OUT/f'sod-{method}.npz')
        for ax,k,label in zip(axes,[0,2,3],[r'$\rho$',r'$u$',r'$p$']):
            ax.plot(data['x'],data['primitive'][:,k],lw=1.2,color=COLORS[method],label=NAMES[method])
            ax.set_xlabel('x')
            ax.set_ylabel(label)
    for ax,k in zip(axes,[0,2,3]):
        ax.plot(data['x'],data['exact'][:,k],'k--',lw=1.5,label='Exact Riemann solution')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=5,fontsize=9,frameon=False)
    save(fig,'sod-comparison','Independent shock, contact and rarefaction benchmark',
         'Sod tube · 200 FV cells · t = 0.2 · all eight methods retain positive density and pressure')
    fig,axes = plt.subplots(1,3,figsize=(13,5.5))
    for ax,key in zip(axes,['nonuniform_cweno_rates','fv_diffusion_rates','smooth_viscous_angular_operator']):
        for name,rows in verification[key].items():
            ax.loglog([r['cells'] for r in rows],[r['l2_error'] for r in rows],'o-',label=f'{name} · last order {rows[-1]["order"]:.2f}')
        ax.set_xlabel('FV cells')
        ax.set_ylabel('L₂ operator error')
        ax.legend(fontsize=8)
    for ax,title in zip(axes,['Nonuniform CWENO','Interior FV diffusion','Angular viscous operator']):
        ax.set_title(title)
    save(fig,'spatial-verification','Measured smooth spatial accuracy',
         'Analytic cell averages and derivatives · wall closure is second order · full shock-case order is not inferred')
    # Export actual profiles for reuse; Kelvin temperature is only defined for viscous configurations.
    all_rows = euler+char+visc+continued
    for r in all_rows:
        d = field(r)
        q = d['primitive']
        vals = properties(q,r['config']['mach'])
        np.savetxt(OUT/(r['file']+'.csv'),np.c_[np.degrees(d['theta']),q[:,0],q[:,1],q[:,2],*vals[::3]],
                   delimiter=',',header='theta_deg,rho_ratio,vr_ratio,vtheta_ratio,pressure_ratio,temperature_ratio',comments='',fmt='%.12e')
    summary = {'euler_runs':len(euler),'euler_converged':sum(r['converged'] for r in euler),
               'viscous_effective_runs':len(visc_effective),'viscous_converged':sum(r['converged'] for r in visc_effective),
               'high_order_runs':sum(r['config']['gradient_order']==4 for r in visc_effective),
               'high_order_converged':sum(r['converged'] and r['config']['gradient_order']==4 for r in visc_effective),
               'verification_passed':verification['passed']}
    primary = [r for r in visc_effective if r['config']['cells']==60 and r['config']['reconstruction']=='primitive-minmod']
    summary.update(primary_viscous_runs=len(primary),primary_viscous_converged=sum(r['converged'] for r in primary))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    def table(rows,viscous=False):
        body = ''
        for r in sorted(rows,key=lambda a:(a['config'].get('thermal_wall',''),a['config']['mach'],a['config']['cells'],a['config']['reconstruction'],METHODS.index(a['config']['flux']))):
            c = r['config']
            status = '<span class="pass">Converged</span>' if r['converged'] else '<span class="fail">Not converged</span>'
            if r.get('steady_residual_converged') and not r.get('explicit_drift_check_passed',True):
                status = '<span class="fail">Explicit drift failed</span>'
            extra = f'<td>{r["skin_friction_coefficient"]:.6g}</td><td>{r["wall_heat_flux_into_fluid"]:.6g}</td>' if viscous else f'<td>{r["wall_pressure_error_percent"]:.4f}%</td>'
            body += f'<tr><td><a href="{html.escape(r["file"])}.json">{NAMES[c["flux"]]}</a></td><td>{c.get("thermal_wall",c["mach"])}</td><td>{c["cells"]}</td><td>{c["reconstruction"]}</td><td>{status}</td><td>{r["residual_rms_max_equation"]:.2e}</td><td>{r["wall_pressure_ratio"]:.6f}</td>{extra}</tr>'
        return '<div class="scroll"><table><thead><tr><th>Flux</th><th>Case</th><th>Cells</th><th>Reconstruction</th><th>Status</th><th>Residual</th><th>pwall/p∞</th>'+('<th>Cf</th><th>qwall into fluid</th>' if viscous else '<th>Wall p error</th>')+'</tr></thead><tbody>'+body+'</tbody></table></div>'
    figures = [('near-wall-profiles','Velocity decreases to zero at a stationary wall. The adiabatic wall heats through viscous dissipation; the isothermal wall removes heat.'),
               ('grid-and-wall-transport','Grid refinement is part of the assessment. Open markers are not accepted steady solutions.'),
               ('steady-convergence','Shared coarse seeds make these convergence histories reproducible, but they favor the seed flux in iteration count. They are not runtime rankings.'),
               ('viscous-properties-adiabatic','Local angular profiles at an adiabatic cone station.'),('viscous-properties-isothermal','Local angular profiles at a cold isothermal cone station.'),
               ('euler-properties-m2.35','All eight component-wise CWENO3 runs converged on the Mach 2.35 grid.'),
               ('euler-properties-m7.95','Dotted inviscid curves did not meet the residual gate; their apparent pressure accuracy is not accepted validation.'),
               ('euler-convergence','Characteristic reconstruction cures the recorded SLAU2 and EC+LF coarse-grid inviscid failures; AUSM+-up2 still requires additional robustness work.'),
               ('sod-comparison','Independent planar exact-solution check separates flux correctness from cone geometry.'),
               ('spatial-verification','Order tests use exact control-volume averages and analytic derivatives. They establish smooth operator accuracy, not global shock or wall convergence order.')]
    gallery = ''.join(f'<section><h2>{name.replace("-"," ").capitalize()}</h2><p>{text}</p><a href="{name}.pdf">Vector PDF</a><img src="{name}.png" alt="{html.escape(text)}"></section>' for name,text in figures)
    page = f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>AUSM and viscous cone comparison</title><style>body{{margin:0;background:#f2f5f7;color:#223745;font:16px/1.65 system-ui}}main{{max-width:1220px;margin:auto;padding:34px}}header{{background:#173d50;color:white;padding:38px;border-radius:16px}}h1{{font-size:36px;line-height:1.2}}h2{{line-height:1.3}}section{{background:white;margin-top:24px;padding:28px;border-radius:12px}}img{{width:100%;display:block;margin-top:20px}}a{{color:#087782}}header a{{color:#8fdddd}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{padding:9px;border-bottom:1px solid #dce4e9;text-align:left}}th{{background:#edf3f5}}.scroll{{overflow:auto}}.pass{{color:#08724f}}.fail{{color:#ae3d40}}.cards{{display:flex;gap:18px;flex-wrap:wrap}}.cards div{{background:#edf4f6;padding:18px;border-radius:10px;flex:1}}code,pre{{background:#eff3f5}}pre{{padding:18px;overflow:auto;font-size:13px}}</style><main><header><p>PROJECT 1 · PYTHON / NUMPY · EXECUTED 7 OCTOBER 2026</p><h1>AUSM-family fluxes and viscous cone modernization</h1><p>Measured flow properties, heat transfer, wall shear, convergence and grid sensitivity.</p><a href="../README.md">Source and migration guide</a> · <a href="verification.json">Verification evidence</a> · <a href="../results/report.html">Earlier inviscid study</a></header>
    <section><div class="cards"><div><strong>{summary['euler_converged']}/{summary['euler_runs']}</strong><br>component-wise inviscid runs converged</div><div><strong>{summary['primary_viscous_converged']}/{summary['primary_viscous_runs']}</strong><br>primary 60-cell viscous comparison converged</div><div><strong>{summary['high_order_converged']}/{summary['high_order_runs']}</strong><br>CWENO5-characteristic / G4 runs converged</div></div><p>Including all refinement attempts, {summary['viscous_converged']} of {summary['viscous_effective_runs']} corrected-wall viscous states passed the residual gate.</p><h2>What was implemented</h2><p>AUSM (1993), AUSM+ (1996), AUSM+-up (2006), AUSM+-up2 and SLAU2 (2013), and a Chandrashekar entropy-conservative flux with local LF dissipation (2013). These are established research algorithms newly implemented in this course code; they are not algorithms invented in 2026. Roe and HLLC remain baselines. The same inviscid face routines serve Euler and viscous equations.</p><p>Viscous physics include Newtonian stresses, Sutherland viscosity, heat conduction, Pr = 0.72, Re = 420,000, no-slip walls and separate adiabatic/isothermal conditions. A quadratic fit to two primitive control-volume averages supplies the second-order wall gradient; source quadrature uses the same wall boundary polynomial. Interior fourth-order diffusion and characteristic CWENO3/5 are optional.</p></section>
    <section><h2>Physical scope and acceptance</h2><p>The original angular viscous model is a <strong>local fixed-radial-station approximation</strong>. Primitive radial gradients are suppressed while spherical stress geometry is retained. Re changes with radius in real flow; this solver does not predict downstream boundary-layer growth. The result is not full axisymmetric Navier–Stokes cone validation. No turbulence, transition, real-gas chemistry or external CFD software is included.</p><p>The inviscid reference is independently integrated Taylor–Maccoll flow. It is not a viscous reference. A run passes only when the maximum RMS residual over all four conservation equations is below 2×10⁻⁸. Every accepted Newton solution also passes a 20-step explicit SSPRK drift check. Constitutive checks, exact Sod flow, Couette heating and analytic smooth operators give independent evidence. Viscous benchmark validation against experiment or a full 2-D reference remains open.</p><p>The quadratic FV wall fit is formally second order with accurate primitive averages. G2 infers primitives directly from conservative means; G4 deconvolves them. This distinction can affect wall heat/shear accuracy; a coupled global wall order is not claimed. Early explicit and cold-Newton failures are preserved separately. The first exploratory wall derivative was first order; only filenames starting <code>fv2-steady-</code> or their continuations use the corrected quadratic FV wall closure. Raw exploratory metadata predates that correction and is not part of the accepted comparison. Nonconverged states remain visible with status markers.</p></section>
    {gallery}<section><h2>Inviscid run table</h2>{table(euler)}<h2>Corrected-wall viscous run table</h2>{table(visc_effective,True)}</section>
    <section><h2>Reproduce</h2><p>Run from the conical directory. NumPy and Matplotlib are the only third-party dependencies. Shared seeds and their configuration are included. Dense Newton matrices are for teaching-size meshes; this is not a production large-mesh performance claim.</p><pre>python verify_extension.py
python shock_tube.py
python extension_study.py euler --grids 60 120
python extension_study.py characteristic --grids 60
python extension_study.py steady-warm --grids 60
python extension_study.py steady-refine --grids 120 240
python extension_study.py steady-high-order --grids 60 120 240
python extension_report.py</pre><p>See <a href="../VISCOUS.md">equations, source audit and implementation choices</a>. Fields are supplied as NPZ and CSV; per-run JSON records all physical parameters, convergence, solver details and wall transport. The ZIP manifest records final source hashes. Concurrent seeded timings are not speed rankings.</p><h2>Primary references</h2><ul><li><a href="https://ntrs.nasa.gov/citations/19920016565">Liou and Steffen: original AUSM</a></li><li><a href="https://ntrs.nasa.gov/citations/19940020698">Liou: AUSM+</a></li><li><a href="https://doi.org/10.1016/j.jcp.2005.09.020">Liou 2006: AUSM+-up</a></li><li><a href="https://doi.org/10.1016/j.jcp.2013.02.046">Kitamura and Shima 2013: pressure flux, AUSM+-up2 / SLAU2</a></li><li><a href="https://arxiv.org/abs/1209.4994">Chandrashekar: entropy and kinetic energy preserving fluxes</a></li><li><a href="https://arxiv.org/abs/1607.07319">Cravero et al.: CWENO finite-volume reconstruction</a></li><li><a href="https://github.com/su2code/SU2/blob/v8.5.0/SU2_CFD/src/numerics/flow/convection/ausm_slau.cpp">Pinned SU2 implementation: formula cross-check only; source not copied or solver run</a></li></ul></section></main></html>'''
    (OUT/'report.html').write_text(page,encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
