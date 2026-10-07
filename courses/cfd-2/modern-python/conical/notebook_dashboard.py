"""Show every method's executed fields, residuals and acceptance status in Colab."""
from pathlib import Path
import html
import json
import numpy as np
from flux import FLUXES

ROOT = Path(__file__).resolve().parent
NAMES = dict(zip(FLUXES,['Roe','Van Leer','HLLC','HLLE','Rusanov','AUSM','AUSM+',
                       'AUSM+-up','AUSM+-up2','SLAU2','EC + LF']))
PALETTE = ['#334155','#ca8a04','#2563eb','#7c3aed','#be185d','#c2410c',
           '#0891b2','#dc2626','#047857','#65a30d','#ea580c']


def load_study():
    catalog = json.loads((ROOT/'results-notebook/catalog.json').read_text(encoding='utf-8'))
    if not catalog['complete']:
        raise RuntimeError('Published study is incomplete')
    runs = []
    for item in catalog['runs']:
        row = json.loads((ROOT/item['json_path']).read_text(encoding='utf-8'))
        with np.load(ROOT/item['npz_path']) as data:
            row.update({k:data[k].copy() for k in ('theta','primitive','conserved_average','history')})
        row.update(group=item['group'],record_path=item['json_path'])
        runs.append(row)
    return runs


def select(runs,group='flux',physics='inviscid',mach=None):
    return [r for r in runs if r['group']==group and
            (('thermal_wall' not in r['config']) if physics=='inviscid' else r['config'].get('thermal_wall')==physics) and
            (mach is None or r['config']['mach']==mach)]


def table_html(runs):
    rows = []
    viscous = bool(runs and 'thermal_wall' in runs[0]['config'])
    for r in runs:
        c = r['config']
        status = 'Accepted' if r['converged'] else 'Gate not met'
        color = '#08704f' if r['converged'] else '#b23d34'
        extra = (f'<td>{r["skin_friction_coefficient"]:.6g}</td><td>{r["wall_heat_flux_into_fluid"]:.6g}</td>' if viscous else '')
        reconstruction = 'DG degree '+str(c['degree']) if 'degree' in c else c['reconstruction']
        rows.append(f'<tr><td>{html.escape(NAMES[c["flux"]])}</td><td>{html.escape(reconstruction)}</td>'
                    f'<td>{c["cells"]}</td><td style="color:{color};font-weight:bold">{status}</td>'
                    f'<td>{r["residual_rms_max_equation"]:.3e}</td><td>{r["wall_pressure_ratio"]:.7f}</td>{extra}</tr>')
    return ('<div style="overflow:auto"><table style="border-collapse:collapse;font:13px system-ui;width:100%" cellpadding="8">'
            '<thead style="background:#eaf2f6"><tr><th>Flux</th><th>Reconstruction</th><th>Cells</th><th>Status</th>'
            '<th>Final residual</th><th>p wall / p∞</th>'+('<th>Cf</th><th>q wall into fluid</th>' if viscous else '')+
            '</tr></thead><tbody>'+''.join(rows)+'</tbody></table></div>')


def show_table(runs):
    text = table_html(runs)
    try:
        from IPython.display import HTML, display
    except ImportError:
        print('Methods:', ', '.join(NAMES[r['config']['flux']] for r in runs))
        print('Accepted:',sum(r['converged'] for r in runs),'/',len(runs))
    else:
        display(HTML(text))
    return text


def style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'axes.grid':True,'grid.alpha':.15,'pdf.fonttype':42})


def label(row, group):
    c = row['config']
    if group=='reconstruction':
        name = c['reconstruction'].upper().replace('MUSCL','MUSCL-MC')
    elif group=='dg':
        name = NAMES[c['flux']]+' / degree '+str(c['degree'])
    elif group=='viscous-high-order':
        name = NAMES[c['flux']]+f' / N={c["cells"]}'
    else:
        name = NAMES[c['flux']]
    return name+(' †' if not row['converged'] else '')


def color(row,index,group):
    import matplotlib.pyplot as plt
    return PALETTE[list(FLUXES).index(row['config']['flux'])] if group=='flux' else plt.cm.tab10(index%10)


def finish(fig,name,title,subtitle,runs,group):
    import matplotlib.pyplot as plt
    out = ROOT/'notebook-figures'
    out.mkdir(exist_ok=True)
    fig.suptitle(title,fontsize=19,fontweight='bold',y=.985)
    fig.text(.5,.927,subtitle,ha='center',fontsize=10,color='#526473')
    handles,labels = fig.axes[0].get_legend_handles_labels()
    if name=='all-eleven-methods':
        handles,labels = handles[:1],labels[:1]
    fig.legend(handles,labels,loc='lower center',ncol=4,frameon=False,fontsize=9)
    fig.tight_layout(rect=[.02,.14,.98,.9],h_pad=2,w_pad=2)
    fig.savefig(out/(name+'.png'),dpi=145)
    fig.savefig(out/(name+'.pdf'))
    plt.show()
    plt.close(fig)
    return fig


def plot_properties(runs,name,title,group='flux'):
    import matplotlib.pyplot as plt
    style()
    fig,axes = plt.subplots(2,3,figsize=(13.5,8.2))
    if not runs:
        raise ValueError('No runs selected')
    c = runs[0]['config']
    if 'thermal_wall' not in c:
        from reference import taylor_maccoll,sample_reference
        from state import primitive
        ref = taylor_maccoll(c['mach'],max_step=1e-4)
        angles = np.radians(np.linspace(c['cone_deg'],c['outer_deg'],800))
        q = primitive(sample_reference(angles,ref))
        vals = properties(q,c)
        for ax,val in zip(axes.flat,vals):
            ax.plot(np.degrees(angles),val,'k--',lw=1.6,label='Taylor–Maccoll')
    for index,r in enumerate(runs):
        for ax,val in zip(axes.flat,properties(r['primitive'],r['config'])):
            ax.plot(np.degrees(r['theta']),val,ls='-' if r['converged'] else ':',
                    color=color(r,index,group),lw=1.5,label=label(r,group))
    for ax,caption in zip(axes.flat,[r'$p/p_\infty$',r'$\rho/\rho_\infty$','Local Mach number',
                                    r'$T/T_\infty$',r'$v_r/U_\infty$',r'$v_\theta/U_\infty$']):
        ax.set(xlabel='Polar angle θ [deg]',ylabel=caption)
    return finish(fig,name,title,f'{len(runs)} executed methods · solid = accepted · dotted † = acceptance gate not met',runs,group)


def properties(q,c):
    pressure = q[:,3]*c['gamma']*c['mach']**2
    return [pressure,q[:,0],np.hypot(q[:,1],q[:,2])/np.sqrt(c['gamma']*q[:,3]/q[:,0]),
            pressure/q[:,0],q[:,1],q[:,2]]


def plot_convergence(runs,name,title,group='flux'):
    import matplotlib.pyplot as plt
    style()
    fig,axes = plt.subplots(1,2,figsize=(13.5,5.8))
    for index,r in enumerate(runs):
        h = np.asarray(r['history'])
        col = color(r,index,group)
        axes[0].semilogy(h[:,0],np.maximum(h[:,1],1e-18),color=col,label=label(r,group),lw=1.4)
        if len(h)==1:
            axes[0].plot(h[:,0],h[:,1],'o',color=col)
        axes[1].semilogy(index,max(r['residual_rms_max_equation'],1e-18),'o',mfc=col if r['converged'] else 'white',mec=col,ms=7)
    for ax in axes:
        ax.axhline(2e-8,color='black',ls='--',lw=1)
        ax.set_ylabel('Maximum RMS residual')
    newton = bool(runs and runs[0].get('steady_solver'))
    axes[0].set_xlabel('Newton attempts; explicit recoveries logged separately' if newton else 'SSPRK steps')
    axes[1].set_xticks(range(len(runs)),[label(r,group) for r in runs],rotation=65,ha='right',fontsize=8)
    axes[1].set_title('Final residual: every method is shown')
    return finish(fig,name,title,'Residual gate = 2 × 10⁻⁸ · seeded Newton histories are not speed rankings' if newton
                  else 'Residual gate = 2 × 10⁻⁸ · full method histories; failures are retained',runs,group)


def plot_wall_transport(runs,name,title):
    import matplotlib.pyplot as plt
    style()
    fig,axes = plt.subplots(1,3,figsize=(13.5,5.6))
    for ax,key,caption in zip(axes,['skin_friction_coefficient','wall_heat_flux_into_fluid','wall_temperature_ratio'],
                              ['Skin friction Cf','Heat flux into fluid','Wall temperature T/T∞']):
        for i,r in enumerate(runs):
            ax.plot(i,r[key],'o',color=color(r,i,'flux'),ms=6)
        ax.set_xticks(range(len(runs)),[NAMES[r['config']['flux']] for r in runs],rotation=65,ha='right',fontsize=8)
        ax.set_ylabel(caption)
    return finish(fig,name,title,'Cold wall: negative heat flux means heat leaves the fluid · local-station angular model',runs,'flux')


def plot_method_cards(runs,name='all-eleven-methods'):
    import matplotlib.pyplot as plt
    from reference import taylor_maccoll,sample_reference
    from state import primitive
    style()
    fig,axes = plt.subplots(3,4,figsize=(14,10.5))
    c = runs[0]['config']
    angle = np.radians(np.linspace(c['cone_deg'],c['outer_deg'],700))
    reference = primitive(sample_reference(angle,taylor_maccoll(c['mach'],max_step=1e-4)))[:,3]*c['gamma']*c['mach']**2
    for i,(ax,r) in enumerate(zip(axes.flat,runs)):
        ax.plot(np.degrees(angle),reference,'k--',lw=1,label='Taylor–Maccoll')
        ax.plot(np.degrees(r['theta']),properties(r['primitive'],r['config'])[0],color=color(r,i,'flux'),lw=1.6,label=NAMES[r['config']['flux']])
        ax.set_title(NAMES[r['config']['flux']],color=color(r,i,'flux'))
        ax.set(xlabel='θ [deg]',ylabel='p/p∞')
        ax.text(.04,.08,f'R = {r["residual_rms_max_equation"]:.2e}\n'+('Accepted' if r['converged'] else 'Gate not met'),transform=ax.transAxes,fontsize=8)
    for ax in list(axes.flat)[len(runs):]:
        ax.axis('off')
    return finish(fig,name,'Every numerical flux has a result',
                  '11 separate pressure panels · common 60-cell CWENO3 case · Mach 2.35',runs,'flux')
