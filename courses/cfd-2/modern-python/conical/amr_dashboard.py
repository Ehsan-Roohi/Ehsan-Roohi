"""Physical cell-average plots and measured local-refinement comparisons."""
import json,hashlib
import numpy as np
from shock_amr import ROOT,CASES,sod_features,riemann_sample
from notebook_dashboard import NAMES
from shock_dashboard import style,finish


def load_amr():
    target=ROOT/'results-amr';catalog=json.loads((target/'catalog.json').read_text());rows=[]
    for entry in catalog['runs']:
        for suffix in ['json','npz']:
            p=target/(entry['name']+'.'+suffix)
            assert hashlib.sha256(p.read_bytes()).hexdigest()==entry[suffix+'_sha256']
        row=json.loads((target/(entry['name']+'.json')).read_text())
        with np.load(target/(entry['name']+'.npz')) as data:row.update({k:data[k].copy() for k in data.files})
        rows.append(row)
    return rows


def label(r):
    c=r['config']
    return f'AMR base {c["base_cells"]}, level {c["max_level"]}' if c['max_level'] else f'Uniform {c["base_cells"]}'


def table_html(rows):
    text='<details><summary>All '+str(len(rows))+' mesh comparisons (click to expand)</summary><table style="border-collapse:collapse;font-size:12px"><tr>'
    headings=['Flux','Mesh','Final / mean cells','Density L1','Contact width','Shock width','Conservation defect','Status']
    text+=''.join('<th style="padding:6px;background:#e7eff7">'+h+'</th>' for h in headings)+'</tr>'
    for r in rows:
        widths=r['jump_widths']
        values=[NAMES[r['config']['flux']],label(r),f'{r["active_cells"]} / {r["mean_active_cells"]:.1f}',
                f'{r.get("errors",{}).get("density",{}).get("L1",float("nan")):.5g}',
                *[f'{widths[k]["width_10_90"]:.5g}' if widths.get(k,{}).get('width_10_90') is not None else 'Not resolved' for k in ['contact','shock']],
                f'{r["conservation_budget_defect"]:.2e}','Passed' if r['integrity_passed'] else r['termination']]
        text+='<tr>'+''.join('<td style="padding:5px;border-bottom:1px solid #ddd">'+str(v)+'</td>' for v in values)+'</tr>'
    return text+'</table></details>'


def show_amr_table(rows):
    from IPython.display import display,HTML
    text=table_html(rows);display(HTML(text));return text


def profiles(rows,flux='hllc'):
    import matplotlib.pyplot as plt
    style();fig,axes=plt.subplots(2,3,figsize=(16,10),layout='constrained')
    chosen=sorted([r for r in rows if r['config']['flux']==flux],key=lambda r:(r['config']['max_level'],r['config']['base_cells']))
    x=np.linspace(0,1,2401);exact=riemann_sample(CASES['sod']['left'],CASES['sod']['right'],(x-.5)/.2)
    colors=['#8d99ae','#457b9d','#2a9d8f','#e9c46a','#e63946']
    for r,color in zip(chosen,colors):
        for ax,k in zip(axes[0],[0,2,3]):ax.stairs(r['primitive'][:,k],r['faces'],baseline=None,color=color,lw=1.25,label=label(r))
        for ax in axes[1,:2]:ax.stairs(r['primitive'][:,0],r['faces'],baseline=None,color=color,lw=1.25)
        axes[1,2].stairs(np.diff(r['faces']),r['faces'],baseline=None,color=color,lw=1.25)
    for ax,k,title in zip(axes[0],[0,2,3],['Density','Velocity','Pressure']):
        ax.plot(x,exact[:,k],'k--',lw=1.5,label='Exact pointwise solution');ax.set(xlabel='Position x',ylabel=title)
    for ax,name in zip(axes[1,:2],['contact','shock']):
        f=sod_features()[name]
        ax.plot(x,exact[:,0],'k--',lw=1.5);ax.set(xlim=(f['x']-.055,f['x']+.055),ylim=(f['low']-.02,f['high']+.02),xlabel='Position x',ylabel='Density',title=name.title()+' detail · raw FV stairs')
    axes[1,2].set(xlabel='Position x',ylabel='Cell width Δx',yscale='log',title='Actual final cell widths')
    h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle('Sod · '+NAMES[flux]+' · local refinement sharpens resolved jumps',fontsize=18,fontweight='bold')
    return finish(fig,'amr-physical-properties')


def mesh_evolution(rows):
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import BoundaryNorm
    style();r=next(r for r in rows if r['config']['flux']=='hllc' and r['config']['max_level'])
    fig,axes=plt.subplots(2,1,figsize=(15,8),layout='constrained');polygons=[];levels=[];times=[]
    for i in range(5):
        faces=r[f'snapshot_{i}_faces'];t=float(r[f'snapshot_{i}_time']);times.append(t)
        level=np.rint(np.log2(1/(r['config']['base_cells']*np.diff(faces)))).astype(int)
        for a,b,k in zip(faces[:-1],faces[1:],level):polygons.append([(a,t-.014),(b,t-.014),(b,t+.014),(a,t+.014)]);levels.append(k)
    norm=BoundaryNorm(np.arange(-.5,4.5),4)
    artist=PolyCollection(polygons,array=np.array(levels),cmap=plt.get_cmap('viridis',4),norm=norm,edgecolors='#ffffff',linewidths=.12)
    axes[0].add_collection(artist);axes[0].set(xlim=(0,1),ylim=(-.02,.22),yticks=times,xlabel='Position x',ylabel='Actual snapshot time',title='Refined cells follow numerical density and pressure gradients')
    axes[0].grid(False);fig.colorbar(artist,ax=axes[0],ticks=[0,1,2,3],label='Dyadic refinement level')
    history=r['mesh_history'];axes[1].plot(history[:,0],history[:,1],color='#e63946',lw=2)
    axes[1].axhline(640,color='#457b9d',ls='--',label='Uniform 640 control');axes[1].set(xlabel='Physical time',ylabel='Active cells',title='Actual mesh size · global stepping · conservative split/merge');axes[1].legend(frameon=False)
    fig.suptitle('Adaptive mesh history · no exact-solution mesh flags',fontsize=18,fontweight='bold')
    return finish(fig,'amr-mesh-evolution')


def all_flux_cards(rows):
    import matplotlib.pyplot as plt
    from flux import FLUXES
    style();fig,axes=plt.subplots(4,4,figsize=(16,12),layout='constrained')
    f=sod_features()['shock'];x=np.linspace(f['x']-.045,f['x']+.045,801)
    exact=riemann_sample(CASES['sod']['left'],CASES['sod']['right'],(x-.5)/.2)
    for ax,flux in zip(axes.flat,FLUXES):
        for r in rows:
            if r['config']['flux']!=flux:continue
            if r['config']['max_level']:color='#e63946';title='AMR'
            elif r['config']['base_cells']==80:color='#457b9d';title='Uniform 80'
            else:continue
            ax.stairs(r['primitive'][:,0],r['faces'],baseline=None,color=color,lw=1.5,label=title)
        ax.plot(x,exact[:,0],'k--',lw=1.2,label='Exact');ax.set(xlim=(x[0],x[-1]),ylim=(.105,.29),xlabel='x',ylabel='Density',title=NAMES[flux])
    axes.flat[-1].axis('off');axes.flat[-1].text(.06,.7,'Blue: uniform 80\nRed: locally refined AMR\nBlack: exact\n\nRaw cell averages shown as stairs.\nNo image sharpening or smoothing.',transform=axes.flat[-1].transAxes,fontsize=12)
    fig.suptitle('Shock detail · every shared flux · genuine local mesh refinement',fontsize=18,fontweight='bold')
    return finish(fig,'amr-all-flux-shock-detail')


def accuracy(rows):
    import matplotlib.pyplot as plt
    from flux import FLUXES
    style();fig,axes=plt.subplots(1,3,figsize=(17,7),layout='constrained')
    for flux,color in zip(FLUXES,plt.cm.turbo(np.linspace(.04,.94,15))):
        selected=[r for r in rows if r['config']['flux']==flux and r['integrity_passed']]
        uniform=sorted([r for r in selected if not r['config']['max_level']],key=lambda r:r['config']['base_cells'])
        adaptive=[r for r in selected if r['config']['max_level']]
        for ax,k in zip(axes,['density','velocity','pressure']):
            ax.loglog([r['events']['rhs_cell_updates'] for r in uniform],[r['errors'][k]['L1'] for r in uniform],'-o',ms=3,color=color,lw=1.1,label=NAMES[flux])
            for r in adaptive:ax.loglog(r['events']['rhs_cell_updates'],r['errors'][k]['L1'],'*',ms=12,color=color,mec='black',mew=.4)
    for ax,k in zip(axes,['Density','Velocity','Pressure']):ax.set(xlabel='RHS cell updates (includes RK stages/retries)',ylabel=k+' width-weighted L1 error')
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=5,fontsize=8,frameon=False)
    fig.suptitle('Accuracy versus numerical RHS work · stars = AMR',fontsize=18,fontweight='bold')
    return finish(fig,'amr-error-versus-work')


def widths(rows):
    import matplotlib.pyplot as plt
    from flux import FLUXES
    style();fig,axes=plt.subplots(1,2,figsize=(15,10),layout='constrained');y=np.arange(15)
    for ax,k in zip(axes,['contact','shock']):
        old=[];new=[]
        for flux in FLUXES:
            old.append(next(r for r in rows if r['config']['flux']==flux and r['config']['base_cells']==80 and not r['config']['max_level'])['jump_widths'][k]['width_10_90'])
            new.append(next(r for r in rows if r['config']['flux']==flux and r['config']['max_level'])['jump_widths'][k]['width_10_90'])
        ax.barh(y-.18,old,height=.34,color='#457b9d',label='Uniform 80');ax.barh(y+.18,new,height=.34,color='#e63946',label='AMR base 80 / level 3')
        ax.set(yticks=y,yticklabels=[NAMES[f] for f in FLUXES],xlabel='Apparent 10–90% density jump width',title=k.title()+' · smaller means sharper');ax.invert_yaxis();ax.legend(frameon=False)
    fig.suptitle('Measured discontinuity thickness · cell-center threshold interpolation',fontsize=18,fontweight='bold')
    return finish(fig,'amr-jump-widths')
