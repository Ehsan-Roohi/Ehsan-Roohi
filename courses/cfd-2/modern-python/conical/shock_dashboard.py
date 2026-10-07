"""Executed shock-tube comparisons, scientific plots and native Colab tables."""
from pathlib import Path
import json,html
import numpy as np
ROOT=Path(__file__).resolve().parent
NAMES={'roe':'Roe','vanleer':'Van Leer','hllc':'HLLC','hlle':'HLL / Davis','rusanov':'Rusanov',
 'ausm':'AUSM','ausm-plus':'AUSM+','ausm-up':'AUSM+-up','ausm-up2':'AUSM+-up2','slau2':'SLAU2','ec-lf':'EC + LF',
 'godunov':'Godunov','roe-nc':'Roe / no fix','steger-warming':'Steger-Warming','global-lf':'Global LF','jst':'JST'}

def load():
    catalog=json.loads((ROOT/'results-shocktube/catalog.json').read_text())
    rows=[]
    for item in catalog['runs']:
        r=json.loads((ROOT/'results-shocktube'/(item['name']+'.json')).read_text())
        p=ROOT/'results-shocktube'/(item['name']+'.npz')
        if p.exists():
            with np.load(p) as data:r.update({k:data[k].copy() for k in data.files})
        r['name']=item['name'];rows.append(r)
    return rows

def select(rows,case='sod',cells=80,flux=None,reconstruction=None,degree=-1):
    return [r for r in rows if r['config']['case']==case and (cells is None or r['config']['cells']==cells)
            and (flux is None or r['config']['flux']==flux) and (reconstruction is None or r['config']['reconstruction']==reconstruction)
            and (degree is None or r['config']['degree']==degree)]

def table_html(rows,collapsible=False):
    body=[]
    for r in rows:
        c=r['config'];error=r.get('errors',{}).get('density',{}).get('L1')
        status='Completed / integrity passed' if r['integrity_passed'] else r.get('termination','Failed')
        body.append('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [NAMES[c['flux']],c['reconstruction'] if c['degree']<0 else 'DG '+str(c['degree']),c['cells'],
          f'{r.get("time",0):.5g}',status,f'{error:.5g}' if error is not None else 'Not eligible',f'{r.get("conservation_budget_defect",float("nan")):.2e}'])+'</tr>')
    text='<div style="overflow:auto"><table cellpadding="7" style="font:12px system-ui;border-collapse:collapse"><thead><tr>'+''.join('<th>'+h+'</th>' for h in ['Flux','Reconstruction','Cells','Actual time','Status','Density L1','Conservation defect'])+'</tr></thead><tbody>'+''.join(body)+'</tbody></table></div>'
    return '<details><summary>All '+str(len(rows))+' run records</summary>'+text+'</details>' if collapsible else text

def show_table(rows,collapsible=False):
    text=table_html(rows,collapsible)
    try:
        from IPython.display import display,HTML
        display(HTML(text))
    except ImportError: print('Runs:',len(rows),'completed:',sum(r['completed'] for r in rows))
    return text

def style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.15,'pdf.fonttype':42})

def label(r,by='flux'):
    c=r['config']
    text=NAMES[c['flux']] if by=='flux' else (c['reconstruction'] if c['degree']<0 else NAMES[c['flux']]+' / DG '+str(c['degree']))
    return text if r['completed'] else text+f' [stopped t={r.get("time",0):.3g}]'

def finish(fig,name):
    import matplotlib.pyplot as plt
    target=ROOT/'shocktube-figures';target.mkdir(exist_ok=True)
    fig.savefig(target/(name+'.png'),dpi=145,bbox_inches='tight');fig.savefig(target/(name+'.pdf'),bbox_inches='tight')
    plt.show();plt.close(fig);return fig

def properties(rows,name,title,by='flux'):
    import matplotlib.pyplot as plt
    style();fig,axes=plt.subplots(2,3,figsize=(15,9),layout='constrained');colors=plt.cm.turbo(np.linspace(.04,.94,len(rows)))
    def values(q):
        a=np.sqrt(1.4*q[:,3]/q[:,0]);return [q[:,0],q[:,2],q[:,3],q[:,3]/(.4*q[:,0]),abs(q[:,2])/a,q[:,3]/q[:,0]]
    ref=next(r for r in rows if 'exact' in r)
    for ax,v in zip(axes.flat,values(ref['exact'])):ax.plot(ref['x'],v,'k--',lw=1.8,label='Exact / final time')
    for r,color in zip(rows,colors):
        if 'primitive' not in r:continue
        for ax,v in zip(axes.flat,values(r['primitive'])):ax.plot(r['x'],v,color=color,ls='-' if r['completed'] else ':',lw=1.15,label=label(r,by))
    for ax,y in zip(axes.flat,['Density','Velocity','Pressure','Specific internal energy','Mach number','Temperature / R=1']):ax.set(xlabel='Position x',ylabel=y)
    h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=5,fontsize=8,frameon=False)
    fig.suptitle(title,fontsize=18,fontweight='bold')
    return finish(fig,name)

def cards(rows):
    import matplotlib.pyplot as plt
    style();fig,axes=plt.subplots(4,4,figsize=(15,12),layout='constrained')
    for ax,r,color in zip(axes.flat,rows,plt.cm.turbo(np.linspace(.04,.94,len(rows)))):
        if 'primitive' in r:
            ax.plot(r['x'],r['exact'][:,0],'k--',lw=1.1);ax.plot(r['x'],r['primitive'][:,0],color=color,lw=1.6)
        e=r.get('errors',{}).get('density',{}).get('L1')
        ax.set(title=NAMES[r['config']['flux']],xlabel='x',ylabel='Density')
        ax.text(.03,.06,f'L1 = {e:.4g}' if e is not None else 'Not completed',transform=ax.transAxes,fontsize=9)
    fig.suptitle('Sod shock tube · every flux has an individual result',fontsize=18,fontweight='bold')
    return finish(fig,'all-flux-cards')

def matrix(rows,fluxes,recons):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    style();fig,ax=plt.subplots(figsize=(17,10),layout='constrained')
    a=np.full((len(fluxes),len(recons)),np.nan)
    for r in rows:
        c=r['config']
        if c['flux'] in fluxes and c['reconstruction'] in recons and r['completed']:
            a[fluxes.index(c['flux']),recons.index(c['reconstruction'])]=r['errors']['density']['L1']
    im=ax.imshow(a,norm=LogNorm(vmin=np.nanmin(a),vmax=np.nanmax(a)),cmap='viridis_r',aspect='auto')
    for i in range(len(fluxes)):
        for j in range(len(recons)):ax.text(j,i,'FAIL' if np.isnan(a[i,j]) else f'{a[i,j]:.3f}',ha='center',va='center',fontsize=7,color='black' if np.isnan(a[i,j]) else 'white')
    ax.set_xticks(range(len(recons)),recons,rotation=55,ha='right',fontsize=9);ax.set_yticks(range(len(fluxes)),[NAMES[f] for f in fluxes]);ax.grid(False)
    ax.set_title('Sod · complete flux × reconstruction matrix · density L1 error',fontsize=17,pad=20)
    fig.colorbar(im,ax=ax,label='Mean absolute density error / lower is better')
    return finish(fig,'sod-full-matrix')

def convergence(rows,name,title,by='flux'):
    import matplotlib.pyplot as plt
    style();fig,axes=plt.subplots(1,3,figsize=(16,6),layout='constrained')
    keys=sorted(set(label(r,by) for r in rows if r['completed']))
    for key,color in zip(keys,plt.cm.turbo(np.linspace(.04,.94,len(keys)))):
        group=sorted([r for r in rows if label(r,by)==key and r['completed']],key=lambda r:r['config']['cells'])
        if len(group)<2:continue
        for ax,prop in zip(axes,['density','velocity','pressure']):ax.loglog([r['config']['cells'] for r in group],[r['errors'][prop]['L1'] for r in group],'-o',ms=4,lw=1.2,color=color,label=key)
    for ax,prop in zip(axes,['Density','Velocity','Pressure']):ax.set(xlabel='Finite-volume cells / elements',ylabel=prop+' L1 error')
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=5,fontsize=8,frameon=False)
    fig.suptitle(title,fontsize=17,fontweight='bold');return finish(fig,name)

def histories(rows):
    import matplotlib.pyplot as plt
    style();fig,axes=plt.subplots(1,3,figsize=(15,6),layout='constrained')
    for r,color in zip(rows,plt.cm.turbo(np.linspace(.04,.94,len(rows)))):
        hist=r.get('history',np.zeros((0,4)))
        if len(hist):
            for ax,k in zip(axes,[1,2,3]):ax.semilogy(hist[:,0],np.maximum(hist[:,k],1e-17),color=color,lw=1.1,label=label(r))
    for ax,y in zip(axes,['Minimum density','Minimum pressure','Conservation budget defect']):ax.set(xlabel='Physical time',ylabel=y)
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=5,fontsize=8,frameon=False)
    fig.suptitle('Transient integrity histories · positivity and conservation',fontsize=17,fontweight='bold');return finish(fig,'transient-integrity')
