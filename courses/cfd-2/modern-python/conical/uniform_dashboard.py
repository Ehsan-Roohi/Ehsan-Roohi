"""Continuous cell-center profiles and measured fixed-grid convergence."""
import hashlib
import json
import numpy as np
from shock_suite import ROOT, CASES
from additional_fluxes import riemann_sample
from uniform_refinement import GRIDS, METHODS, features
from shock_dashboard import NAMES, style, finish, show_table

COLORS=['#a3adb8','#7890ab','#457b9d','#2a9d8f','#e9a33b','#c63243']


def load_uniform():
    target=ROOT/'results-uniform'; catalog=json.loads((target/'catalog.json').read_text());rows=[]
    for entry in catalog['runs']:
        for suffix in ['json','npz']:
            assert hashlib.sha256((target/(entry['name']+'.'+suffix)).read_bytes()).hexdigest()==entry[suffix+'_sha256']
        row=json.loads((target/(entry['name']+'.json')).read_text())
        with np.load(target/(entry['name']+'.npz')) as data:row.update({k:data[k].copy() for k in data.files})
        rows.append(row)
    rows.sort(key=lambda r:(METHODS.index(r['config']['flux']),r['config']['cells']))
    return rows


def reference_curve():
    # Include the analytic jumps explicitly so the reference is vertical there.
    locations=[f['x'] for f in features().values()]
    x=np.sort(np.r_[np.linspace(0,1,12001),*[p+d for p in locations for d in [-1e-12,1e-12]]])
    return x,riemann_sample(CASES['sod']['left'],CASES['sod']['right'],(x-.5)/.2)


def selected(rows,flux):
    return sorted([r for r in rows if r['config']['flux']==flux],key=lambda r:r['config']['cells'])


def colors_for(rows):
    import matplotlib.pyplot as plt
    return COLORS if len(rows)==6 else plt.cm.viridis(np.linspace(.15,.9,len(rows)))


def properties(rows,flux='hllc'):
    import matplotlib.pyplot as plt
    style();fig,axes=plt.subplots(2,3,figsize=(16,10),layout='constrained')
    def values(q):
        return [q[:,0],q[:,2],q[:,3],q[:,3]/(.4*q[:,0]),abs(q[:,2])/np.sqrt(1.4*q[:,3]/q[:,0]),q[:,3]/q[:,0]]
    chosen=selected(rows,flux)
    for r,color in zip(chosen,colors_for(chosen)):
        for ax,v in zip(axes.flat,values(r['primitive'])):
            ax.plot(r['x'],v,color=color,lw=1.4,label=f'N = {r["config"]["cells"]:,}')
    x,q=reference_curve()
    for ax,v,title in zip(axes.flat,values(q),['Density','Velocity','Pressure','Specific internal energy','Mach number','Temperature / R = 1']):
        ax.plot(x,v,'k--',lw=1.3,label='Exact');ax.set(xlabel='Position x',ylabel=title)
    h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=7,frameon=False)
    fig.suptitle('Sod · '+NAMES[flux]+' · fixed uniform-grid refinement',fontsize=18,fontweight='bold')
    return finish(fig,'uniform-physical-properties')


def details(rows,flux='hllc'):
    import matplotlib.pyplot as plt
    style();fig,axes=plt.subplots(2,2,figsize=(15,10),layout='constrained')
    specifications=[('contact',0,'Contact · density'),('shock',0,'Shock · density'),('shock',2,'Shock · velocity'),('shock',3,'Shock · pressure')]
    x,q=reference_curve();feature=features()
    chosen=selected(rows,flux)
    for ax,(name,k,title) in zip(axes.flat,specifications):
        for r,color in zip(chosen,colors_for(chosen)):
            ax.plot(r['x'],r['primitive'][:,k],color=color,lw=1.5,label=f'N = {r["config"]["cells"]:,}')
        ax.plot(x,q[:,k],'k--',lw=1.4,label='Exact')
        ax.set(xlim=(feature[name]['x']-.04,feature[name]['x']+.04),xlabel='Position x',ylabel=['Density','','Velocity','Pressure'][k],title=title)
        if k==0:ax.set_ylim(feature[name]['low']-.015,feature[name]['high']+.015)
        elif k==2:ax.set_ylim(-.04,1.01)
        else:ax.set_ylim(.085,.32)
    h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=7,frameon=False)
    fig.suptitle(NAMES[flux]+' · thinner numerical transitions as Δx decreases',fontsize=18,fontweight='bold')
    return finish(fig,'uniform-shock-contact-detail')


def cards(rows):
    import matplotlib.pyplot as plt
    style();fig,axes=plt.subplots(4,4,figsize=(16,13),layout='constrained')
    x,q=reference_curve();f=features()['shock']
    for ax,flux in zip(axes.flat,METHODS):
        chosen=selected(rows,flux)
        for r,color in zip(chosen,colors_for(chosen)):
            ax.plot(r['x'],r['primitive'][:,0],color=color,lw=1.3,label=f'N = {r["config"]["cells"]:,}')
        ax.plot(x,q[:,0],'k--',lw=1.2,label='Exact')
        ax.set(xlim=(f['x']-.035,f['x']+.035),ylim=(.11,.28),xlabel='Position x',ylabel='Density',title=NAMES[flux])
    h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=7,frameon=False)
    fig.suptitle('Shock detail · all 15 fluxes and JST · fixed uniform grids',fontsize=18,fontweight='bold')
    return finish(fig,'uniform-all-method-shock-detail')


def convergence(rows):
    import matplotlib.pyplot as plt
    style();fig,axes=plt.subplots(1,3,figsize=(17,7),layout='constrained')
    for flux,color in zip(METHODS,plt.cm.turbo(np.linspace(.04,.94,len(METHODS)))):
        chosen=selected(rows,flux)
        if not chosen:continue
        for ax,k in zip(axes,['density','velocity','pressure']):
            ax.loglog([r['config']['cells'] for r in chosen],[r['errors'][k]['L1'] for r in chosen],'-o',ms=3,lw=1.2,color=color,label=NAMES[flux])
    for ax,k in zip(axes,['Density','Velocity','Pressure']):ax.set(xlabel='Uniform cell count N',ylabel=k+' L1 error',xticks=GRIDS,xticklabels=[str(n) for n in GRIDS])
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=4,frameon=False,fontsize=9)
    fig.suptitle('Actual mesh convergence · Δx = 1/N · same CFL and final time',fontsize=18,fontweight='bold')
    return finish(fig,'uniform-grid-error-convergence')


def widths(rows):
    import matplotlib.pyplot as plt
    style();fig,axes=plt.subplots(1,2,figsize=(15,7),layout='constrained')
    for flux,color in zip(METHODS,plt.cm.turbo(np.linspace(.04,.94,len(METHODS)))):
        chosen=selected(rows,flux)
        if not chosen:continue
        for ax,k in zip(axes,['contact','shock']):
            pairs=[(r['config']['cells'],r['jump_widths'][k]['width_10_90']) for r in chosen if r['jump_widths'][k]['width_10_90'] is not None]
            if pairs:ax.loglog(*zip(*pairs),'-o',ms=3,lw=1.2,color=color,label=NAMES[flux])
    for ax,k in zip(axes,['Contact','Shock']):ax.set(xlabel='Uniform cell count N',ylabel='10–90% density transition width',title=k,xticks=GRIDS,xticklabels=[str(n) for n in GRIDS])
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=4,frameon=False,fontsize=9)
    fig.suptitle('Measured rounding decreases through mesh refinement',fontsize=18,fontweight='bold')
    return finish(fig,'uniform-transition-width-convergence')
