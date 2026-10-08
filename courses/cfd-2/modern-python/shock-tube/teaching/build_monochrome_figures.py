"""Lecture-only Times New Roman plots from verified, existing numerical records.

No solver execution, smoothing, interpolation of numerical profiles, or changes
to notebook figures. Lines join recorded cell means; the exact curve is analytic.
"""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT.parent / 'conical'
sys.path.insert(0, str(CORE))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from uniform_dashboard import load_uniform, reference_curve, selected
from uniform_refinement import METHODS, GRIDS, features
from shock_dashboard import NAMES

OUT = ROOT / 'teaching/figures'
OUT.mkdir(parents=True, exist_ok=True)
for name in ['times.ttf', 'timesbd.ttf', 'timesi.ttf', 'timesbi.ttf']:
    font_manager.fontManager.addfont('C:/Windows/Fonts/' + name)
plt.rcParams.update({'font.family': 'Times New Roman', 'font.size': 11,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Times New Roman',
    'mathtext.it': 'Times New Roman:italic', 'mathtext.bf': 'Times New Roman:bold',
    'mathtext.sf': 'Times New Roman', 'mathtext.tt': 'Times New Roman',
    'mathtext.cal': 'Times New Roman', 'mathtext.fallback': None,
    'axes.labelsize': 11, 'axes.titlesize': 12, 'xtick.labelsize': 9,
    'ytick.labelsize': 9, 'legend.fontsize': 10, 'axes.linewidth': .65,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'text.color': 'black', 'axes.labelcolor': 'black', 'axes.edgecolor': 'black'})
assert 'times.ttf' in font_manager.findfont('Times New Roman').lower()
rows = load_uniform()
xref, qref = reference_curve()
dash = [':', '--', '-.', (0, (5, 1, 1, 1)), (0, (2, 1)), '-']
markers = ['o', 's', '^', 'v', 'D', None]
outputs = []

def profile(ax, flux, k, transform=None):
    transform = transform or (lambda q: q[:, k])
    ax.plot(xref, transform(qref), color='black', lw=1.8,
            ls=(0, (8, 3)), label='Exact', zorder=1)
    for j, row in enumerate(selected(rows, flux)):
        ax.plot(row['x'], transform(row['primitive']), color='black',
                lw=1 if j < 5 else 1.25, ls=dash[j], marker=markers[j],
                ms=2.6, markerfacecolor='white', markeredgewidth=.6,
                markevery=max(1, len(row['x']) // 35),
                label=f'N = {row["config"]["cells"]:,}', zorder=2+j)
    ax.set_xlabel('Position x')

def finish(fig, axes, name, ncol=4):
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=ncol,
               frameon=False, handlelength=3, columnspacing=1.5)
    path = OUT / name
    fig.savefig(path, dpi=240, facecolor='white')
    plt.close(fig)
    outputs.append({'path': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})

transforms = [lambda q: q[:, 0], lambda q: q[:, 2], lambda q: q[:, 3],
    lambda q: q[:, 3]/(.4*q[:, 0]),
    lambda q: abs(q[:, 2])/np.sqrt(1.4*q[:, 3]/q[:, 0]),
    lambda q: q[:, 3]/q[:, 0]]
fig, axes = plt.subplots(2, 3, figsize=(8, 6), layout='constrained')
for ax, fn, title in zip(axes.flat, transforms,
        ['Density', 'Velocity', 'Pressure', 'Specific internal energy', 'Mach number', 'Temperature (R = 1)']):
    profile(ax, 'hllc', 0, fn)
    ax.set_title(title)
finish(fig, axes, 'uniform-physical-properties.png')

fig, axes = plt.subplots(2, 2, figsize=(8, 5.7), layout='constrained')
feature = features()
for ax, (name, k, title) in zip(axes.flat,
        [('contact', 0, 'Contact: density'), ('shock', 0, 'Shock: density'),
         ('shock', 2, 'Shock: velocity'), ('shock', 3, 'Shock: pressure')]):
    profile(ax, 'hllc', k)
    ax.set(title=title, xlim=(feature[name]['x']-.04, feature[name]['x']+.04))
    if k == 0: ax.set_ylim(feature[name]['low']-.015, feature[name]['high']+.015)
    elif k == 2: ax.set_ylim(-.04, 1.01)
    else: ax.set_ylim(.085, .32)
finish(fig, axes, 'uniform-shock-contact-detail.png')

shock = feature['shock']
for start in range(0, len(METHODS), 4):
    fig, axes = plt.subplots(2, 2, figsize=(8, 7), layout='constrained')
    for ax, flux in zip(axes.flat, METHODS[start:start+4]):
        profile(ax, flux, 0)
        ax.set(title=NAMES[flux], ylabel='Density',
               xlim=(shock['x']-.035, shock['x']+.035), ylim=(.11, .28))
    finish(fig, axes, f'uniform-all-method-shock-detail-{start//4+1}.png')

fig, axes = plt.subplots(1, 3, figsize=(8, 5.3), layout='constrained')
symbols = ['o','s','^','v','D','<','>','p','h','x','+','*','1','2','3','4']
for j, flux in enumerate(METHODS):
    chosen = selected(rows, flux)
    for ax, key in zip(axes, ['density', 'velocity', 'pressure']):
        ax.loglog([r['config']['cells'] for r in chosen],
            [r['errors'][key]['L1'] for r in chosen],
            color='black', lw=.85, ls=dash[j % 6], marker=symbols[j],
            ms=4, markerfacecolor='white', markeredgewidth=.7, label=NAMES[flux])
for ax, title in zip(axes, ['Density', 'Velocity', 'Pressure']):
    ax.set(title=title, xlabel='Uniform cell count N', ylabel='L1 error',
           xticks=[80, 320, 1280, 2560], xticklabels=['80','320','1280','2560'])
    ax.tick_params(axis='x', labelrotation=35)
finish(fig, axes, 'uniform-grid-error-convergence.png')

(OUT/'figure-provenance.json').write_text(json.dumps({
    'font': 'Times New Roman', 'style': 'monochrome',
    'source_catalog_sha256': hashlib.sha256((CORE/'results-uniform/catalog.json').read_bytes()).hexdigest(),
    'verified_records': len(rows), 'simulations_rerun': False,
    'numerical_profiles': 'Unmodified recorded cell means connected by straight segments',
    'outputs': outputs}, indent=2)+'\n')
print(json.dumps({'verified_records': len(rows), 'figures': len(outputs), 'font': 'Times New Roman'}))
