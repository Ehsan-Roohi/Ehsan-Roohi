"""Render three course case figures from preserved numerical data (no CFD solve).

Run from any directory: python tools/render_case_figures.py
Requires NumPy and Matplotlib. Figure provenance and checks are saved beside plots.
"""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.patches import Polygon

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'case-studies'
BG, PANEL, TEXT, MUTED, ACCENT = '#081725', '#102536', '#eef6fa', '#9db6c7', '#52dfd0'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
    'text.color': TEXT, 'axes.labelcolor': MUTED, 'xtick.color': MUTED,
    'ytick.color': MUTED, 'axes.edgecolor': '#345064', 'axes.facecolor': PANEL,
    'figure.facecolor': BG, 'savefig.facecolor': BG})
evidence = {'purpose': 'Visualization of preserved data; no solver was executed by this script.', 'cases': {}}

def fingerprint(path):
    p = ROOT / path
    return {'path': path, 'bytes': p.stat().st_size,
            'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def canvas(number, title, subtitle, facts):
    fig = plt.figure(figsize=(17, 9.6), dpi=150)
    fig.text(.055, .95, f'CFD II  /  MAZAHERI                                      CASE {number:02d}',
             color=ACCENT, fontsize=11, weight='bold')
    fig.text(.055, .885, title, fontsize=30, weight='bold')
    fig.text(.055, .846, subtitle, fontsize=12, color=MUTED)
    for x, (label, value) in zip([.055, .335, .62], facts):
        fig.text(x, .785, value, fontsize=18, weight='bold')
        fig.text(x, .754, label.upper(), fontsize=9, color=MUTED)
    return fig

def style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.tick_params(labelsize=9)
    ax.spines[['top', 'right']].set_visible(False)

def colorbar(fig, artist, ax, label):
    cb = fig.colorbar(artist, ax=ax, fraction=.032, pad=.025, shrink=.88)
    cb.set_label(label, fontsize=10)
    cb.ax.tick_params(labelsize=9)
    cb.outline.set_edgecolor('#345064')

def save(fig, slug, note):
    fig.text(.055, .045, note, fontsize=10, color=MUTED)
    folder = OUT / slug; folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / 'figure.png', dpi=150)
    plt.close(fig)

def cone():
    source = 'results/01-conical-flow/ROE FINAL.plt'
    lines = (ROOT / source).read_text().splitlines()
    assert re.findall(r'"([^"]+)"', lines[0]) == ['TETA', 'MACH', 'UR', 'UT', 'MACHT', 'RHO', 'P', 'T', 'ANGLE']
    values = np.fromstring(' '.join(lines[1:]).replace('D', 'E'), sep=' ')
    assert values.size % 9 == 0
    a = values.reshape(-1, 9)
    assert np.isfinite(a).all() and (np.diff(a[:, 0]) > 0).all() and (a[:, 6] > 0).all()
    # TETA is exported as degrees above the cone surface, not absolute angle.
    # Display reconstructs a self-similar sector: pressure is constant on each ray.
    half_angle = 10.0
    theta = np.deg2rad(half_angle + a[:, 0])
    radius = np.linspace(.05, 3.15, 360)
    x = np.cos(theta[:, None]) * radius
    y = np.sin(theta[:, None]) * radius
    pressure = np.broadcast_to(a[:, 6, None], x.shape)
    fig = canvas(1, 'Conical flow', 'Compression around a cone, viewed through the preserved angular pressure profile.',
        [('stored angular samples', str(len(a))), ('display cone half-angle', '10\N{DEGREE SIGN}'),
         ('visualized quantity', 'Pressure ratio')])
    ax = fig.add_axes([.07, .20, .56, .48])
    levels = np.linspace(pressure.min(), pressure.max(), 80)
    cf = ax.contourf(x, y, pressure, levels=levels, cmap='turbo')
    ax.contourf(x, -y, pressure, levels=levels, cmap='turbo')
    wall = 3.0 * np.tan(np.deg2rad(half_angle))
    ax.add_patch(Polygon([(0, 0), (3, wall), (3, -wall)], facecolor='#263e50', edgecolor='#f3f8fa', lw=1.2, zorder=5))
    ax.text(1.55, -.035, 'CONE', ha='center', fontsize=13, weight='bold', color=TEXT, zorder=6)
    ax.set(xlim=(-.05, 3.01), ylim=(-.83, .83), aspect='equal')
    style(ax, 'Axial distance / display length', 'Radial distance / display length')
    ax.set_title('Self-similar pressure reconstruction', loc='left', pad=15, fontsize=13, color=TEXT)
    colorbar(fig, cf, ax, r'Stored pressure ratio  $p/p_\infty$')
    side = fig.add_axes([.75, .29, .20, .33])
    side.plot(a[:, 0], a[:, 6], color=ACCENT, lw=2.5)
    side.fill_between(a[:, 0], 1, a[:, 6], color=ACCENT, alpha=.12)
    side.axhline(1, color=MUTED, lw=.8, ls='--')
    side.set_title('Original 1D profile', loc='left', pad=15, color=TEXT, fontsize=13)
    style(side, 'Angle above cone surface (deg)', r'$p/p_\infty$')
    side.grid(alpha=.12)
    fig.text(.75, .21, 'Warm colors: higher pressure\nCool colors: approach to farfield', color=MUTED, fontsize=10, linespacing=1.8)
    save(fig, 'conical-flow', 'HISTORICAL DATA  |  Angular profile extended along rays; mirrored for display. Not a newly validated CFD solution.')
    evidence['cases']['conical-flow'] = {'inputs': [fingerprint(source)], 'samples': len(a),
        'pressure_ratio_min': float(a[:, 6].min()), 'pressure_ratio_max': float(a[:, 6].max()),
        'angular_offset_degrees': [float(a[0, 0]), float(a[-1, 0])],
        'display_cone_half_angle_degrees': half_angle,
        'transformation': 'Angular pressure profile reconstructed along rays and mirrored about axis. Radial scale is arbitrary.',
        'configuration_status': 'Original run configuration and convergence not independently established. Filename does not prove flux choice.'}

def airfoil():
    source = 'results/02-airfoil-85-medium/output.plt'
    lines = (ROOT / source).read_text().splitlines()
    names = re.findall(r'"([^"]+)"', lines[0])
    assert names == ['X', 'Y', 'U', 'V', 'P', 'MACH', 'CPP', 'S']
    n = int(re.search(r'N=\s*(\d+)', lines[1]).group(1))
    e = int(re.search(r'E=\s*(\d+)', lines[1]).group(1))
    raw = np.fromstring(' '.join(lines[2:]).replace('D', 'E'), sep=' ')
    assert raw.size == n * len(names) + e * 3
    a = raw[:n*8].reshape(n, 8)
    connectivity = raw[n*8:].reshape(e, 3)
    assert np.isfinite(a).all() and np.equal(connectivity, np.round(connectivity)).all()
    tri = connectivity.astype(int) - 1
    assert tri.min() >= 0 and tri.max() < n
    xy = a[:, :2]
    edge1 = xy[tri[:, 1]] - xy[tri[:, 0]]
    edge2 = xy[tri[:, 2]] - xy[tri[:, 0]]
    areas = (edge1[:, 0] * edge2[:, 1] - edge1[:, 1] * edge2[:, 0]) / 2
    assert (np.abs(areas) > 1e-14).all()
    mesh = mtri.Triangulation(a[:, 0], a[:, 1], triangles=tri)
    edges = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    boundary = unique[counts == 1]
    centers = xy[boundary].mean(axis=1)
    body_edges = boundary[(centers[:, 0] >= -1e-7) & (centers[:, 0] <= 1.000001) & (np.abs(centers[:, 1]) < .2)]
    body = np.unique(body_edges)
    assert len(body) > 10
    fig = canvas(2, 'Compressible airfoil', 'A triangular-mesh view of acceleration and compression around the archived airfoil.',
        [('mesh nodes', f'{n:,}'), ('triangular cells', f'{e:,}'), ('source Mach setting', '0.85')])
    ax = fig.add_axes([.07, .20, .56, .48])
    cf = ax.tricontourf(mesh, a[:, 5], levels=np.linspace(a[:, 5].min(), a[:, 5].max(), 90), cmap='turbo')
    ax.tricontour(mesh, a[:, 5], levels=[1.0], colors=['#f9f9f9'], linewidths=.9)
    for edge in body_edges:
        ax.plot(xy[edge, 0], xy[edge, 1], color=TEXT, lw=1.0)
    ax.set(xlim=(-.25, 1.55), ylim=(-.55, .55), aspect='equal')
    style(ax, 'x / chord', 'y / chord')
    ax.set_title('Stored Mach field  |  white isoline: M = 1', loc='left', pad=15, fontsize=13, color=TEXT)
    colorbar(fig, cf, ax, 'Stored MACH (dimensionless)')
    side = fig.add_axes([.75, .29, .20, .33])
    for sign, name, col in [(1, 'Upper surface', ACCENT), (-1, 'Lower surface', '#ffbc7a')]:
        ids = body[sign * a[body, 1] >= -1e-9]
        ids = ids[np.argsort(a[ids, 0])]
        side.plot(a[ids, 0], a[ids, 6], color=col, lw=2, marker='.', markersize=3, label=name)
    side.invert_yaxis()
    side.set_title('Surface pressure coefficient', loc='left', pad=15, color=TEXT, fontsize=13)
    style(side, 'x / chord', 'Stored CPP')
    side.legend(frameon=False, fontsize=9, labelcolor=TEXT)
    side.grid(alpha=.12)
    fig.text(.75, .195, 'Warm colors: higher local Mach\nSurface curves: archived CPP values', color=MUTED, fontsize=10, linespacing=1.8)
    save(fig, 'compressible-airfoil', 'HISTORICAL DATA  |  Original triangular connectivity retained. Archived fields are not a newly validated solution.')
    evidence['cases']['compressible-airfoil'] = {'inputs': [fingerprint(source)], 'nodes': n, 'triangles': e,
        'mach_min': float(a[:, 5].min()), 'mach_max': float(a[:, 5].max()), 'surface_nodes': len(body),
        'transformation': 'Piecewise-linear contours on original triangles; surface CPP from topological boundary nodes.',
        'configuration_status': '0.85 is the preserved source setting; original output convergence not independently established.'}

def viscous():
    base = 'projects/03-viscous-airfoil/roohi-code/Viscous/'
    sources = [base + name for name in ['node.dat', 'connectivity1.DAT', 'connectivity2.dat']]
    raw_nodes = np.loadtxt(ROOT / sources[0])
    assert np.equal(raw_nodes[:, 0], np.arange(1, len(raw_nodes) + 1)).all()
    nodes = raw_nodes[:, 1:3]
    cells = []
    for path, count in zip(sources[1:], [4, 3]):
        raw = np.loadtxt(ROOT / path, dtype=int)
        assert (raw[:, 2] == count).all()
        cells.extend(row[3:3+count] - 1 for row in raw)
    assert len(nodes) == 838 and len(cells) == 1374
    polygons = [nodes[cell] for cell in cells]
    areas = np.array([abs(np.dot(p[:, 0], np.roll(p[:, 1], -1)) - np.dot(p[:, 1], np.roll(p[:, 0], -1))) / 2 for p in polygons])
    assert np.isfinite(areas).all() and (areas > 1e-14).all()
    norm = LogNorm(vmin=areas.min(), vmax=areas.max())
    fig = canvas(3, 'Viscous airfoil', 'Near-wall resolution and the transition from quadrilateral layers to triangular cells.',
        [('source Reynolds setting', '1,000'), ('quadrilaterals + triangles', '232 + 1,142'), ('source Mach setting', '1.2')])
    ax = fig.add_axes([.07, .20, .56, .48])
    collection = PolyCollection(polygons, array=areas, cmap='viridis', norm=norm, edgecolors='#aac4ce', linewidths=.12)
    ax.add_collection(collection)
    ax.set(xlim=(-.3, 1.6), ylim=(-.57, .57), aspect='equal')
    style(ax, 'x / chord', 'y / chord')
    ax.set_title('Actual mixed mesh  |  color represents cell area', loc='left', pad=15, fontsize=13, color=TEXT)
    colorbar(fig, collection, ax, r'Cell area / chord$^2$  (log scale)')
    side = fig.add_axes([.75, .29, .20, .33])
    zoom = PolyCollection(polygons, array=areas, cmap='viridis', norm=norm, edgecolors=TEXT, linewidths=.45)
    side.add_collection(zoom)
    side.set(xlim=(-.012, .14), ylim=(-.012, .13), aspect='equal')
    side.set_title('Leading-edge mesh detail', loc='left', pad=15, color=TEXT, fontsize=13)
    style(side, 'x / chord', 'y / chord')
    fig.text(.75, .21, 'Dark colors: smaller cells\nBright colors: larger cells', color=MUTED, fontsize=10, linespacing=1.8)
    save(fig, 'viscous-airfoil', 'MESH DATA  |  Colors show geometric cell area, not velocity, pressure or a computed viscous-flow field.')
    evidence['cases']['viscous-airfoil'] = {'inputs': [fingerprint(p) for p in sources],
        'nodes': len(nodes), 'cells': len(cells), 'quadrilaterals': 232, 'triangles': 1142,
        'area_min_chord_squared': float(areas.min()), 'area_max_chord_squared': float(areas.max()),
        'transformation': 'Shoelace area of original polygons; no triangulation or flow-field interpolation.',
        'configuration_status': 'Reynolds 1000 and Mach 1.2 are source settings; this visualization contains no flow solution.'}

if __name__ == '__main__':
    cone(); airfoil(); viscous()
    (OUT / 'figure-provenance.json').write_text(json.dumps(evidence, indent=2) + '\n', encoding='utf8')
    print(json.dumps({'figures_rendered': 3, 'cases': list(evidence['cases']), 'inputs_finite': True, 'geometry_checks_passed': True}))
