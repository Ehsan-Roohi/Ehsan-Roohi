"""Check the archived co-located airfoil meshes and compiled source dimensions."""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def numbers(file):
    return [[float(x.replace('D', 'E').replace('d', 'e')) for x in line.split()]
            for line in file.read_text('latin1').splitlines() if line.strip()]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', default='docs/validation/mesh-validation.json')
    a = p.parse_args()
    report = {'scope': 'Geometry/connectivity and source-dimension checks; not flow validation',
              'cases': []}
    for node_file in sorted((ROOT / 'projects').rglob('*')):
        if node_file.name.lower() != 'node.dat':
            continue
        directory = node_file.parent
        input_files = sorted(f for f in directory.iterdir()
                             if f.name.lower().startswith('connectivity') and f.suffix.lower() == '.dat')
        if not input_files:
            continue
        nodes = {int(row[0]): (row[1], row[2]) for row in numbers(node_file)}
        polygons = []
        issues = []
        for file in input_files:
            for row in numbers(file):
                ids = [int(x) for x in row[3:]]
                if len(ids) not in (3, 4):
                    issues.append(f'Unsupported connectivity record in {file.name}')
                else:
                    polygons.append(ids)
        areas, edges, quad_y_errors = [], [], []
        for ids in polygons:
            if any(i not in nodes for i in ids):
                issues.append('Connectivity references an absent node')
                continue
            if len(set(ids)) != len(ids):
                issues.append('Repeated node in a cell')
            xy = [nodes[i] for i in ids]
            signed_area = 0.5 * sum(xy[k][0]*xy[(k+1) % len(xy)][1]
                                    - xy[(k+1) % len(xy)][0]*xy[k][1] for k in range(len(xy)))
            areas.append(signed_area)
            edges.extend(tuple(sorted((ids[k], ids[(k+1) % len(ids)]))) for k in range(len(ids)))
            if len(ids) == 4:
                correct = sum(pt[1] for pt in xy) / 4
                source_formula = (xy[0][1] + xy[1][1] + xy[2][1] + xy[3][0]) / 4
                quad_y_errors.append(abs(source_formula-correct))
        multiplicity = Counter(edges)
        counts = Counter(len(ids) for ids in polygons)
        entry = {'case': directory.relative_to(ROOT).as_posix(), 'nodes': len(nodes),
                 'cells': len(polygons), 'triangle_cells': counts[3], 'quad_cells': counts[4],
                 'positive_signed_areas': sum(x > 0 for x in areas),
                 'negative_signed_areas': sum(x < 0 for x in areas),
                 'zero_area_cells': sum(abs(x) < 1e-14 for x in areas),
                 'minimum_absolute_area': min(map(abs, areas)),
                 'boundary_edges': sum(v == 1 for v in multiplicity.values()),
                 'non_manifold_edges': sum(v > 2 for v in multiplicity.values()),
                 'input_hashes': {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                                  for f in [node_file, *input_files]}, 'sources': [],
                 'issues': sorted(set(issues))}
        for source in sorted(f for f in directory.iterdir() if f.suffix.lower() in ('.for', '.f90')):
            text = source.read_text('latin1')
            nn = re.search(r'(?i)parameter\s*\(\s*NN\s*=\s*(\d+)', text)
            ne = re.search(r'(?i)parameter\s*\(\s*NE\s*=\s*(\d+)', text)
            row = {'file': source.name, 'NN': int(nn[1]) if nn else None,
                   'NE': int(ne[1]) if ne else None}
            row['dimensions_match_inputs'] = row['NN'] == len(nodes) and row['NE'] == len(polygons)
            bad_centroid = bool(re.search(r'(?i)Y_CELL\(I\)\s*=.*YNODE\(N3\)\+XNODE\(N4\)', text))
            row['quad_Y_centroid_uses_X4'] = bad_centroid
            if bad_centroid and quad_y_errors:
                row['maximum_quad_Y_centroid_error'] = max(quad_y_errors)
                row['affected_quad_cells'] = sum(x > 1e-14 for x in quad_y_errors)
            entry['sources'].append(row)
        entry['geometry_checks_passed'] = not (entry['issues'] or entry['zero_area_cells']
                                                or entry['non_manifold_edges'])
        report['cases'].append(entry)
    path = ROOT / a.output
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+'\n', 'utf8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
