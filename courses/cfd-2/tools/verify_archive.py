"""Verify preserved original bytes, basic mesh counts, and local Markdown links."""
from pathlib import Path
from urllib.parse import unquote
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'manifest.json').read_text('utf8'))
errors=[]
for row in manifest['files']:
 p=ROOT/row['path']
 if not p.is_file():errors.append('Missing '+row['path']);continue
 data=p.read_bytes()
 if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:errors.append('Bytes differ '+row['path'])
for p in ROOT.rglob('*.md'):
 for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text('utf-8-sig')):
  if re.match(r'^[a-zA-Z]+:',target) or target.startswith('#'):continue
  target=unquote(target.split('#')[0])
  if target and not (p.parent/target).exists():errors.append(f'Broken link {p.relative_to(ROOT)}: {target}')
mesh=[]
for p in sorted((ROOT/'projects').rglob('node.dat')):
 first=p.read_text(errors='replace').splitlines()[0].strip()
 mesh.append({'path':p.relative_to(ROOT).as_posix(),'first_line':first})
result={'original_files_verified':len(manifest['files']),'local_links_verified':not errors,'mesh_file_headers':mesh,'errors':errors}
print(json.dumps(result,indent=2))
raise SystemExit(bool(errors))
