"""Create an isolated legacy-code copy without active Windows deletion calls."""
from pathlib import Path
import argparse,re,shutil,json
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('destination');a=p.parse_args()
 source=(ROOT/a.source).resolve();dest=(ROOT/a.destination).resolve()
 if not source.is_relative_to(ROOT/'projects') or not source.is_dir():p.error('Source must be a directory in this archive projects folder')
 if source==dest or source in dest.parents or dest==ROOT or dest in ROOT.parents:p.error('Destination must be separate from source and archive')
 if dest.exists():p.error('Destination already exists; choose a new run directory')
 dest.parent.mkdir(parents=True,exist_ok=True);shutil.copytree(source,dest)
 changes=[]
 for f in dest.rglob('*'):
  if f.suffix.lower() not in ['.for','.f90']:continue
  data=f.read_bytes();lines=data.splitlines(keepends=True);new=[];removed=[]
  for i,line in enumerate(lines,1):
   stripped=line.lstrip()
   if re.match(rb'(?i)^CALL\s+SYSTEM\s*\(\s*[\x22\x27]del\s',stripped):
    new.append(b'! Cleanup command removed from working copy.'+(b'\r\n' if line.endswith(b'\r\n') else b'\n'));removed.append(i)
   else:new.append(line)
  if removed:f.write_bytes(b''.join(new));changes.append({'file':str(f.relative_to(dest)),'removed_cleanup_lines':removed})
 (dest/'working-copy-changes.json').write_text(json.dumps({'source':a.source,'changes':changes},indent=2)+'\n','utf8')
 print(json.dumps({'destination':str(dest),'changed_sources':len(changes)},indent=2))
if __name__=='__main__':main()
