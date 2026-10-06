"""Reproduce isolated behaviors of the pinned upstream code, not a full CFD run."""
from pathlib import Path
from types import SimpleNamespace
import argparse, ast, importlib.util, json, tempfile, urllib.request
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import os
from contextlib import nullcontext

COMMIT='6a1d90684c9462aa48416894f147910ca20e238f'
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--source-dir');parser.add_argument('--output');args=parser.parse_args()
 source=Path(args.source_dir) if args.source_dir else Path('.review-cache')
 source.mkdir(exist_ok=True)
 if not args.source_dir:
  for name,remote in [('interpolations.py','src/numerics/interpolations.py'),('helpers.py','src/helpers.py')]:
   url=f'https://raw.githubusercontent.com/jyoun35/SIMPLE/{COMMIT}/{remote}'
   source.joinpath(name).write_bytes(urllib.request.urlopen(url,timeout=30).read())
 spec=importlib.util.spec_from_file_location('upstream_interpolations',source/'interpolations.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 def face(d,pn):
  mesh=SimpleNamespace(cell_centroids=np.array([[0.,0.],d]))
  return m.rhie_and_chow_interp(mesh,0,1,1.,0.,1.,1.,0.,0.,0.,0.,0.,pn,1.,1.,1.,1.,1.,1.)
 orth=face([1.,0.],1.);skew=face([1.,1.],2.)
 assert np.allclose(orth,[0.,0.])
 assert np.isclose(skew[0],-(2.-np.sqrt(2.)))
 result={'upstream_commit':COMMIT,'numpy_version':np.__version__,'matplotlib_version':matplotlib.__version__,'scope':'Isolated source-function diagnostics; no full flow simulation','orthogonal_affine_pressure_velocity':list(orth),'skew_affine_pressure_velocity':list(skew),'skew_expected_velocity_for_affine_pressure':[0.,0.]}
 # Extract the two unmodified function bodies; avoid importing SciPy-dependent helpers.
 tree=ast.parse((source/'helpers.py').read_text('utf-8-sig'))
 chosen=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['save_checkpoint','plot_pressure_distribution']]
 env=dict(os=os,np=np,plt=plt,PolyCollection=PolyCollection)
 exec(compile(ast.Module(body=chosen,type_ignores=[]),'upstream_helpers.py','exec'),env)
 mesh=SimpleNamespace(node_coords=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]]),cells=[[0,1,2,3]],boundary_faces=[(0,1,0,'slip-wall'),(2,3,0,'slip-wall')])
 soln=SimpleNamespace(u=np.array([1.]),v=np.array([0.]),p=np.array([101325.]))
 bcs={'pressure-outlet':{'p':101325.},'velocity-inlet':{'u':1.,'v':0.}}
 captures=[];original_close=plt.close
 def capture(fig=None):
  if fig is not None and hasattr(fig,'axes'):
   for ax in fig.axes:
    if ax.get_title().startswith('Pressure Relative'):
     captures.extend([tuple(c.get_clim()) for c in ax.collections])
  original_close(fig)
 plt.close=capture
 dest=source/'checkpoint-diagnostic';dest.mkdir(exist_ok=True)
 with nullcontext(str(dest)) as dest:
  try:
   env['save_checkpoint'](dest,1,mesh,soln,np.array([0.1]),np.array([[0.,.1,.1,.1]]),bcs,{'x':np.array([0.,1.]),'Cp':np.array([0.,0.])})
   result['checkpoint_exception']=None
  except Exception as exc:
   result['checkpoint_exception']={'type':type(exc).__name__,'message':str(exc)}
 plt.close=original_close;plt.close('all')
 result['pressure_plot_declared_clim']=[500.,-500.];result['pressure_plot_observed_clim']=captures
 print(json.dumps(result,indent=2))
 if args.output:Path(args.output).write_text(json.dumps(result,indent=2)+'\n','utf8')
if __name__=='__main__':main()

