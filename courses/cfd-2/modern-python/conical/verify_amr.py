"""Independent nonuniform moments, remapping, interface balance and contact tests."""
from pathlib import Path
import json
import numpy as np
from shock_amr import AdaptiveTube,AMRConfig,LeafMesh,SUPPORTED
from shock_suite import TubeSolver,TubeConfig
from state import conserved,admissible
from flux import FLUXES

ROOT=Path(__file__).resolve().parent


def main():
    solver=AdaptiveTube(AMRConfig(base_cells=16,max_level=3))
    # Nonuniform polynomial input contains exact FV means, including ghosts.
    m=solver.mesh; a=m.ext_faces[:-1];b=m.ext_faces[1:]
    averages=np.stack([np.ones_like(a),.5*(a+b),(a*a+a*b+b*b)/3],axis=-1)
    fitted=m.fit(averages,(-1,0,1))
    x=m.center[2:-2];h=m.width[2:-2]
    expected=np.zeros_like(fitted)
    expected[:,0]=np.stack([np.ones_like(x),x,x*x],-1)
    expected[:,1]=np.stack([0*x,h,2*x*h],-1)
    expected[:,2,2]=h*h
    moment_error=float(np.max(abs(fitted-expected)))
    assert moment_error<1e-12
    # Admissible linear conservative field; split only away from outflow ghosts.
    solver.start=np.arange(16)*8;solver.span=np.full(16,8);solver.update_mesh()
    solver.u=conserved(np.column_stack([1+.02*solver.x,0*solver.x,.1+0*solver.x,1+0*solver.x]))
    old=solver.u.copy();integral=solver.integral().copy()
    split=np.zeros(16,bool);split[4:12]=True;solver.prolong(split)
    split_error=float(np.max(abs(solver.integral()-integral)))
    assert split_error<1e-13 and admissible(solver.u)
    solver.adapt()
    merge_error=float(np.max(abs(solver.integral()-integral)))
    assert merge_error<1e-13 and len(solver.u)==16
    assert np.max(abs(solver.u-old))<1e-13
    # Shared interface fluxes telescope on a nonuniform partition for all fluxes.
    balances={};constant_errors={}
    for flux in FLUXES:
        c=AMRConfig(flux=flux,base_cells=16,max_level=3)
        s=AdaptiveTube(c);u=conserved(np.tile([1.,.05,.2,1.],(len(s.u),1)))
        rhs,budget=s.rhs(u);constant_errors[flux]=float(np.max(abs(rhs)))
        balances[flux]=float(np.max(abs(np.einsum('n,nv->v',s.dx,rhs)-budget)))
        assert constant_errors[flux]<1e-10 and balances[flux]<1e-12
        rhs,budget=s.rhs(s.u)
        balances[flux]=float(np.max(abs(np.einsum('n,nv->v',s.dx,rhs)-budget)))
        assert balances[flux]<1e-12
    # First-order uniform path must reproduce the existing shared FV operator.
    a=AdaptiveTube(AMRConfig(base_cells=40,max_level=0,reconstruction='first'))
    b=TubeSolver(TubeConfig(cells=40,reconstruction='first'))
    uniform_error=float(np.max(abs(a.rhs(a.u)[0]-b.rhs(b.u)[0])))
    assert uniform_error<1e-10
    stationary=AdaptiveTube(AMRConfig(case='stationary-contact',base_cells=40,max_level=3))
    r=stationary.run();assert r['integrity_passed']
    assert r['errors']['density']['L1']<1e-12
    assert np.allclose(np.diff(r['faces']),1/(40*8)*stationary.span)
    ratio=stationary.span[1:]/stationary.span[:-1]
    assert np.max(ratio)<=2 and np.min(ratio)>=.5
    assert stationary.start[0]==0 and stationary.start[-1]+stationary.span[-1]==stationary.finest
    pilots={}
    for reconstruction in SUPPORTED:
        p=AdaptiveTube(AMRConfig(base_cells=40,max_level=2,reconstruction=reconstruction)).run()
        assert p['integrity_passed']
        pilots[reconstruction]={'integrity_passed':True,'density_L1':p['errors']['density']['L1'],
                                'conservation_defect':p['conservation_budget_defect']}
    result={'passed':True,'nonuniform_polynomial_moment_error':moment_error,
            'split_conservation_error':split_error,'merge_conservation_error':merge_error,
            'uniform_first_order_rhs_error':uniform_error,'constant_rhs_max_errors':constant_errors,
            'nonuniform_flux_balance_errors':balances,'stationary_contact_density_L1':r['errors']['density']['L1'],
            'stationary_contact_conservation_defect':r['conservation_budget_defect'],
            'dyadic_2_to_1_partition_verified':True,'supported_reconstruction_Sod_pilots':pilots}
    (ROOT/'amr-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
