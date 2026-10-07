"""Audit new FV components independently of nonlinear cone convergence."""
from pathlib import Path
import json,numpy as np
from flux import FLUXES,numerical_flux
from state import conserved,physical_flux
from additional_fluxes import riemann_sample,jst_flux
from advanced_reconstruction import RECONSTRUCTIONS,matrices
from reconstruction import coefficients,evaluate
from shock_tube import exact
ROOT=Path(__file__).resolve().parent

def main():
    rng=np.random.default_rng(17)
    q=np.c_[rng.uniform(.2,2,20),rng.uniform(-.4,.4,20),rng.uniform(-1,1,20),rng.uniform(.1,2,20)]
    u=conserved(q);equal={f:float(np.max(abs(numerical_flux(f,u,u)-physical_flux(u)))) for f in FLUXES}
    assert max(equal.values())<3e-12
    ref,ps,us=exact(np.linspace(0,1,400))
    got=riemann_sample([1,0,0,1],[.125,0,0,.1],(np.linspace(0,1,400)-.5)/.2)
    difference=float(np.max(abs(got-ref)));assert difference<1e-12
    vacuum=riemann_sample([1,0,-6,.4],[1,0,6,.4],0.)
    assert vacuum[0]==0 and vacuum[3]==0
    smooth={}
    for kind in RECONSTRUCTIONS:
        errors=[]
        for n in [20,40,80,160]:
            h=1/n;x=(np.arange(n+8)-3.5)*h
            density=1+.2*np.sin(2*np.pi*x)*np.sinc(h)
            avg=conserved(np.c_[density,np.zeros(n+8),np.ones(n+8),np.ones(n+8)])
            pol=coefficients(avg,kind,h)
            face=evaluate(pol,.5)[1:-1,0]
            # Upwind positive unit speed: the face polynomial from the left cell.
            derivative=np.diff(face[:-1])/h
            xc=(np.arange(n)+.5)*h
            exact_derivative=.2*2*np.pi*np.cos(2*np.pi*xc)*np.sinc(h)
            errors.append(float(np.mean(abs(derivative-exact_derivative))))
        smooth[kind]={'cells':[20,40,80,160],'density_rhs_L1':errors,
                      'observed_orders':np.log2(np.array(errors[:-1])/errors[1:]).tolist()}
    record={'passed':True,'equal_state_max_errors':equal,'independent_Sod_max_error':difference,
            'vacuum_exact_sample_passed':True,'smooth_FV_operator':smooth,
            'note':'Observed smooth orders are measurements; discontinuities and cone source/wall closure reduce problem order.'}
    (ROOT/'additional-verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'passed':True,'fluxes':len(FLUXES),'reconstructions':len(RECONSTRUCTIONS),
                      'fine_grid_orders':{k:round(v['observed_orders'][-1],3) for k,v in smooth.items()}},indent=2))
if __name__=='__main__':main()
