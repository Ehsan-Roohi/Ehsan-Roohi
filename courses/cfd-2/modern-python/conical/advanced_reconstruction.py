"""Finite-volume face reconstructions from cell averages, independently derived.

WENO-JS/Z (3/5), WENO-JS7, targeted TENO5; conservative ENO2 and
MUSCL-THINC-BVD. Source quadrature uses a separate cell-average CWENO
polynomial: order 3 for WENO3, order 5 for WENO5/7 and TENO5.
Thus WENO7 is seventh-order at faces; the cone source limits overall order.
"""
from functools import lru_cache
import numpy as np
from reconstruction import average_matrix, smoothness_matrix

FACE_METHODS = ['eno2','weno3-js','weno3-z','weno5-js','weno5-z','weno7-js','teno3','teno5','teno7','muscl-thinc-bvd']
LIMITERS = ['muscl-minmod','muscl-vanleer','muscl-superbee','muscl-vanalbada']
RECONSTRUCTIONS = ['first','muscl','cweno3','cweno5','cweno3-char','cweno5-char']+LIMITERS+FACE_METHODS


class FacePolynomial:
    def __init__(self, polynomial, left, right):
        self.polynomial,self.left,self.right = polynomial,left,right
    def __getitem__(self,index):
        return FacePolynomial(self.polynomial[index],self.left[index],self.right[index])


@lru_cache(None)
def matrices(order,side):
    r = (order+1)//2
    offsets = [tuple(range(-r+1+j,j+1)) for j in range(r)]
    full = tuple(range(-r+1,r))
    opt = np.array([side**k for k in range(order)])@np.linalg.inv(average_matrix(full,order-1))
    candidates = []
    for stencil in offsets:
        row = np.array([side**k for k in range(r)])@np.linalg.inv(average_matrix(stencil,r-1))
        big = np.zeros(order)
        for k,v in zip(stencil,row): big[full.index(k)] = v
        candidates.append(big)
    d = np.linalg.lstsq(np.array(candidates).T,opt,rcond=None)[0]
    assert np.max(abs(np.array(candidates).T@d-opt)) < 2e-12
    assert np.min(d)>0
    return offsets,[np.linalg.inv(average_matrix(s,r-1)) for s in offsets],d


def weno_faces(ext,kind,h,gamma):
    from characteristics import basis
    order = int(kind[4])
    r = (order+1)//2
    # Three physical ghosts suffice for WENO7; boundary ghost reconstructions
    # use edge extension. No claim of seventh-order boundary closure is made.
    padded = np.pad(ext,((max(0,r-3),)*2,(0,0)),mode='edge')
    shift = max(0,r-3)
    center = ext[2:-2]
    start,stop = 2+shift,len(padded)-2-shift
    right_basis,left_basis = basis(center,gamma)
    offsets,mats,_ = matrices(order,.5)
    pols=[]
    for stencil,mat in zip(offsets,mats):
        vals=np.stack([padded[start+k:stop+k] for k in stencil],axis=1)
        vals=np.einsum('nij,npj->npi',left_basis,vals)
        pols.append(np.einsum('ij,njk->nik',mat,vals))
    pols=np.stack(pols,axis=1)
    beta=np.maximum(np.einsum('nsik,ij,nsjk->nsk',pols,smoothness_matrix(r-1),pols),0)
    scale=np.maximum(np.max(abs(np.einsum('nij,nj->ni',left_basis,center)),axis=0),1e-3)**2
    eps=1e-12*scale
    if kind.endswith('-z') or kind=='teno5':
        tau=abs(beta[:,0]-beta[:,-1])
    outputs=[]
    for side in [-.5,.5]:
        _,_,d=matrices(order,side)
        if kind.startswith('teno'):
            # TENO5 uses tau5; TENO3/7 use the inverse-beta cutoff detector
            # used by CAELUM. These are documented FV variants, not TENO-A/LAD.
            logg=6*np.log1p(tau[:,None,:]/(beta+eps)) if kind=='teno5' else -6*np.log(beta+eps)
            g=np.exp(logg-np.max(logg,axis=1,keepdims=True))
            keep=g/np.sum(g,axis=1,keepdims=True)>=1e-5
            alpha=d[None,:,None]*keep
        elif kind.endswith('-z'):
            alpha=d[None,:,None]*(1+(tau[:,None,:]/(beta+eps))**2)
        else:
            alpha=d[None,:,None]/(beta+eps)**2
        weight=alpha/np.sum(alpha,axis=1,keepdims=True)
        val=np.einsum('nsik,i->nsk',pols,np.array([side**k for k in range(r)]))
        out=np.einsum('nsk,nsk->nk',weight,val)
        outputs.append(np.einsum('nij,nj->ni',right_basis,out))
    return outputs


def thinc_faces(ext):
    center,prev,nxt=ext[2:-2],ext[1:-3],ext[3:-1]
    low=np.minimum(prev,nxt); span=abs(nxt-prev)
    direction=np.sign(nxt-prev)
    monotone=(center-prev)*(nxt-center)>0
    c=np.clip((center-low)/(span+1e-30),1e-12,1-1e-12)
    beta=1.6
    b=np.exp(beta*direction*(2*c-1))
    a=(b/np.cosh(beta)-1)/np.tanh(beta)
    left=low+span/2*(1+direction*a)
    right=low+span/2*(1+direction*(np.tanh(beta)+a)/(1+a*np.tanh(beta)))
    return np.where(monotone,left,center),np.where(monotone,right,center),monotone


def reconstruct(ext,kind,h,gamma):
    from reconstruction import coefficients,evaluate
    if kind=='eno2':
        center=ext[2:-2];dl=center-ext[1:-3];dr=ext[3:-1]-center
        slope=np.where(abs(dl)<=abs(dr),dl,dr)
        return np.stack((center,slope),axis=1)
    if kind in LIMITERS:
        center=ext[2:-2];a=center-ext[1:-3];b=ext[3:-1]-center
        same=a*b>0
        if kind=='muscl-minmod': slope=np.sign(a)*np.minimum(abs(a),abs(b))
        elif kind=='muscl-vanleer': slope=2*a*b/(a+b+np.where(a+b>=0,1,-1)*1e-30)
        elif kind=='muscl-superbee': slope=np.sign(a)*np.maximum(np.minimum(2*abs(a),abs(b)),np.minimum(abs(a),2*abs(b)))
        else: slope=(a*b*(a+b))/(a*a+b*b+1e-30)
        return np.stack((center,np.where(same,slope,0)),axis=1)
    if kind=='muscl-thinc-bvd':
        polynomial=coefficients(ext,'muscl',h,gamma)
        ml,mr=evaluate(polynomial,-.5),evaluate(polynomial,.5)
        tl,tr,mono=thinc_faces(ext)
        def variation(l,r):
            before=np.r_[r[:1],r[:-1]];after=np.r_[l[1:],l[-1:]]
            return abs(l-before)+abs(r-after)
        use=(variation(tl,tr)<variation(ml,mr))&mono
        # Component-wise one-stage BVD selection; independent educational variant.
        left,right=np.where(use,tl,ml),np.where(use,tr,mr)
        return FacePolynomial(polynomial,left,right)
    left,right=weno_faces(ext,kind,h,gamma)
    source_kind='cweno3' if kind.startswith(('weno3','teno3')) else 'cweno5'
    return FacePolynomial(coefficients(ext,source_kind,h,gamma),left,right)


def limit_faces(obj,means,gamma,samples):
    from reconstruction import limit_admissibility,evaluate
    from state import primitive
    obj.polynomial,count=limit_admissibility(obj.polynomial,means,gamma,samples)
    both=np.stack((obj.left,obj.right),axis=1)
    q=primitive(both,gamma)
    good=np.all((q[...,0]>1e-12)&(q[...,3]>1e-12)&np.all(np.isfinite(q),axis=-1),axis=1)
    if np.all(good): return obj,count
    bad=~good;mean=means[bad];delta=both[bad]-mean[:,None,:]
    lo,hi=np.zeros(len(mean)),np.ones(len(mean))
    for _ in range(36):
        mid=(lo+hi)/2
        trial=mean[:,None,:]+mid[:,None,None]*delta
        ok=np.all((trial[...,0]>1e-12)&(primitive(trial,gamma)[...,3]>1e-12),axis=1)
        lo,hi=np.where(ok,mid,lo),np.where(ok,hi,mid)
    both[bad]=mean[:,None,:]+(.999*lo)[:,None,None]*delta
    obj.left,obj.right=both[:,0],both[:,1]
    return obj,count+int(sum(bad))
