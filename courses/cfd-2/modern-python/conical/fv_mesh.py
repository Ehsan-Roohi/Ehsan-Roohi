"""Nonuniform cell-average CWENO reconstruction and face diffusion operators.

All polynomial fits integrate over control volumes; inputs are cell averages.
The diffusion fit is independently selectable as linear (2) or cubic (4).
"""
import numpy as np
from reconstruction import smoothness_matrix


class Mesh:
    def __init__(self, lower, upper, cells, stretch=0):
        s = np.linspace(0, 1, cells+1)
        z = s if stretch == 0 else np.expm1(stretch*s)/np.expm1(stretch)
        self.faces = lower+(upper-lower)*z
        widths = np.diff(self.faces)
        self.ext_faces = np.r_[lower-np.cumsum(widths[:3])[::-1], self.faces,
                               upper+np.arange(1, 4)*widths[-1]]
        self.width = np.diff(self.ext_faces)
        self.center = .5*(self.ext_faces[1:]+self.ext_faces[:-1])
        self.n = cells
        self.fits = {}
        for offsets in [(-1, 0, 1), (-1, 0), (0, 1), (-2, -1, 0, 1, 2),
                        (-2, -1, 0), (0, 1, 2)]:
            indices = np.arange(2, cells+4)[:, None]+np.array(offsets)
            lo = (self.ext_faces[indices]-self.center[2:-2, None])/self.width[2:-2, None]
            hi = (self.ext_faces[indices+1]-self.center[2:-2, None])/self.width[2:-2, None]
            degree = len(offsets)-1
            mat = np.stack([(hi**(k+1)-lo**(k+1))/((k+1)*(hi-lo)) for k in range(degree+1)], -1)
            self.fits[offsets] = (indices, np.linalg.inv(mat))
        self.face_fits = {}
        # Face f separates cells f-1 and f in the extended mesh.
        f = np.arange(3, cells+4)
        scale = .5*(self.width[f-1]+self.width[f])
        for order in (2, 4):
            offsets = np.arange(-order//2, order//2)
            indices = f[:, None]+offsets
            lo = (self.ext_faces[indices]-self.ext_faces[f, None])/scale[:, None]
            hi = (self.ext_faces[indices+1]-self.ext_faces[f, None])/scale[:, None]
            mat = np.stack([(hi**(k+1)-lo**(k+1))/((k+1)*(hi-lo)) for k in range(order)], -1)
            inv = np.linalg.inv(mat)
            self.face_fits[order] = (indices, inv[:, 0], inv[:, 1]/scale[:, None])

    def fit(self, ext, offsets):
        indices, inv = self.fits[offsets]
        return np.einsum('nij,njk->nik', inv, ext[indices])

    def reconstruct(self, ext, kind='cweno3', gamma=1.4):
        center = ext[2:-2]
        if kind == 'primitive-minmod':
            from state import primitive
            qp = self.reconstruct(primitive(ext, gamma), 'minmod')
            a, b = qp[:, 0], qp[:, 1]
            rho, vr, vt, p = np.moveaxis(a, -1, 0)
            dr, dvr, dvt, dp = np.moveaxis(b, -1, 0)
            result = np.zeros((len(center), 4, 4))
            result[:, 0] = center
            result[:, 1, 0] = dr
            for component, velocity, slope in [(1, vr, dvr), (2, vt, dvt)]:
                result[:, 1, component] = rho*slope+dr*velocity
                result[:, 2, component] = dr*slope
            v2, vdv, dv2 = vr*vr+vt*vt, 2*(vr*dvr+vt*dvt), dvr*dvr+dvt*dvt
            result[:, 1, 3] = dp/(gamma-1)+.5*(rho*vdv+dr*v2)
            result[:, 2, 3] = .5*(rho*dv2+dr*vdv)
            result[:, 3, 3] = .5*dr*dv2
            # Exact mean correction after primitive-to-conservative polynomial
            # conversion. No loss of conservative FV cell averages.
            result[:, 0] -= result[:, 2]/12
            return result
        if kind == 'first':
            return center[:, None].copy()
        if kind in ('muscl', 'minmod'):
            dl = (center-ext[1:-3])/(self.center[2:-2]-self.center[1:-3])[:, None]
            dr = (ext[3:-1]-center)/(self.center[3:-1]-self.center[2:-2])[:, None]
            slope = np.where(dl*dr > 0, np.sign(dl)*np.minimum(abs(dl), abs(dr)), 0)
            if kind == 'muscl':
                dc = (ext[3:-1]-ext[1:-3])/(self.center[3:-1]-self.center[1:-3])[:, None]
                slope = np.where(dl*dr > 0, np.sign(dc)*np.minimum(np.minimum(2*abs(dl), 2*abs(dr)), abs(dc)), 0)
            return np.stack((center, slope*self.width[2:-2, None]), axis=1)
        order = int(kind[5])
        r = (order-1)//2
        high = self.fit(ext, tuple(range(-r, r+1)))
        d = np.array([.5, .25, .25]) if order == 3 else np.array([.5, .125, .25, .125])
        candidates = [np.zeros_like(high)]
        for j in range(r+1):
            low = np.zeros_like(high)
            low[:, :r+1] = self.fit(ext, tuple(range(-r+j, j+1)))
            candidates.append(low)
        candidates[0] = (high-sum(w*p for w, p in zip(d[1:], candidates[1:])))/d[0]
        candidates = np.stack(candidates, axis=1)
        if kind.endswith('-char'):
            from characteristics import basis
            right, left = basis(center, gamma)
            candidates = np.einsum('nij,nspj->nspi', left, candidates)
            scale_center = np.einsum('nij,nj->ni', left, center)
        else:
            scale_center = center
        beta = np.maximum(np.einsum('nsik,ij,nsjk->nsk', candidates, smoothness_matrix(2*r), candidates), 0)
        scale = np.maximum(abs(scale_center), 1e-3)**2
        eps = (self.width[2:-2, None]**2+1e-14)*scale
        alpha = d[None, :, None]/(beta+eps[:, None, :])**2
        result = np.einsum('nsk,nsik->nik', alpha/alpha.sum(axis=1, keepdims=True), candidates)
        if kind.endswith('-char'):
            result = np.einsum('nij,npj->npi', right, result)
        return result

    def face_value_gradient(self, averages, order=2):
        indices, value, gradient = self.face_fits[order]
        data = averages[indices]
        return (np.einsum('ni,nik->nk', value, data),
                np.einsum('ni,nik->nk', gradient, data))
