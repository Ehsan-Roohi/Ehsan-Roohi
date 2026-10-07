"""Cell-average polynomial FV reconstruction, including CWENO3 and CWENO5.

One nonlinear polynomial per cell supplies both faces and source quadrature.
Independent implementation of the CWENO framework; no upstream code copied.
"""
from functools import lru_cache
from math import factorial
import numpy as np


def average_matrix(offsets, degree):
    return np.array([[( (k+.5)**(m+1)-(k-.5)**(m+1))/(m+1)
                      for m in range(degree+1)] for k in offsets])


@lru_cache(None)
def smoothness_matrix(degree):
    b = np.zeros((degree+1, degree+1))
    for derivative in range(1, degree+1):
        for i in range(derivative, degree+1):
            for j in range(derivative, degree+1):
                exponent = i+j-2*derivative
                integral = (.5**(exponent+1)-(-.5)**(exponent+1))/(exponent+1)
                b[i, j] += factorial(i)/factorial(i-derivative)*factorial(j)/factorial(j-derivative)*integral
    return b


@lru_cache(None)
def operators(order):
    r = (order-1)//2
    opt = np.linalg.inv(average_matrix(range(-r, r+1), 2*r))
    stencils = [tuple(range(-r+j, j+1)) for j in range(r+1)]
    low = [np.linalg.inv(average_matrix(s, r)) for s in stencils]
    weights = np.array([.5, .25, .25]) if order == 3 else np.array([.5, .125, .25, .125])
    return opt, stencils, low, weights


def coefficients(ext, kind, h, gamma=1.4):
    """Return polynomials for ext[2:-2]; physical cells have three ghosts."""
    center = ext[2:-2]
    if kind == "first":
        return center[:, None, :].copy()
    if kind == "muscl":
        dl, dr = center-ext[1:-3], ext[3:-1]-center
        dc = .5*(dl+dr)
        slope = np.where((dl*dr) > 0, np.sign(dc)*np.minimum(np.minimum(2*abs(dl), 2*abs(dr)), abs(dc)), 0)
        return np.stack((center, slope), axis=1)
    order = int(kind[5])
    opt, stencils, low, d = operators(order)
    r = (order-1)//2
    data = np.stack([ext[2+k:len(ext)-2+k] for k in range(-r, r+1)], axis=1)
    high = np.einsum("ij,njk->nik", opt, data)
    candidates = [np.zeros_like(high)]
    for offsets, matrix in zip(stencils, low):
        vals = np.stack([ext[2+k:len(ext)-2+k] for k in offsets], axis=1)
        c = np.zeros_like(high)
        c[:, :r+1] = np.einsum("ij,njk->nik", matrix, vals)
        candidates.append(c)
    candidates[0] = (high-sum(w*c for w, c in zip(d[1:], candidates[1:])))/d[0]
    candidates = np.stack(candidates, axis=1)
    if kind.endswith('-char'):
        from characteristics import basis
        right, left = basis(center, gamma)
        candidates = np.einsum('nij,nspj->nspi', left, candidates)
        data = np.einsum('nij,npj->npi', left, data)
    beta = np.maximum(np.einsum("nsik,ij,nsjk->nsk", candidates, smoothness_matrix(2*r), candidates), 0)
    # epsilon proportional to h^2; component-scaled to avoid unit imbalance.
    scale = np.maximum(np.max(abs(data), axis=1), 1e-3)**2
    epsilon = (h*h+1e-14)*scale
    alpha = d[None, :, None]/(beta+epsilon[:, None, :])**2
    weights = alpha/np.sum(alpha, axis=1, keepdims=True)
    result = np.einsum("nsk,nsik->nik", weights, candidates)
    if kind.endswith('-char'):
        result = np.einsum('nij,npj->npi', right, result)
    return result


def evaluate(c, s):
    if np.ndim(s) == 0:
        return np.einsum("nik,i->nk", c, np.array([s**i for i in range(c.shape[1])]))
    return np.einsum("nik,qi->nqk", c, np.stack([np.asarray(s)**i for i in range(c.shape[1])], axis=1))


def limit_admissibility(c, means, gamma, sample_points):
    """Scale each polynomial toward its admissible mean, preserving cell average.

    This enforces the sampled reconstruction states, not a global positivity theorem.
    """
    from state import primitive
    values = evaluate(c, sample_points)
    bad = np.any((values[..., 0] <= 1e-12) | (primitive(values, gamma)[..., 3] <= 1e-12), axis=1)
    if not np.any(bad):
        return c, 0
    delta = c[bad].copy()
    delta[:, 0] -= means[bad]
    lo, hi = np.zeros(sum(bad)), np.ones(sum(bad))
    mean = means[bad]
    for _ in range(35):
        mid = .5*(lo+hi)
        v = mean[:, None, :]+mid[:, None, None]*(values[bad]-mean[:, None, :])
        good = np.all((v[..., 0] > 1e-12) & (primitive(v, gamma)[..., 3] > 1e-12), axis=1)
        lo, hi = np.where(good, mid, lo), np.where(good, hi, mid)
    c[bad] = delta*(.999*lo)[:, None, None]
    c[bad, 0] += mean
    return c, int(sum(bad))
