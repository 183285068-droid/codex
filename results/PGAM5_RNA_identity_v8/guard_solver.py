"""Fixed robust simplex fitting with optional total PGAM5 RNA accounting."""
import numpy as np
from scipy.optimize import nnls, minimize

def simplex(A,b,weights=None,initial=None,budget=False,pg_index=None,target_mask=None):
    p = A[pg_index]*(np.asarray(target_mask) if target_mask is not None else 1) if budget else None
    if budget:
        assert pg_index is not None and np.min(p)<=b[pg_index]+1e-12
    keep = np.flatnonzero(p<=1e-12) if budget and b[pg_index]<=1e-12 else np.arange(A.shape[1])
    B = A[:,keep]
    weights = np.ones(len(b)) if weights is None else weights
    gram = B.T@(weights[:,None]*B)/len(b)+.001*np.eye(len(keep))
    rhs = B.T@(weights*b)/len(b)
    if initial is None:
        start = nnls(B,b,maxiter=6000)[0]
        start = start/start.sum() if start.sum() else np.ones(len(keep))/len(keep)
    else:
        start = initial[keep].copy()
        start = start/start.sum() if start.sum() else np.ones(len(keep))/len(keep)
    constraints = [{'type':'eq','fun':lambda c:c.sum()-1,'jac':lambda c:np.ones(len(c))}]
    if budget and len(keep)==A.shape[1]:
        p = p[keep]
        limit = b[pg_index]
        if p@start>limit:
            origin = np.zeros(len(keep));origin[np.argmin(p)] = 1
            weight = np.clip((limit-p@origin)/(p@start-p@origin),0,1)
            start = origin+weight*(start-origin)
        constraints.append({'type':'ineq','fun':lambda c:limit-p@c,'jac':lambda c:-p})
    result = minimize(lambda c:.5*c@gram@c-rhs@c,start,jac=lambda c:gram@c-rhs,
        method='SLSQP',bounds=[(0,1)]*len(keep),constraints=constraints,
        options={'ftol':1e-11,'maxiter':700})
    assert result.success,result.message
    c = np.zeros(A.shape[1]);c[keep] = np.maximum(result.x,0)
    assert abs(c.sum()-1)<1e-7
    if budget:
        original_p = A[pg_index]*(np.asarray(target_mask) if target_mask is not None else 1)
        assert original_p@c<=b[pg_index]+1e-7
    return c

def fit(A,b,budget=False,pg_index=None,target_mask=None):
    c = simplex(A,b,budget=budget,pg_index=pg_index,target_mask=target_mask)
    weights = np.ones(len(b))
    for _ in range(3):
        weights = 1/(1+(A@c-b)**2)
        c = simplex(A,b,weights,c,budget,pg_index,target_mask)
    return c,weights,float(np.linalg.norm(A@c-b)/max(np.linalg.norm(b),1e-12))
