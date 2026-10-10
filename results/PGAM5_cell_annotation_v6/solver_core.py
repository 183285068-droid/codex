"""Fixed v5 nonnegative simplex solver; never selected from v6 outcomes."""
import numpy as np
from scipy.optimize import nnls, minimize


def simplex(A,b,weights=None,initial=None):
    weights = np.ones(len(b)) if weights is None else weights
    gram = A.T@(weights[:,None]*A)/len(b)+.001*np.eye(A.shape[1])
    rhs = A.T@(weights*b)/len(b)
    if initial is None:
        initial = nnls(A,b,maxiter=3000)[0]
        initial = initial/initial.sum() if initial.sum() else np.ones(A.shape[1])/A.shape[1]
    result = minimize(lambda c:.5*c@gram@c-rhs@c,initial,jac=lambda c:gram@c-rhs,
        method='SLSQP',bounds=[(0,1)]*A.shape[1],
        constraints={'type':'eq','fun':lambda c:c.sum()-1,'jac':lambda c:np.ones(len(c))},
        options={'ftol':1e-11,'maxiter':500})
    assert result.success,result.message
    return np.maximum(result.x,0)


def fit(A,b):
    c = simplex(A,b)
    weights = np.ones(len(b))
    for _ in range(3):
        weights = 1/(1+(A@c-b)**2)
        c = simplex(A,b,weights,c)
    return c,weights,float(np.linalg.norm(A@c-b)/max(np.linalg.norm(b),1e-12))
