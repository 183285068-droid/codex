import numpy as np
from scipy.stats import spearmanr
from analyze import rho_batch
rng=np.random.default_rng(3);x=np.r_[np.zeros(70),np.arange(1,31)/30];xp=np.array([rng.permutation(x) for _ in range(999)])
# A score made only from the shared predictor must yield no excess signal.
null=rho_batch(xp,xp/4);assert np.allclose(null,1)
# Independent scipy calculations must agree with the batched rank calculation.
other=rng.random(len(x));v=rho_batch(xp,(xp+other)/4)
for i in [0,7,101,998]:assert np.isclose(v[i],spearmanr(xp[i],(xp[i]+other)/4).statistic)
# Opposing, varying RNA programs can have a negative excess even when the shared PGAM5 term remains.
y=(x+3*(1-x))/4;observed=spearmanr(x,y).statistic;null=rho_batch(xp,(xp+3*(1-x))/4);assert observed < np.median(null)
print('Self-only null, independent scipy correlation, and opposing-program sanity checks passed.')
