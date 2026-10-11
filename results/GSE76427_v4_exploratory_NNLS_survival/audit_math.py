import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

def bh(p):
    p=np.array(p,dtype=float);order=np.argsort(p)
    adjusted=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.empty(len(p));out[order]=np.minimum(1,adjusted)
    return out

def univariate_cox(data):
    time=data.time_months.to_numpy();event=data.event.to_numpy();x=data.high.to_numpy()
    times=np.unique(time[event==1])
    def score_information(beta):
        score,information=0.,0.
        weights=np.exp(beta*x)
        for t in times:
            risk=time>=t;deaths=(time==t)&(event==1);d=deaths.sum()
            r0=weights[risk].sum();r1=(weights[risk]*x[risk]).sum()
            d0=weights[deaths].sum();d1=(weights[deaths]*x[deaths]).sum()
            score+=x[deaths].sum()
            for fraction in np.arange(d)/d:
                denominator=r0-fraction*d0
                first=(r1-fraction*d1)/denominator
                score-=first
                information+=first-first**2  # x is binary, so x**2 == x
        return score,information
    beta=brentq(lambda b:score_information(b)[0],-20,20,xtol=1e-12)
    information=score_information(beta)[1];se=1/np.sqrt(information)
    return np.exp(beta),np.exp(beta-1.959963984540054*se),np.exp(beta+1.959963984540054*se),2*norm.sf(abs(beta/se))
