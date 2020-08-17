import glob,os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
from numpy import *
import re
import scipy.stats as st
from os.path import expanduser
import cmath
import emcee
import scipy.optimize as op
import corner
import time

start = time.time()

home = os.path.expanduser("~")

def lnlike(theta, t, f, f_err):
    t0, u0, tE,f_s = theta
    u = np.sqrt(u0**2+((t-t0)/tE)**2)
    A = ((u**2)+2)/(u*np.sqrt(u**2+4))
    model = (f_s * (A-1)) +1
    inv_sigma2 = 1.0/(f_err**2)
    return -0.5*(np.sum((f-model)**2*inv_sigma2))
    
tempdata = home+'/Library/Mobile Documents/com~apple~CloudDocs/OSU trip/Matt/cassan_32_4_1923.det.lc'

t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)

t = ma.masked_where(code != 0, t)
f = ma.masked_where(code != 0  , f)
f_err = ma.masked_where(code != 0 , f_err)

f_s_true = 0.82659
u0_true = f_s_true/(max(f)-1+f_s_true)
#n = [i for i, elem in enumerate(f) if max(f)==elem]
t0_true = 538 #t[int(n[0])]
tE_true = 6

nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [t0_true, u0_true, tE_true,f_s_true], args=(t, f, f_err),method = 'Nelder-Mead',tol = 1e-6)
t0_ml, u0_ml, tE_ml,f_s_ml = result['x']
print result

def lnprior(theta):
    t0, u0, tE ,f_s = theta
    if 1 < t0 < 2000 and -3 < u0 < 3 and -1 < np.log10(tE) < 3 and 0 < f_s <1:
        return 0.0
    return -np.inf

def lnprob(theta, t, f, f_err):
   lp = lnprior(theta)
   if not np.isfinite(lp):
       return -np.inf
   return lp + lnlike(theta, t, f, f_err)  

#ndim, nwalkers = 4, 2
#pos = [result['x'] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(t, f, f_err))
#sampler.run_mcmc(pos, 500)
#samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
#fig = corner.corner(samples, labels=["$t0$", "$u0$", "$tE$"],
 #                     truths=[t0_true, u0_true, tE_true])
#fig.savefig("triangle.png")

#plt.show()

#samples[:, 2] = np.exp(samples[:, 2])
#t0_mcmc, u0_mcmc, tE_mcmc , f_s_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
 #                            zip(*np.percentile(samples, [16, 50, 84],
  #                                              axis=0)))
                                                
#print t0_mcmc[0], u0_mcmc[0], tE_mcmc[0], f_s_mcmc[0]

u = np.sqrt(u0_ml**2+((t-t0_ml)/tE_ml)**2)
A = ((u**2)+2)/(u*np.sqrt(u**2+4))
f_model = (f_s_ml * (A-1)) +1

plt.close()
plt.plot(t,f-f_model,'r',t,f_true-f_model,'b')
#plt.plot(t,f_model)
end = time.time()
print end-start
plt.show()
