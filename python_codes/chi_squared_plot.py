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
import gzip
from scipy.stats import chi2
from scipy import signal

start = time.time()

home = os.path.expanduser("~")

def fwhm(valuelist, peakpos,base):
    peakvalue = valuelist[peakpos]-base
    phalf = (peakvalue / 2.0)+base

    # go left and right, starting from peakpos
    ind1 = peakpos
    ind2 = peakpos   

    while ind1>2 and valuelist[ind1]>phalf:
        ind1=ind1-1
    while ind2<len(valuelist)-1 and valuelist[ind2]>phalf:
        ind2=ind2+1  
    return ind1,ind2

def fun (t0,u0,tE,f_s):
    u = np.sqrt(u0**2+((t-t0)/tE)**2)
    A = ((u**2)+2)/(u*np.sqrt(u**2+4))
    return (f_s * (A-1)) +1
        
def lnlike(theta, t, f, f_err):
    t0, u0, tE,f_s = theta
    model = fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    return -0.5*(np.sum((f-model)**2*inv_sigma2))

def chisq_cum(theta, t, f, f_err):
    t0, u0, tE,f_s = theta
    model = fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    chii = (f-model)**2*inv_sigma2
    
    c = (np.cumsum(chii))
    #c = (c)-N
    #cu_chi= np.cumsum(c)
    #cu_chi = cu_chi/max(cu_chi)
    return c
    
def chisq(theta, t, f, f_err):
    t0, u0, tE,f_s = theta
    model = fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    chii = (f-model)**2*inv_sigma2
        #c = (c)-N
    #cu_chi= np.cumsum(c)
    #cu_chi = cu_chi/max(cu_chi)
    return chii

name = 'cassan_32_47_2531.det.lc.gz'  
tempdata = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/lcsample/'+str(name)
#'/Desktop/alllc1/'+str(name) #

t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)

t = ma.masked_where(code != 0, t)
f = ma.masked_where(code != 0  , f)
f_err = ma.masked_where(code != 0 , f_err)
file_content = []
if tempdata.endswith('.gz'):
    fname = gzip.open(str(tempdata), 'rb')
    f_s_true = fname.readline().split(' ')[1]
    t,f,f_err,f_true,code = np.loadtxt(str(tempdata),usecols=(0,1,2,3,5),unpack=True)
    t = ma.masked_where(code != 0, t)
    f = ma.masked_where(code != 0  , f)
    f_err = ma.masked_where(code != 0 , f_err)
    u0_true = float(f_s_true)/(max(f)-1+float(f_s_true))
    t0_true = t[f.argmax()]
    ind1, ind2 = fwhm(f,f.argmax(),min(f))
    tE_true = t[ind2]-t[ind1]
    print tE_true
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [t0_true, u0_true, tE_true,f_s_true], args=(t, f, f_err),method = 'Nelder-Mead',tol = 1e-6)
    t0_ml, u0_ml, tE_ml,f_s_ml = result['x']
    #print result
print t0_ml, u0_ml, tE_ml,f_s_ml


f_ris = f-fun(t0_ml, u0_ml, tE_ml,f_s_ml)
t0_true = t[(f-f_ris).argmax()]

#plt.close()
#plt.title(str(name))
#plt.plot(t,f,'r*',alpha=0.9)
#plt.plot(t,fun(t0_ml, u0_ml, tE_ml,f_s_ml),'b')
#end = time.time()
#print end-start
#plt.show()

N = np.arange(1,len(t)+1)

    
plt.close()
chi_2_cum = chisq_cum([t0_ml, u0_ml, tE_ml,f_s_ml],t,f,f_err)
chi_2 = chisq([t0_ml, u0_ml, tE_ml,f_s_ml],t,f,f_err)

t_chi = np.column_stack((t,(np.abs((chi_2)-N))/max(np.abs((chi_2)-N))))


half_chi_ = [t_chi[(f_ris).argmax(),1],t_chi[(f_ris).argmin(),1]]
half_chi = half_chi_[np.asarray([np.abs(l - 0.5) for l in half_chi_]).argmin()]


t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)

t = ma.masked_where(code != 0, t)
f = ma.masked_where(code != 0  , f)
f_err = ma.masked_where(code != 0 , f_err)


plt.title('$\chi^2$ - N')
plt.plot(t,(chi_2_cum-N)/max(chi_2_cum-N),'k.')#,t,chi_2,'b.')
#np.abs(chi_2_cum-N)/max(np.abs(chi_2_cum-N))
#plt.plot(t,N,'b.')
plt.xlabel('time')
plt.ylabel('$\chi^2$ - N')
#plt.plot(t,chi_2-N,'r.')
#plt.plot(t,chi_N,'k-')
#plt.axvline(x=t[chi_2.argmax()], color='r',linewidth=2)
plt.axhline(y=half_chi, color='r',linewidth=2)
#plt.plot(t,f,'b.',t,fun(t0_ml, u0_ml, tE_ml,f_s_ml),'r')
end = time.time()
print end-start
plt.show()
#for i in arange(len(t)):
    #chi_sq.append(chisq([t0_ml, u0_ml, tE_ml,f_s_ml],t_,f_,f_err_))
    
#plt.close()
#plt.plot(t,chi_sq)
    