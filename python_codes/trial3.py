import glob,os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
from numpy import *
import re
import scipy.stats as st
from os.path import expanduser
import cmath
import scipy.optimize as op
import time
import gzip

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

def lnlike2(theta, t, f, f_err):
    mean, sigma, amp,t0, u0, tE,f_s = theta
    model = amp*(1/sigma*np.sqrt(2*pi))*np.exp(-((t-mean)**2)/(2*(sigma**2)))+fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    return -0.5*(np.sum((f-model)**2*inv_sigma2))
    
name = 'cassan_32_4_3072.det.lc.gz'  
tempdata = home+'/Library/Mobile Documents/com~apple~CloudDocs/OSU trip/Matt/lcsample/'+str(name)

t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)

t = ma.masked_where(code != 0, t)
f = ma.masked_where(code != 0  , f)
f_err = ma.masked_where(code != 0 , f_err)
file_content = []
if tempdata.endswith('.gz'):
    fname = gzip.open(str(tempdata), 'rb')
    f_s_true = fname.readline().split(' ')[1]
    x = fname.readlines()[4].split(' ')
    q = x[5]
    s = x[6]
    print q,s
    t,f,f_err,f_true,code = np.loadtxt(str(tempdata),usecols=(0,1,2,3,5),unpack=True)
    t = ma.masked_where(code != 0, t)
    f = ma.masked_where(code != 0  , f)
    f_err = ma.masked_where(code != 0 , f_err)
    
    
    u0_true = float(f_s_true)/(max(f)-1+float(f_s_true))
    t0_true = t[f.argmax()]
    ind1, ind2 = fwhm(f,f.argmax(),min(f))
    tE_true = t[ind2]-t[ind1]
    
    
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [t0_true, u0_true, tE_true,f_s_true], args=(t, f, f_err),method = 'Nelder-Mead',tol = 1e-6)
    t0_ml, u0_ml, tE_ml,f_s_ml = result['x']
    #print t0_ml, u0_ml, tE_ml,f_s_ml


f_ris = f-fun(t0_ml, u0_ml, tE_ml,f_s_ml)
f_ris_true = f_true-fun(t0_ml, u0_ml, tE_ml,f_s_ml)


t_mean = t[f_ris.argmax()]
ind1, ind2 = fwhm(f_ris,f_ris.argmax(),min(f_ris))
sigma = t[ind2]-t[ind1]
amp = max(f_ris)-min(f_ris)

print t_mean,sigma,amp


nll = lambda *args: -lnlike2(*args)
result = op.minimize(nll, [t_mean, sigma, amp,t0_ml, u0_ml, tE_ml,f_s_ml], args=(t, f, f_err),method = 'Nelder-Mead',tol = 1e-6)
mean_ml, sigma_ml, amp_ml,t0_ml2, u0_ml2, tE_ml2,f_s_ml2 = result['x']
print mean_ml, sigma_ml, amp_ml,t0_ml2, u0_ml2, tE_ml2,f_s_ml2

print (sigma_ml/tE_ml2)**2, np.abs(mean_ml-t0_ml2)/tE_ml2

plt.close()
plt.figure(1)
plt.plot(t, fun(t0_ml, u0_ml, tE_ml,f_s_ml))
plt.plot(t,f,'r')
plt.figure(2)
plt.plot(t, amp_ml*(1/sigma_ml*np.sqrt(2*pi))*np.exp(-((t-mean_ml)**2)/(2*(sigma_ml**2)))+fun(t0_ml, u0_ml, tE_ml,f_s_ml))
plt.plot(t,f,'r')
plt.figure(3)
plt.plot(t, f-(amp_ml*(1/sigma_ml*np.sqrt(2*pi))*np.exp(-((t-mean_ml)**2)/(2*(sigma_ml**2)))+fun(t0_ml, u0_ml, tE_ml,f_s_ml)))
#plt.plot(t,f,'r')
plt.show()


#(amp_ml*(1/sigma_ml*np.sqrt(2*pi))*np.exp(-((t-mean_ml)**2)/(2*(sigma_ml**2)))+
