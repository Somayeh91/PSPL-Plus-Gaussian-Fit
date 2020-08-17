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

direc = os.listdir(".")

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
    model = amp*(1/np.sqrt(2*pi*(sigma**2)))*np.exp(-((t-mean)**2)/(2*(sigma**2)))+fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    return -0.5*(np.sum((f-model)**2*inv_sigma2))
    
mean_s = [[],[]]
sigma_s = [[],[]]
c = len(direc)
for i in range(40):

    
    file_content = []
    if direc[i].endswith('.gz'):
        t,f,f_err,f_true,code = np.loadtxt(direc[i],usecols=(0,1,2,3,5),unpack=True)

        t = ma.masked_where(code != 0, t)
        f = ma.masked_where(code != 0  , f)
        f_err = ma.masked_where(code != 0 , f_err)
        fname = gzip.open(str(direc[i]), 'rb')
        f_s_true = fname.readline().split(' ')[1]
        x = fname.readlines()[4].split(' ')
        q = x[5]
        s = x[6]
        t,f,f_err,f_true,code = np.loadtxt(str(direc[i]),usecols=(0,1,2,3,5),unpack=True)
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



        f_ris = f-fun(t0_ml, u0_ml, tE_ml,f_s_ml)
        f_ris_true = f_true-fun(t0_ml, u0_ml, tE_ml,f_s_ml)
        
# making sure we are finding time of maximum absolute deviation:
        if np.abs(max(f_ris))>np.abs(min(f_ris)): 
            t_mean = t[f_ris.argmax()]
        elif np.abs(min(f_ris))>np.abs(max(f_ris)): 
            t_mean = t[f_ris.argmin()]
            
        ind1, ind2 = fwhm(f_ris,f_ris.argmax(),min(f_ris))
        sigma = t[ind2]-t[ind1]
        amp = max(f_ris)
        if (len(direc)-i)% 20 == 0:
            print (len(direc)-i)
    
    
    

        nll = lambda *args: -lnlike2(*args)
        result = op.minimize(nll, [t_mean, sigma, amp,t0_true, u0_true, tE_true,f_s_true], args=(t, f, f_err),method = 'Nelder-Mead',tol = 1e-6)
        mean_ml, sigma_ml, amp_ml,t0_ml, u0_ml, tE_ml,f_s_ml = result['x']
        
        mean_s[0].append(np.abs(mean_ml-t0_ml)/tE_ml)
        mean_s[1].append(s)
        sigma_s[0].append((sigma_ml/tE_ml)**2)
        sigma_s[1].append(q)



plt.close()
plt.figure(1)
plt.title('Separation between the planet and the lens star')
plt.plot(mean_s[1],mean_s[0],'o')
plt.plot(mean_s[1],mean_s[1],'k-')
#plt.axis([0,max(mean_s[1]), 0, max(mean_s[1])])
plt.xlabel('True Separation',size=15)
plt.ylabel('Estimated Separation',size=15)
#fig = plt.gcf()
#fig.set_size_inches(12.0,6.0)
#fig.savefig('./../lc_check/'+str(i.split('.det')[0])+'_1.png')
plt.figure(2)
plt.title('Mass Ratio')
plt.plot(sigma_s[1],sigma_s[0],'o')
plt.plot(sigma_s[1],sigma_s[1],'k-')
#plt.axis([0,max(sigma_s[1]), 0, max(sigma_s[1])])
plt.xlabel('True Mass Ratio',size=15)
plt.ylabel('Estimated Mass Ratio',size=15)
#fig = plt.gcf()
#fig.set_size_inches(12.0,6.0)
#fig.savefig('./../lc_check/'+str(i.split('.det')[0])+'_2.png') 

plt.show()
#(amp_ml*(1/sigma_ml*np.sqrt(2*pi))*np.exp(-((t-mean_ml)**2)/(2*(sigma_ml**2)))+
