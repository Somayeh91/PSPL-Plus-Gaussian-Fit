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
from scipy.interpolate import UnivariateSpline


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
    

for i in direc:
    file_content = []
    if i.endswith('.gz'):
        fname = gzip.open(str(i), 'rb')
        f_s_true = fname.readline().split(' ')[1]
        t,f,f_err,f_true,code = np.loadtxt(str(i),usecols=(0,1,2,3,5),unpack=True)
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
        
        plt.close()
        plt.figure(1)
        #plt.plot(t,f-fun(t0_ml, u0_ml, tE_ml,f_s_ml),'r',alpha=0.9)
        plt.plot(t,f_true-fun(t0_ml, u0_ml, tE_ml,f_s_ml),'b')
        plt.title(str(i.split('.det')[0]))
        plt.ylabel('Residuals')
        plt.axis([t0_ml-tE_ml*1.1,t0_ml+tE_ml*1.1, min(f_true-fun(t0_ml, u0_ml, tE_ml,f_s_ml))*1.1, max(f_true-fun(t0_ml, u0_ml, tE_ml,f_s_ml))*1.1])
        fig = plt.gcf()
        fig.set_size_inches(12.0,6.0)
        fig.savefig('./../residuals/'+str(i.split('.det')[0])+'.png')
        #plt.plot(t,f_model)
        
        
        #plt.figure(2)
        #plt.plot(t,f,'k',t,fun(t0_ml, u0_ml, tE_ml,f_s_ml),'r')
        #plt.axis([t0_ml-tE_ml*2,t0_ml+tE_ml*2, min(f_true), max(f_true)])
        #fig = plt.gcf()
        #fig.set_size_inches(12.0,6.0)
        #fig.savefig('./../residuals/'+str(i.split('.det')[0])+'1'+'.png')
    

