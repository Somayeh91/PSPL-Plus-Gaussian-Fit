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
from scipy.stats import chi2
from scipy import signal

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

def chisq(theta, t, f, f_err):
    t0, u0, tE,f_s = theta
    model = fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    chii = ((f-model)**2*inv_sigma2)
    
    c = (np.cumsum(chii))
    return c
def subtract(chi):
    c=0
    r = 0
    for i in range (len(chi)-1):
        if ((chi[i+1]-chi[i])) > c:
            c = chi[i+1]-chi[i]
            r = i
        
    return chi[r+1]
    
def slope(n):
    return float(n[3]-n[1])/float(n[2]-n[0]), n[1]-(float(n[3]-n[1])/float(n[2]-n[0]))*n[0]
    
def lin_fit (line_,y_):
    return (y_ - line_[1])/line_[0]
    
def x_y_finder (t_chi,y):

    
    x_m = t_chi[np.where(t_chi[:,1]>y)]
    if shape(x_m)==(0,2):
        x2,y2 = nan,nan
    else:
        x2 , y2 =   x_m[0][0], x_m [0][1]
    x_l = t_chi[np.where(t_chi[:,1]<y)]
    x_l_ = x_l [x_l[:,0]<x2]
    if shape(x_l_)==(0,2):
        x1,y1 = nan,nan
    else:
        x1 , y1 =   x_l_ [ len(x_l_)-1 ][0], x_l_ [ len(x_l_)-1][1]
    if (x1 == 'nan') or x2 == 'nan':
        return 0,0,0,1
    else:
        return x1,y1,x2,y2

percentile_, mean_tEp, perc = [], [[],[]], []
jj = [0.15,0.25,0.35,0.45]
    
for i in range(40):
    
    tempdata = direc[i]
    print tempdata

    if tempdata.endswith('.gz'):
        t,f,f_err,f_true,code = np.loadtxt(direc[i],usecols=(0,1,2,3,5),unpack=True)

        t = ma.masked_where(code != 0, t)
        f = ma.masked_where(code != 0  , f)
        f_err = ma.masked_where(code != 0 , f_err)
        fname = gzip.open(str(direc[i]), 'rb')
        x_0 = fname.readlines()[0:7]
        f_s_true = x_0[0].split(' ')[1]
        q = x_0[5].split(' ')[5]
        s = x_0[5].split(' ')[6]
        tE_theo = x_0[6].split(' ')[4]
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

    
        N = np.arange(1,len(t)+1)
        chi_2 = chisq([t0_ml, u0_ml, tE_ml,f_s_ml],t,f,f_err)
        t_chi = np.column_stack((t,(np.abs((chi_2)-N))/max(np.abs((chi_2)-N))))
    
    
         
        duration = [(lin_fit(slope(x_y_finder(t_chi,0.5+j)),0.5+j)-lin_fit(slope(x_y_finder(t_chi,0.5-j)),0.5-j))/2 for j in jj]
    
        percentile_.append( jj[(np.abs(duration-(np.sqrt(float(q))*float(tE_theo)))).argmin()])
        
        #perc.append((np.abs(duration[(np.abs(duration-np.sqrt(float(q))*float(tE_theo))).argmin()]-(np.sqrt(float(q))*float(tE_theo)))/(np.sqrt(float(q))*float(tE_theo)))*100)
        mean_tEp[1].append(min(np.median(duration),np.mean(duration)))
        mean_tEp[0].append(np.sqrt(float(q))*float(tE_theo))
        #print chi_15,chi_50,chi_85
        #print 'tE,p exp is:'+ str(np.mean(duration)), 'tE,p theo is: '+str(np.sqrt(float(q))*float(tE_theo))

print (percentile_)    
plt.close()
plt.plot(mean_tEp[0], mean_tEp[1], 'r.')
plt.plot(mean_tEp[0],mean_tEp[0], 'b-')
#plt.axhline(y=np.sqrt(float(q))*tE_ml, color='red',linewidth=2)
#plt.plot(perc,'b.')
end = time.time()
print (end-start)
plt.show()
           