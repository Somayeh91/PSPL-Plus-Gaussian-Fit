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

def chisq(theta, t, f, f_err):
    t0, u0, tE,f_s = theta
    model = fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    chii = ((f-model)**2*inv_sigma2)
    
    c = (np.cumsum(chii))
    #c = (c)-N
    #cu_chi= np.cumsum(c)
    #cu_chi = cu_chi/max(cu_chi)
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

    
    #more = (t_chi[np.where(t_chi[:,1]>y)[0][0]:10+np.where(t_chi[:,1]>y)[0][0],1]-y).argmin()
    x_m = t_chi[np.where(t_chi[:,1]>y)]
    x2 , y2 =   x_m[0][0], x_m [0][1]
    #t_chi[np.where(t_chi[:,1]>y)[0][0]:10+np.where(t_chi[:,1]>y)[0][0],0][more], t_chi[np.where(t_chi[:,1]>y)[0][0]:10+np.where(t_chi[:,1]>y)[0][0],1][more]
    #less = np.abs(t_chi[np.where(t_chi[:,1]<y)[0][shape(np.where(t_chi[:,1]<y))[1]-1]-10:np.where(t_chi[:,1]<y)[0][shape(np.where(t_chi[:,1]<y))[1]-1],1]-y).argmin()
    x_l = t_chi[np.where(t_chi[:,1]<y)]
    x_l_ = x_l [x_l[:,0]<x2]
    x1 , y1 =   x_l_ [ len(x_l_)-1 ][0], x_l_ [ len(x_l_)-1][1]
    #t_chi[np.where(t_chi[:,1]<y)[0][shape(np.where(t_chi[:,1]<y))[1]-1]-10:np.where(t_chi[:,1]<y)[0][shape(np.where(t_chi[:,1]<y))[1]-1],0][less], t_chi[np.where(t_chi[:,1]<y)[0][shape(np.where(t_chi[:,1]<y))[1]-1]-10:np.where(t_chi[:,1]<y)[0][shape(np.where(t_chi[:,1]<y))[1]-1],1][less]
    return x1,y1,x2,y2
    
    
    
name = 'cassan_32_46_225.det.lc.gz'  
tempdata = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/lcsample/'+str(name)

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
    t,f,f_err,f_true,code = np.loadtxt(str(tempdata),usecols=(0,1,2,3,5),unpack=True)
    t = ma.masked_where(code != 0, t)
    f = ma.masked_where(code != 0  , f)
    f_err = ma.masked_where(code != 0 , f_err)
    u0_true = float(f_s_true)/(max(f)-1+float(f_s_true))
    t0_true = t[f.argmax()]
    ind1, ind2 = fwhm(f,f.argmax(),min(f))
    tE_true = t[ind2]-t[ind1]
#    print tE_true
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [t0_true, u0_true, tE_true,f_s_true], args=(t, f, f_err),method = 'Nelder-Mead',tol = 1e-6)
    t0_ml, u0_ml, tE_ml,f_s_ml = result['x']
    #print result

#print t0_ml, u0_ml, tE_ml,f_s_ml

#plt.close()
#plt.title(str(name))
#plt.plot(t,f,'r*',alpha=0.9)
#plt.plot(t,fun(t0_ml, u0_ml, tE_ml,f_s_ml),'b')
#end = time.time()
#print end-start
#plt.show()

N = np.arange(1,len(t)+1)
chi_2 = chisq([t0_ml, u0_ml, tE_ml,f_s_ml],t,f,f_err)
t_chi = np.column_stack((t,((chi_2)-N)/max((chi_2)-N)))

#sl = [slope(t_chi[i][0],t_chi[i][1],t_chi[i+1][0],t_chi[i+1][1]) for i in range(len(f)-1)]
#sl = np.append(sl, sl[len(sl)-1])
#chi_low , chi_high, duration = [], [], []

#chi_50 = lin_fit(slope(x_y_finder(t_chi,0.5)),0.5)
#for j in arange(10,45,5):
    #chi_low.append(lin_fit(slope(x_y_finder(t_chi,0.5-(j/100))),0.5-(j/100)))
    #chi_high.append(lin_fit(slope(x_y_finder(t_chi,0.5+(j/100))),0.5+(j/100)))
    #duration.append((lin_fit(slope(x_y_finder(t_chi,0.5+(j/100))),0.5+(j/100))-lin_fit(slope(x_y_finder(t_chi,0.5-(j/100))),0.5-(j/100)))/2)
    
jj = np.linspace(0.1,0.49,30,endpoint=True)
duration = [(lin_fit(slope(x_y_finder(t_chi,0.5+j)),0.5+j)-lin_fit(slope(x_y_finder(t_chi,0.5-j)),0.5-j))/2 for j in jj]

#print jj[(np.abs(duration-np.sqrt(float(q))*tE_ml)).argmin()]

#print chi_15,chi_50,chi_85
print 'tE,p exp is:'+ str(np.mean(duration))
print 'tE,p theo is: '+str(np.sqrt(float(q))*tE_ml)


#myRange = [4136,44354]
#filtered = [x for x in t_chi if myRange[0] <= x[1] <= myRange[1]]

#print np.median([i[1] for i in filtered])

plt.close()

#plt.plot(t,chi_2,'r.')
#plt.plot(t,N,'b.')sl[0]
plt.plot(jj, duration, 'b.')
#plt.plot(t,slope(t,chi_2), 'r.', t,slope(t,N), 'b.')
#plt.plot(t,x,'r.')
#plt.plot(t,sl,'b-')
#plt.plot(t,((chi_2)-N)/max((chi_2)-N),'b.')
#plt.axvline(x=initial, color='#1a9641',linewidth=2)
#plt.axvline(x=final, color='#1a9641',linewidth=2)
#plt.plot(t,f_true,'k.')
plt.axhline(y=np.sqrt(float(q))*tE_ml, color='red',linewidth=2)
#plt.axhline(y=0.5, color='red',linewidth=2)
#plt.axhline(y=0.5+0.345, color='red',linewidth=2)
#plt.axhline(y=0.5-0.345, color='red',linewidth=2)
#plt.plot(t,f,'b',t,fun(t0_ml, u0_ml, tE_ml,f_s_ml),'r-')
end = time.time()
print (end-start)
plt.show()
#for i in arange(len(t)):
    #chi_sq.append(chisq([t0_ml, u0_ml, tE_ml,f_s_ml],t_,f_,f_err_))
    
#plt.close()
#plt.plot(t,chi_sq)
    