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

def tE_finder (t,f,f_s):
    
    A_base = (float(f_s)*0.34)+1
    
    t_max = t[f.argmax()]
    if t[t>t_max] != []:
        tE_right = np.abs(t[t>t_max][(np.abs(f[t>t_max]-(A_base))).argmin()]-t_max)
    else:
        tE_right = 1
    if t[t<t_max] != []:
        tE_left = np.abs(t[t<t_max][(np.abs(f[t<t_max]-(A_base))).argmin()]-t_max)
    else:
        tE_left = 1
    
    return min([tE_right,tE_left])

def fun (t0,u0,tE,f_s):
    u = np.sqrt(u0**2+((t-t0)/tE)**2)
    A = ((u**2)+2)/(u*np.sqrt(u**2+4))
    return (f_s * (A-1)) +1
        
def fun2 (mean, sigma,amp, t0,u0,tE,f_s):
    u = np.sqrt(u0**2+((t-t0)/tE)**2)
    A = (((amp/np.sqrt(2*pi*(sigma**2)))*np.exp(-((t-mean)**2)/(2*(sigma**2)))))+((u**2)+2)/(u*np.sqrt((u**2)+4))
    return (f_s * (A-1)) +1

def lnlike(theta, t, f, f_err):
    t0, u0, tE,f_s = theta
    model = fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    return -0.5*(np.sum((f-model)**2*inv_sigma2))

def lnlike2(theta, t, f, f_err):
    mean, sigma,amp, t0,u0,tE,f_s = theta
    model = fun2(mean, sigma,amp, t0,u0,tE,f_s)
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
        
ot = open(home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/result_file/'+str(len(direc))+'.dat', 'w')
        
s_ = [[],[]]
q_ = [[],[]]

for i in range(len(direc)):
    
    tempdata = direc[i]
    print tempdata

    if tempdata.endswith('.gz'):
        t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)
        t = ma.masked_where(code != 0, t)
        f = ma.masked_where(code != 0  , f)
        f_err = ma.masked_where(code != 0 , f_err)
        fname = gzip.open(tempdata, 'rb')
        x_0 = fname.readlines()[0:7]
        f_s_true = x_0[0].split(' ')[1]
        q = x_0[5].split(' ')[5]
        s = x_0[5].split(' ')[6]
        tE_theo = x_0[6].split(' ')[4]

        t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)
    
        t = ma.masked_where(code != 0, t)
        f = ma.masked_where(code != 0  , f)
        f_err = ma.masked_where(code != 0 , f_err)
        
    
        u0_true = float(f_s_true)/(max(f)-1+float(f_s_true))
        t0_true = t[f.argmax()]
        ind1, ind2 = fwhm(f,f.argmax(),min(f))
        #tE_true = t[ind2]-t[ind1]
        tE_true = [tE_finder (t,f, f_s_true),t[ind2]-t[ind1]]
        #print 'tE_true = '+ str( tE_true)

        tE_ = [[],[]]    
        for i in tE_true:
               
            nll = lambda *args: -lnlike(*args)
            result = op.minimize(nll, [t0_true, u0_true, i,f_s_true], args=(t, f, f_err),method = 'Nelder-Mead')
            t0_ml, u0_ml, tE_ml,f_s_ml = result['x']
            tE_[0].append(lnlike([t0_ml, u0_ml, tE_ml,f_s_ml],t,f,f_err))
            tE_[1].append([t0_ml, u0_ml, tE_ml,f_s_ml])
    
        #print tE_[1]
        mm = np.asarray( tE_[0])
        tE__ = tE_[1][mm.argmax()]

        t0_ml, u0_ml, tE_ml,f_s_ml = tE__[0],tE__[1],tE__[2],tE__[3]

        f_ris = f-fun(t0_ml, u0_ml, tE_ml,f_s_ml)
        f_ris_true = f_true-fun(t0_ml, u0_ml, tE_ml,f_s_ml)

        u0_true = float(f_s_true)/(max(f-f_ris)-1+float(f_s_true))
        t0_true = t[(f-f_ris).argmax()]

        #print 't0_true = ' + str(t0_true)
        

        
        N = np.arange(1,len(t)+1)
        chi_2 = chisq([t0_ml, u0_ml, tE_ml,f_s_ml],t,f,f_err)
        t_chi = np.column_stack((t,(np.abs((chi_2)-N))/max(np.abs((chi_2)-N))))
        
        half_chi_ = [t_chi[(f_ris).argmax(),1],t_chi[(f_ris).argmin(),1]]
        t_chi_ =  [t_chi[(f_ris).argmax(),0],t_chi[(f_ris).argmin(),0]] 
        f_ris__ = [(f_ris).max(),(f_ris).min()]
        half_chi = half_chi_[np.asarray([np.abs(l - 0.5) for l in half_chi_]).argmin()]
        t_half_chi = t_chi_[np.asarray([np.abs(l - 0.5) for l in half_chi_]).argmin()]
        f_ris_max = f_ris__[np.asarray([np.abs(l - 0.5) for l in half_chi_]).argmin()] 

        jj = [0.25,0.35] #np.linspace(0.15,1-half_chi,4)
        #print jj

        ttt =t_half_chi    
        duration = [min([np.abs((lin_fit(slope(x_y_finder(t_chi,half_chi+j)),half_chi+j))-ttt ), np.abs(ttt-(lin_fit(slope(x_y_finder(t_chi,half_chi-j)),half_chi-j)))]) for j in jj]
        
        min_model__ = [[],[]]
        min_model_ = [[],[]]
                
        t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)

        t = ma.masked_where(code != 0, t)
        f = ma.masked_where(code != 0  , f)
        f_err = ma.masked_where(code != 0 , f_err)

        

        # making sure we are finding time of maximum absolute deviation:
                

                        

        if np.asarray([np.abs(l - 0.5) for l in half_chi_]).argmin() == 1:
            amp = 'neg'
        else:
            amp = 'pos'
                        

        for sigma in duration:
            
            amp_ = f_ris_max
            t_mean_ = t_half_chi
            amp_ = amp_ * sigma * np.sqrt(2*pi)
            #print 'sigma = '+str(sigma) + ', amp = '+ str(amp_)
            
            nll = lambda *args: -lnlike2(*args)
            result = op.minimize(nll, [t_mean_,sigma,amp_,t0_ml, u0_ml, tE_ml,f_s_ml], args=(t,f, f_err),method = 'Nelder-Mead')
            mean_mll, sigma_mll,amp_mll,t0_mll, u0_mll, tE_mll,f_s_mll = result['x']
            #print result['x']
            min_model_[0].append(lnlike2([mean_mll, sigma_mll,amp_mll,t0_mll, u0_mll, tE_mll,f_s_mll],t,f,f_err))
            min_model_[1].append([mean_mll, sigma_mll,amp_mll,t0_mll, u0_mll, tE_mll,f_s_mll])

        

        mmm_ = np.asarray( min_model_[0])



        final_param = min_model_[1][mmm_.argmax()]
    
        #print 'second fit =' + str(final_param)
    

    

        x_c = np.sqrt((final_param[4]**2)+((final_param[3]-final_param[0])/final_param[5])**2 )

        if amp == 'pos':
            s_exp = (x_c+np.sqrt((x_c**2)+4))/2   
        elif amp == 'neg':
            s_exp = np.abs((x_c-np.sqrt((x_c**2)+4))/2)
        #print 'theo q = ' +q +',theo s = ' +s
        #print 'exp q = '+str((final_param[1]/final_param[5])**2)+ ', exp s = ' +str(s_exp)
        #print '\n'
        
        if np.log10((final_param[1]/final_param[5])**2) != 'nan' and np.log10(s_exp) != 'nan':
            # theoretical values of mass ratio
            q_[1].append(np.log10(float(q)))
            ot.write(str(np.log10(float(q))))
            ot.write(' ')
            # experimental values of mass ratio
            q_[0].append(np.log10((final_param[1]/final_param[5])**2))
            ot.write(str(np.log10((final_param[1]/final_param[5])**2)))
            ot.write(' ')
            # theoretical values of separation
            s_[1].append(np.log10(float(s)))
            ot.write(str(np.log10(float(s))))
            ot.write(' ')
            # experimental values of separation
            s_[0].append(np.log10(s_exp))
            ot.write(str(np.log10(s_exp)))
            ot.write(' ')
            ot.write('\n')

print q_[0] , q_[1]
print s_[0] , s_[1]

q_theo = np.asarray(q_[1])
s_theo = np.asarray(s_[1])
q_exp = np.asarray(q_[0])
s_exp = np.asarray(s_[0])
    
color = [ '#e41a1c', '#f781bf', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628','#e7298a', '#e6ab02']                                                                                                                                                                                        

color_code = np.linspace(min(q_[1]),max(q_[1]),10)
cdict = {
  'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
  'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
  'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
}

cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
        
plt.close()
plt.figure(1)
plt.title('Mass Ratio')
plt.plot(q_[1],q_[1],'g-',q_[1],q_[0],'b.')
plt.xlabel('Theoretical Mass Ratio - log',size=15)
plt.ylabel('Experimental Mass Ratio - log',size=15)
plt.figure(2)
plt.title('Projected Separation')
plt.xlabel('Theoretical Separation - log',size=15)
plt.ylabel('Experimental Separation - log',size=15)
plt.plot(s_[1],s_[1],'g-')
plt.pcolor(s_[1],s_[0],q_[1],cmap=cm)
plt.colorbar()
end = time.time()
print (end-start)/60
plt.show()







 