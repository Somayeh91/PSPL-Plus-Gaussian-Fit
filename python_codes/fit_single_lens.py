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
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import (mark_inset,inset_axes,InsetPosition)
from tqdm import tqdm

start = time.time()

home = os.path.expanduser("~")

name = 'allc1_info.CSV'  
tempdata = home+'/Desktop/alllc1/'+str(name)
dff = pd.read_csv(tempdata)

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
    df = pd.DataFrame({'t' : t, 'f' : f})
    
    
    A_base = (float(f_s)*0.34)+1
    
    t_max = df['t'][df['f'].argmax()]
    t_right = df['t'][df['t']>t_max]
    t_left = df['t'][df['t']<t_max]
    
    if len(t_right) != 0:
        tE_right = np.abs(t_right[(np.abs(df['f'][df['t']>t_max]-(A_base))).argmin()]-t_max)
    else:
        tE_right = 1
    if len(t_left) != 0:
        tE_left = np.abs(t_left[(np.abs(df['f'][df['t']<t_max]-(A_base))).argmin()]-t_max)
    else:
        tE_left = 1
    
    return min([tE_right,tE_left])

def fun (t0,u0,tE,f_s):
    u = np.sqrt(u0**2+((df['t']-t0)/tE)**2)
    A = ((u**2)+2)/(u*np.sqrt(u**2+4))
    return (f_s * (A-1)) +1
        

def lnlike(theta, t, f, f_err):
    t0, u0, tE,f_s = theta
    model = fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    return -0.5*(np.sum((f-model)**2*inv_sigma2))

        

name_lc = []
chi_2_1 = []


for i in tqdm( range(len(dff['name']))):
    
    tempdata1 = home+'/Desktop/alllc1/'+dff['name'][i]
    #print dff['name'][i]

    if dff['name'][i].endswith('.gz'):
        t,f,f_err,f_true,code = np.loadtxt(tempdata1,usecols=(0,1,2,3,5),unpack=True)
        df = pd.DataFrame({'t':t , 'f':f , 'f_err' : f_err , 'f_true': f_true, 'code':code})
        df = df[df['code']==0]
        
        #df['f'] = df['f_true']
        
        fname = gzip.open(tempdata1, 'rb')
        x_0 = fname.readlines()[0:7]
        f_s_true = x_0[0].split(' ')[1]
        q = x_0[5].split(' ')[5]
        s = x_0[5].split(' ')[6]
        tE_theo = x_0[6].split(' ')[4]
        t0_theo = x_0[6].split(' ')[3]
        u0_theo = x_0[6].split(' ')[1]
        
    
        A_max = 1.0/(float(f_s_true)/(max(df['f'])-1+float(f_s_true)))
        u0_true = np.sqrt( ( (1+np.sqrt(1+16*(A_max**2)))/(2* A_max) ) - 2 )
        t0_true = df['t'][df['f'].argmax()]
        ind1, ind2 = fwhm(df['f'],df['f'].argmax(),min(df['f']))
        #tE_true = t[ind2]-t[ind1]
        tE_true = [tE_finder (df['t'],df[ 'f'], f_s_true),t[ind2]-t[ind1]]
        #print 'tE_true = '+ str( tE_true)

        tE_ = [[],[]]    
        for j in tE_true:
               
            nll = lambda *args: -lnlike(*args)
            result = op.minimize(nll, [t0_true, u0_true, j,f_s_true], args=(df['t'],df[ 'f'], df['f_err']),method = 'Nelder-Mead')
            t0_ml, u0_ml, tE_ml,f_s_ml = result['x']
            tE_[0].append(lnlike([t0_ml, u0_ml, tE_ml,f_s_ml],df['t'],df[ 'f'], df['f_err']))
            tE_[1].append([t0_ml, u0_ml, tE_ml,f_s_ml])
    
        #print tE_[1]
        mm = np.asarray( tE_[0])
        tE__ = tE_[1][mm.argmax()]

        print dff['name'][i], (mm.max()*-2)/len(df['t'])
        
        # name
        name_lc.append(dff['name'][i])
        # minimum chi squared 1
        chi_2_1.append(mm.max())

result = pd.DataFrame({'name': name_lc , 'chi_2_1' : chi_2_1})
                    
#failure_list = pd.DataFrame({'name' : failure})                    
                    
result.to_csv(home+'/Desktop/trial_runs/alllc1_result_method3.CSV')
#failure_list.to_csv(home+'/Desktop/trial_runs/alllc1_fails.CSV')


end = time.time()
print (end-start)/(60*60)




 