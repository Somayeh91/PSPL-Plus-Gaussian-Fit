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
tempdata = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/allc/'+str(name)
dff = pd.read_csv(tempdata)

direc = os.listdir(".")


def fun (t0,u0,tE,f_s):
    u = np.sqrt(u0**2+((df['t']-t0)/tE)**2)
    A = ((u**2)+2)/(u*np.sqrt(u**2+4))
    return (f_s * (A-1)) +1
        
        
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
        tE = x_0[6].split(' ')[4]
        t0 = x_0[6].split(' ')[3]
        u0 = x_0[6].split(' ')[1]
        del_t = (np.sqrt(np.abs((float(s)-(1/float(s)))**2 - (float(u0))**2))) * float(tE)

        F_t = fun(t0,u0,tE,f_s)    
        
        f_ris = df['f']-fun(t0_ml, u0_ml, tE_ml,f_s_ml)
        f_ris_true = df['f_true']-F_t

        #u0_true = float(f_s_true)/(max(df['f']-f_ris)-1+float(f_s_true))
        if np.abs(max(f_ris_true))> np.abs(min(f_ris_true)):
            t_s = df['t'][f_ris_true.argmax()]
        else:
            t_s = df['t'][f_ris_true.argmin()]


       
            
            
result = pd.DataFrame({'name': name_lc ,'q_true_log': q_[0],'q_fitted_log': q_[1],'s_true_log': s_[0],'s_fitted_log': s_[1],
                    'q_true' : 10**np.asarray((q_[0])) , 'q_fitted' : 10**np.asarray((q_[1])), 's_true' : 10**np.asarray(s_[0]), 's_fitted' : 10**np.asarray(s_[1]),
                    'tp_fitted' : tp_ , 'tEp_true' : tEp_[0], 'tEp_fitted' : tEp_[1], 'ampl_fitted' : ampli_, 'del_t' : del_t_, 's_>_1' : s1,
                    'tE_true' : tE_list[0], 'tE_fitted' : tE_list[1], 't0_true' : t0_[0], 't0_fitted' : t0_[1], 'x_c' : x_c_, 's_<_1' : s2,
                    'u0_true' : u0_[0], 'u0_fitted' : u0_[1], 'f_s_true' : f_s_[0], 'f_s_fitted': f_s_[1], 'chi_2_2' : chi_2_2, 'chi_2_1' : chi_2_1})
                    
failure_list = pd.DataFrame({'name' : failure})                    
                    
result.to_csv(home+'/Desktop/trial_runs/alllc1_result.CSV')
failure_list.to_csv(home+'/Desktop/trial_runs/alllc1_fails.CSV')


plt.close()
plt.figure(3)
plt.title('Mass Ratio_Method1')
plt.plot(result['q_true_log'],result['q_true_log'],'g-',result['q_true_log'],result['q_fitted_log'],'b.')
plt.xlabel('True Mass Ratio - log',size=15)
plt.ylabel('Fitted Mass Ratio - log',size=15)
fig = plt.gcf()
fig.set_size_inches(12.0,6.0)
fig.savefig(home+'/Desktop/trial_runs/alllc1/Mass Ratio.png')
plt.close()
plt.figure(4)
plt.title('Projected Separation_Method1')
plt.plot(result['s_true_log'],result['s_true_log'],'g-',result['s_true_log'],result['s_fitted_log'],'b.')
plt.xlabel('True Separation - log',size=15)
plt.ylabel('Fitted Separation - log',size=15)
fig = plt.gcf()
fig.set_size_inches(12.0,6.0)
fig.savefig(home+'/Desktop/trial_runs/alllc1/Projected Separation.png')



end = time.time()
print (end-start)/(60*60)




 