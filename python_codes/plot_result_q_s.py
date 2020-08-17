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


def med_med (true,fitted):
    temp = fitted - true
    return (np.median(np.abs(temp-np.median(temp))))

def rms (true,fitted):
    temp = fitted - true
    return np.sqrt((np.sum(temp**2))/len(temp))
    
start = time.time()

home = os.path.expanduser("~")

direc = os.listdir(".")
'''
name = 'alllc_second_run_A_>_1.CSV'  
#Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/result_file/
tempdata = home+'/Desktop/trial_runs/alllce1_A_>_1_fails/'+str(name)
'''

name = 'alllc_result_v2.CSV'  
#Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/result_file/
tempdata = home+'/Desktop/trial_runs/'+str(name)

df = pd.read_csv(tempdata)
df['u0_true'] = np.abs(df['u0_true'])
df['u0_fitted'] = np.abs(df['u0_fitted'])

#(df['chi_2_2']>-22500)& (df['s_fitted']<5) & 
#df = df[((df['f_s_true']*( (2 + df['u0_true']**2) / (df['u0_true']*np.sqrt(4 + df['u0_true']**2)) ) + (1-df['f_s_true']) )>1.5) & (df['u0_true']>0.1) ]

#df = df[ (df['chi_2_2']>-25000) & (df['s_fitted']<5)&((df['f_s_true']*( (2 + df['u0_true']**2) / (df['u0_true']*np.sqrt(4 + df['u0_true']**2)) ) + (1-df['f_s_true']) )>1.1)]

color = ['#e41a1c', '#f781bf', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628','#e7298a', '#e6ab02']
print len(df)

df['q_fitted_log'][df['q_fitted_log']>0] = np.log10( 1/df['q_fitted'][df['q_fitted_log']>0] )
df['q_fitted'][df['q_fitted_log']>0] = ( 1/df['q_fitted'][df['q_fitted_log']>0] )

#df = df[((-2*df['chi_2_1']/41039)<1.1)][ (np.abs((-2*df['chi_2_1']/41039)-(-2*df['chi_2_2']/41039))>0.02)] #was one of earlier constraints that doesn't look good now!
df_bad = df[df['s_fitted']<10][np.abs(df['u0_fitted'])<0.045] # 02/20/2018 looks better than most
df = df[df['s_fitted']<10][np.abs(df['u0_fitted'])>0.045] # 02/20/2018 looks better than most
#df = df[(df['ampl_fitted']<0.001)&(df['ampl_fitted']>-0.001)] # 06/12/2018
print len(df)

df_new_del_t = df[np.abs(df['u0_fitted'])>0.045]
df_new_ampl = df[(np.abs(df['t0_fitted']-df['tp_fitted'])>1)]
df_new_ampl_del_t = df_new_del_t[(df_new_del_t['ampl_fitted']>0.016) | (df_new_del_t['ampl_fitted']<-0.016)]
# cut-offs for s_fitted
low_cut = 0.28
high_cut = 5
err_1 = 0.25
err_2 = 0.1
err_1_q = 1



plt.close('all')

f, axarr = plt.subplots(1, 2)
#f.suptitle('Plots of Fitted physical parameters "q" and "s" Versus True parameters for '+str(len(df))+' targets',size=20)

#axarr[0].set_title('Projected Separation' ,size=26)

axarr[0].plot (df['s_true'],df['s_fitted'],'b.',label='_nolegend_',markersize=8 , alpha = 0.4)
#axarr[0].plot (df['s_true'][(df['ampl_fitted']>0.5)|(df['ampl_fitted']<-0.5)],df['s_fitted'][(df['ampl_fitted']>0.5)|(df['ampl_fitted']<-0.5)],'r.',label='_nolegend_',markersize=8 , alpha = 0.4)

axarr[0].plot ((0,5),(0,5),'g-',label='_nolegend_')
axarr[0].set_xlim((0,5))
axarr[0].set_ylim((0,5))
#plt.legend()
axarr[0].set_xlabel('True Projected Separation ($R_E$)',size=22)
axarr[0].set_ylabel('Fitted Projected Separation ($R_E$)',size=22)
axarr[0].grid()


#axarr[1].set_title('Mass Ratio ',size=26)
axarr[1].loglog (df['q_true'],df['q_fitted'],'b.',markersize=10,label='_nolegend_',alpha=0.4)
#axarr[1].plot (df['q_true_log'][(df['ampl_fitted']>0.5)|(df['ampl_fitted']<-0.5)],df['q_fitted_log'][(df['ampl_fitted']>0.5)|(df['ampl_fitted']<-0.5)],'r.',markersize=10,label='_nolegend_',alpha=0.4)

axarr[1].loglog (df['q_true'],df['q_true'],'g-',label='_nolegend_')
axarr[1].set_xlim((10**-9,1))
axarr[1].set_ylim(( 10**-9,1))

axarr[1].set_xlabel('True Mass Ratio',size=22)
axarr[1].set_ylabel('Fitted Mass Ratio',size=22)
#plt.legend()
axarr[1].grid()





f.set_size_inches(15.0,7)
f.savefig(home+'/Desktop/result_q_s.png')
plt.show()
