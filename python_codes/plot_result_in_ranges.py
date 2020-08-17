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
import matplotlib.lines as mlines



def med_med (true,fitted):
    temp = fitted - true
    return (np.median(np.abs(temp-np.median(temp))))


start = time.time()

home = os.path.expanduser("~")

direc = os.listdir(".")

name = 'alllc_result_v2.CSV'  
tempdata = home+'/Desktop/trial_runs/'+str(name)
df = pd.read_csv(tempdata)
df['u0_true'] = np.abs(df['u0_true'])
df['u0_fitted'] = np.abs(df['u0_fitted'])

print len(df)
#df = df[ (df['chi_2_2']>-25000) & (df['s_fitted']<5) & ((df['f_s_true']*( (2 + df['u0_true']**2) / (df['u0_true']*np.sqrt(4 + df['u0_true']**2)) ) + (1-df['f_s_true']) )>1.1)]
#df = df[(np.abs(df['t0_fitted']-df['tp_fitted'])>1) & ((-2*df['chi_2_2']/41039)>1.003) ]

df['q_fitted_log'][df['q_fitted_log']>0] = np.log10( 1/df['q_fitted'][df['q_fitted_log']>0] )
df['q_fitted'][df['q_fitted_log']>0] = ( 1/df['q_fitted'][df['q_fitted_log']>0] )

#df = df[df['s_fitted']<5][(np.abs(df['t0_fitted']-df['tp_fitted'])>1) & ((-2*df['chi_2_2']/41039)>1.003) ]
df = df[df['s_fitted']<10][np.abs(df['u0_fitted'])>0.045] #02/20/2018

df_new_del_t = df

print(len(df))

color = ['#e41a1c', '#f781bf', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628','#e7298a', '#e6ab02']
temp_time = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/allc/time_stamps.CSV'
time_ = pd.read_csv(temp_time)

q_range =  np.array([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,1])
s_range =  np.array([0,0.4,0.8, 0.9 , 1, 1.1, 1.2 , 1.5, 2.5, 3.5, 5])
p_range =  [0.5,5,20,30,40,50,500]


plt.close('all')

f, axarr = plt.subplots(2, 2)
#f.suptitle('Plots of Fitted physical parameters "q" and "s" Versus True parameters',size=25)
i = 7
j= 1
axarr[0, 0].set_title('Mass Ratio',size=26)
#axarr[0, 0].set_xlabel('Log (True Mass Ratio)',size=15)
axarr[0, 0].set_ylabel('Log (Fitted Mass Ratio)',size=22)
axarr[0, 0].plot((-8,0),(-8,0),'g-',label='_nolegend_')
axarr[0, 0].plot(df['q_true_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],df['q_fitted_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],'o',color='g',alpha=0.6,label= '('+str((s_range[i]))+')-('+str((s_range[i+1]))+')')
#plt.axis([min(df['q_fitted_log']),max(df['q_fitted_log']),min(df['q_fitted_log']),max(df['q_fitted_log'])])
axarr[0, 0].set_xlim((-8,0))
axarr[0, 0].tick_params(axis='x',which='both', bottom='off', top='off',labelbottom='off')
axarr[0, 0].set_ylim((-8,0))
axarr[0, 0].text(-7.7, -0.75, 'Separation of $1.5$ to $2.5$',size=15,
        bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
#plt.text(-7.5,-0.5, 'Median Absolute Deviation = '+str(round(med_med(df['q_true_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],df['q_fitted_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])]),6)))
#axarr[0, 0].legend(title='Separation', loc=2,fontsize='x-small')   
axarr[0, 0].grid()
#fig = plt.gcf()
#fig.savefig(home+'/Desktop/'+'('+str((s_range[i]))+')-('+str((s_range[i+1]))+').png')

axarr[0, 1].set_title('Projected Separation',size=26)
#axarr[0, 1].set_xlabel('True Separation ',size=15)
axarr[0, 1].set_ylabel('Fitted Separation ',size=22) 
axarr[0, 1].plot((0,5),(0,5),'g-',label='_nolegend_')
axarr[0, 1].plot(df['s_true'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])],df['s_fitted'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])],'o',color='g',alpha=0.6,label= '('+str((q_range[j]))+')-('+str((q_range[j+1]))+')')
#plt.axis([min(df['s_true']),max(df['s_true']),min(df['s_true']),max(df['s_true'])])
axarr[0, 1].set_xlim((0,5))
axarr[0, 1].set_ylim((0,5))
axarr[0, 1].tick_params(axis='x',which='both', bottom='off', top='off',labelbottom='off')
axarr[0, 1].text(0.2,4.5, 'Mass Ratio of $10^{-6}$ to $10^{-5}$',size=15,
        bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
#plt.text(4,0.5, 'Median Absolute Deviation = '+str(round(med_med(df['s_true'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])],df['s_fitted'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])]),3)))
#axarr[0, 1].legend(title='Mass Ratio', loc=2,fontsize='x-small')   
axarr[0, 1].grid()   

i = 4
j= 4
#axarr[1, 0].set_title('Mass Ratio',size=15)
axarr[1, 0].set_xlabel('Log (True Mass Ratio)',size=22)
axarr[1, 0].set_ylabel('Log (Fitted Mass Ratio)',size=22)
axarr[1, 0].plot((-8,0),(-8,0),'g-',label='_nolegend_')
axarr[1, 0].plot(df['q_true_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],df['q_fitted_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],'o',color='r',alpha=0.6,label= '('+str((s_range[i]))+')-('+str((s_range[i+1]))+')')
plt.axis([min(df['q_fitted_log']),max(df['q_fitted_log']),min(df['q_fitted_log']),max(df['q_fitted_log'])])
axarr[1, 0].set_xlim((-8,0))
axarr[1, 0].set_ylim((-8,0))

#plt.text(-7.5,-0.5, 'Median Absolute Deviation = '+str(round(med_med(df['q_true_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],df['q_fitted_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])]),6)))
axarr[1, 0].text(-7.7, -0.75, 'Separation of $1.0$ to $1.1$',size=15,
        bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
#axarr[1, 0].legend(title='Separation', loc=2,fontsize='x-small') 
axarr[1, 0].grid()
#fig = plt.gcf()
#fig.savefig(home+'/Desktop/'+'('+str((s_range[i]))+')-('+str((s_range[i+1]))+').png')

#axarr[1, 1].set_title('Projected Separation',size=15)
axarr[1, 1].set_xlabel('True Separation ',size=22)
axarr[1, 1].set_ylabel('Fitted Separation ',size=22) 
axarr[1, 1].plot((0,5),(0,5),'g-',label='_nolegend_')
axarr[1, 1].plot(df['s_true'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])],df['s_fitted'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])],'o',color='r',alpha=0.6,label= '('+str((q_range[j]))+')-('+str((q_range[j+1]))+')')
#plt.axis([min(df['s_true']),max(df['s_true']),min(df['s_true']),max(df['s_true'])])
axarr[1, 1].set_xlim((0,5))
axarr[1, 1].set_ylim((0,5))
#plt.text(4,0.5, 'Median Absolute Deviation = '+str(round(med_med(df['s_true'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])],df['s_fitted'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])]),3)))
axarr[1, 1].text(0.2, 4.55, 'Mass Ratio of $0.001$ to $0.01$',size=15,
        bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
#axarr[1, 1].legend(title='Mass Ratio', loc=2,fontsize='x-small')   
axarr[1, 1].grid()  

plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=None, hspace=0.04)  

f.set_size_inches(15.0,10.0)
f.savefig(home+'/Desktop/result_range.png')
end = time.time()
print end-start
plt.show()