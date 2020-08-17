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

u = df['u0_true']
df['max_A'] = (2+u**2)/(u*np.sqrt(4+u**2))

df = df[df['s_fitted']<10][np.abs(df['u0_fitted'])>0.045] #02/20/2018

#df_new_del_t = df[(np.abs(df['t0_fitted']-df['tp_fitted'])>1) | ((-2*df['chi_2_1']/41039)>1.04) | (np.abs((-2*df['chi_2_1']/41039)-(-2*df['chi_2_2']/41039))>0.02)]
##df_new_del_t = df[(np.abs(df['t0_fitted']-df['tp_fitted'])>1) & ((-2*df['chi_2_2']/41039)>1.003) ]
df_new_del_t = df[np.abs(df['u0_fitted'])>0.045]
df_new_ampl = df[(np.abs(df['t0_fitted']-df['tp_fitted'])>1)]
df_new_ampl_del_t = df_new_del_t[(df_new_del_t['ampl_fitted']>0.016) | (df_new_del_t['ampl_fitted']<-0.016)]
# cut-offs for s_fitted
low_cut = 0.28
high_cut = 5
err_1 = 0.25
err_2 = 0.1
err_1_q = 1

x1 = df['s_fitted']/df['s_true']
y1 = df['ampl_fitted']
x1_name ='Ratio of fitted separation over true separation' 
y1_name = 'Amplitude of the planetary perturbation'

# Regionalizing the plot of s_fitted vs s_true

s_deg = (np.abs(df['s_true']-np.sqrt((df['s_true']-(1/df['s_true']))**2 + 4)))
df_precise = df[(df['s_fitted']>df['s_true']-0.1) & (df['s_fitted']< df['s_true']+0.1)]
df_s_equal_1 = (df [ (df['s_fitted']<1.05) & (df['s_fitted']>0.95) ][ (df['s_fitted'] > df['s_true']+0.25) | (df['s_fitted'] < df['s_true']-0.25) ]
                   [ (df['s_fitted'] < s_deg - 0.2) | (df['s_fitted'] > s_deg +0.2) ]).reset_index(drop=True)
s_deg_ = (df[ (df['s_fitted'] < s_deg+0.2 ) & (df['s_fitted'] > s_deg-0.2 )] [ (df['s_fitted']>1.05) | (df['s_fitted']<0.95) ]
            [ (df['s_fitted'] > df['s_true']+0.25) | (df['s_fitted'] < df['s_true']-0.25) ]).reset_index(drop=True)
            
s_err_up = (df [ (df['s_fitted'] < df['s_true']+0.25 ) & (df['s_fitted'] > df['s_true']+0.1 )]).reset_index(drop=True)

s_err_down = (df [ (df['s_fitted'] < df['s_true']-0.1 ) & (df['s_fitted'] > df['s_true']-0.25 )]).reset_index(drop=True)
#[ (df['s_fitted']>1.05) | (df['s_fitted']<0.95) ][ (df['s_fitted'] > s_deg+0.2) | (df['s_fitted'] < s_deg-0.2) ]

s_scatter_1 = (df[ (df['s_fitted']> s_deg+0.2) & (df['s_fitted']> df['s_true']+0.25) ]).reset_index(drop=True)
s_scatter_2 = (df[ (df['s_fitted']> s_deg+0.2) & (df['s_fitted']< df['s_true']-0.25) ][(df['s_fitted']>1.05) | (df['s_fitted']<0.95)]).reset_index(drop=True)
s_scatter_3 = (df[ (df['s_fitted']< s_deg-0.2) & (df['s_fitted']< df['s_true']-0.25) ]).reset_index(drop=True)
s_scatter_4 = (df[ (df['s_fitted']< s_deg-0.2) & (df['s_fitted']> df['s_true']+0.25) ][(df['s_fitted']>1.05) | (df['s_fitted']<0.95)]).reset_index(drop=True)


#df['x_c'] = np.sqrt( df['u0_fitted']**2 + np.abs(df['tp']-df['t0_fitted'])**2/(df['tE_fitted']) )

plt.close('all')
plt.figure(1)
plt.rc('axes',edgecolor='black')
plt.rcParams['axes.facecolor'] = 'white'

plt.title('Fitted'+ str(y1_name)+ 'vs True '+str(x1_name)+ 'for ' +str(len(df))+ ' targets',size=20)

plt.plot (x1,y1,'b.',label='_nolegend_',markersize=8 , alpha = 0.4)

#plt.plot (x1,x1,'g-',label='_nolegend_')

#plt.text(max(df['s_true'])-1.5, (high_cut-low_cut)/2, 'Median Absolute Deviation = '+str(round(med_med(df['s_true'],df['s_fitted']),3)),size=15)
#plt.axis([0,5,0,5])

#plt.legend()
plt.xlabel(str(x1_name),size=20)
plt.ylabel(str(y1_name),size=20)
plt.grid()


plt.figure(2) 
plt.title('Fitted'+ str(y1_name)+ 'vs True '+str(x1_name)+ 'for ' +str(len(df))+ ' targets',size=20)
plt.plot (df['q_true_log'],df['q_fitted_log'],'b.',markersize=10,label='_nolegend_',alpha=0.4)
plt.plot (df['q_true_log'],df['q_true_log'],'g-',label='_nolegend_')
#plt.plot (df['q_true_log'],df['q_true_log']+1,'y-.',label = '_nolegend_')
#plt.plot (df['q_true_log'],df['q_true_log']-1,'y-.',label = '_nolegend_')
#plt.plot (df['q_true_log'],df['q_true_log']-2.5,'y--', label = '_nolegend_')

#plt.plot (df['q_true_log'][df['s_fitted']<low_cut],df['q_fitted_log'][df['s_fitted']<low_cut],'.',color='red',markeredgecolor='none',markersize=10, label = '_nolegend_')
#plt.plot (df['q_true_log'][df['s_fitted']>high_cut],df['q_fitted_log'][df['s_fitted']>high_cut],'.',color='red',markeredgecolor='none',markersize=10,label = 'Failure')
#plt.plot (df['q_true_log'][(df['s_fitted']<high_cut) & ( df['s_fitted']> df['s_true']+err_1 ) ],df['q_fitted_log'][(df['s_fitted']<high_cut) & ( df['s_fitted']> df['s_true']+err_1 ) ],'.',color = 'orange',markersize=10, label = 's_true + 0.25 < s_fitted < 2.5')
#plt.plot (df['q_true_log'][(df['s_fitted']<df['s_true']-err_1) & ( df['s_fitted']> low_cut ) ],df['q_fitted_log'][(df['s_fitted']<df['s_true']-err_1) & ( df['s_fitted']> low_cut ) ],'.',color = 'orange',markersize=10, label = '0.28 < s_fitted < s_true - 0.25')
plt.plot (df['q_true'],df['q_true']+err_2,'g--')
#plt.plot (df['q_true_log'][(df['s_fitted']<df['s_true']+err_1) & ( df['s_fitted']> df['s_true']+err_2 ) ],df['q_fitted_log'][(df['s_fitted']<df['s_true']+err_1) & ( df['s_fitted']>df['s_true']+err_2) ],'.',color = '#6BFF33',markersize=10, label = 's_true + 0.1 < s_fitted < s_true + 0.25')
#plt.plot (df['q_true_log'][(df['s_fitted']<df['s_true']-err_2) & ( df['s_fitted']> df['s_true']-err_1 ) ],df['q_fitted_log'][(df['s_fitted']<df['s_true']-err_2) & ( df['s_fitted']>df['s_true']-err_1) ],'.',color = '#6BFF33',markersize=10, label = 's_true - 0.25 < s_fitted < s_true - 0.1')
#plt.text(min(df['q_fitted_log'])+1, (min(df['q_fitted_log'])-max(df['q_fitted_log']))/2, 'Median Absolute Deviation = '+str(round(med_med(df['q_true_log'],df['q_fitted_log']),3)),size=15)

plt.axis([min(df['q_fitted_log']), max(df['q_fitted_log']), min(df['q_fitted_log']), max(df['q_fitted_log'])])

plt.xlabel(str(x1_name),size=15)
plt.ylabel(str(y1_name),size=15)
#plt.legend()
plt.grid()


'''
plt.figure(3)
plt.rc('axes',edgecolor='black')
plt.rcParams['axes.facecolor'] = 'white' 
plt.title('Fitted Pojected Separation vs True Projected Separation for '+str(len(df_new_del_t))+' targets',size=15)
plt.plot (df_new_del_t['s_true'],df_new_del_t['s_fitted'],'b.',label='_nolegend_',markersize=8 , alpha = 0.4)
plt.plot (df_new_del_t['s_true'],df_new_del_t['s_true'],'g-',label='_nolegend_')
plt.axis([0,5,0,5])
plt.grid()


plt.figure(6)
plt.rc('axes',edgecolor='black')
plt.rcParams['axes.facecolor'] = 'white'
plt.title('Fitted Pojected Separation vs True Projected Separation for '+str(len(df_new_del_t))+' targets',size=15)
plt.plot (df_new_del_t['q_true_log'],df_new_del_t['q_fitted_log'],'b.',label='_nolegend_',markersize=8 , alpha = 0.4)
plt.plot (df_new_del_t['q_true_log'],df_new_del_t['q_true_log'],'g-',label='_nolegend_')
plt.axis([min(df_new_del_t['q_fitted_log']), max(df_new_del_t['q_fitted_log']), min(df_new_del_t['q_fitted_log']), max(df_new_del_t['q_fitted_log'])])
plt.grid()



plt.figure(4)
plt.rc('axes',edgecolor='black')
plt.rcParams['axes.facecolor'] = 'white'
plt.title('Fitted Pojected Separation vs True Projected Separation for df_new_ampl',size=15)
plt.plot (df_new_ampl['s_true'],df_new_ampl['s_fitted'],'b.',label='_nolegend_',markersize=8 , alpha = 0.4)
plt.plot (df_new_ampl['s_true'],df_new_ampl['s_true'],'g-',label='_nolegend_')
plt.axis([0,5,0,5])


plt.figure(5)
plt.rc('axes',edgecolor='black')
plt.rcParams['axes.facecolor'] = 'white'
plt.title('Fitted Pojected Separation vs True Projected Separation for df_new_ampl',size=15)
plt.plot (df_new_ampl['q_true_log'],df_new_ampl['q_fitted_log'],'b.',label='_nolegend_',markersize=8 , alpha = 0.4)
plt.plot (df_new_ampl['q_true_log'],df_new_ampl['q_true_log'],'g-',label='_nolegend_')
plt.axis([min(df_new_ampl['q_fitted_log']), max(df_new_ampl['q_fitted_log']), min(df_new_ampl['q_fitted_log']), max(df_new_ampl['q_fitted_log'])])


plt.figure(7)
plt.rc('axes',edgecolor='black')
plt.rcParams['axes.facecolor'] = 'white'
plt.title('Fitted Pojected Separation vs True Projected Separation for df_new_ampl_del_t',size=10)
plt.plot (df_new_ampl_del_t['s_true'],df_new_ampl_del_t['s_fitted'],'b.',label='_nolegend_',markersize=8 , alpha = 0.4)
plt.plot (df_new_ampl_del_t['s_true'],df_new_ampl_del_t['s_true'],'g-',label='_nolegend_')
plt.axis([0,5,0,5])

plt.figure(8)
plt.rc('axes',edgecolor='black')
plt.rcParams['axes.facecolor'] = 'white'
plt.title('Fitted Pojected Separation vs True Projected Separation for df_new_ampl_del_t',size=10)
plt.plot (df_new_ampl_del_t['q_true_log'],df_new_ampl_del_t['q_fitted_log'],'b.',label='_nolegend_',markersize=8 , alpha = 0.4)
plt.plot (df_new_ampl_del_t['q_true_log'],df_new_ampl_del_t['q_true_log'],'g-',label='_nolegend_')
plt.axis([min(df_new_ampl_del_t['q_fitted_log']), max(df_new_ampl_del_t['q_fitted_log']), min(df_new_ampl_del_t['q_fitted_log']), max(df_new_ampl_del_t['q_fitted_log'])])
'''




plt.figure(3)
plt.rc('axes',edgecolor='black')
plt.rcParams['axes.facecolor'] = 'white'

plt.title('Fitted log(Pojected Separation) vs True log(Projected Separation) for ' +str(len(df))+ ' targets',size = 15)

plt.plot (df['s_true'],df['s_fitted'],'b.',label='_nolegend_',markersize=8,alpha=0.4)

plt.plot (df['s_true'],df['s_true'],'g-',label='_nolegend_')

#plt.plot (df['s_true_log'],df['s_true_log']+err_1,'y-.',label='_nolegend_')
#plt.plot (df['s_true_log'],df['s_true_log']-err_1,'y-.',label='_nolegend_')

#plt.plot (df['s_true_log'][df['s_fitted']<low_cut],df['s_fitted_log'][df['s_fitted']<low_cut],'r.',markersize=8,label = '_nolegend_')
#plt.plot (df['s_true_log'][df['s_fitted']>high_cut],df['s_fitted_log'][df['s_fitted']>high_cut],'r.',markersize=8, label = '_nolegend_')

#plt.plot (df['s_true_log'][(df['s_fitted']<high_cut) & ( df['s_fitted']> df['s_true']+err_1 ) ],df['s_fitted_log'][(df['s_fitted']<high_cut) 
#          & ( df['s_fitted']> df['s_true']+err_1 ) ],'.',color = 'orange',markeredgecolor='black',markersize=8,label = 's_true + 0.25 < s_fitted < 2.5')

#plt.plot (df['s_true_log'][(df['s_fitted']<df['s_true']-err_1) & ( df['s_fitted']> low_cut ) ],df['s_fitted_log'][(df['s_fitted']<df['s_true']-err_1) 
#          & ( df['s_fitted']> low_cut ) ],'.',color = 'orange',markeredgecolor='black', markersize=8, label = '0.28 < s_fitted < s_true - 0.25')

#plt.plot (df['s_true_log'],df['s_true_log']-err_2,'y-.',label='_nolegend_')
#plt.plot (df['s_true_log'],df['s_true_log']+err_2,'y-.',label='_nolegend_')

#plt.plot (df['s_true_log'][(df['s_fitted']<df['s_true']+err_1) & ( df['s_fitted']> df['s_true']+err_2 ) ],df['s_fitted_log'][(df['s_fitted']<df['s_true']+err_1) 
#          & ( df['s_fitted']>df['s_true']+err_2) ],'.',markersize=8,color = '#6BFF33',markeredgecolor='black', label = 's_true + 0.1 < s_fitted < s_true + 0.25')

#plt.plot (df['s_true_log'][(df['s_fitted']<df['s_true']-err_2) & ( df['s_fitted']> df['s_true']-err_1 ) ],df['s_fitted_log'][(df['s_fitted']<df['s_true']-err_2) 
#          & ( df['s_fitted']>df['s_true']-err_1) ],'.',markersize=8,color = '#6BFF33',markeredgecolor='black', label = 's_true - 0.25 < s_fitted < s_true - 0.1')

#plt.axis([min(df['s_true']), max(df['s_true']), low_cut, high_cut])

#plt.legend()
plt.xlabel('s_true_log',size = 15)
plt.ylabel('s_fitted_log', size =15)
plt.grid()







plt.figure(4)
plt.plot(df['s_true'],df['s_true'],'g-')
plt.plot(df_s_equal_1['s_true'],df_s_equal_1['s_fitted'],'k.',alpha=0.7)
plt.plot(df['s_true'],s_deg,'g.',alpha=0.4)
plt.plot(s_deg_['s_true'],s_deg_['s_fitted'],'k.',alpha=0.7)
plt.plot(s_err_up['s_true'],s_err_up['s_fitted'],'c.',alpha=0.4 )
plt.plot(s_err_down['s_true'],s_err_down['s_fitted'],'c.',alpha=0.4 )
plt.plot(s_scatter_1['s_true'],s_scatter_1['s_fitted'],'r.',alpha=0.7)
plt.plot(s_scatter_2['s_true'],s_scatter_2['s_fitted'],'r.',alpha=0.7)
plt.plot(s_scatter_3['s_true'],s_scatter_3['s_fitted'],'r.',alpha=0.7)
plt.plot(s_scatter_4['s_true'],s_scatter_4['s_fitted'],'r.',alpha=0.7)

'''
plt.figure(4)
plt.plot(df['q_fitted_log'],np.abs(df['s_fitted_log']),'b.',alpha=0.5)

plt.figure(5)
plt.plot(df['q_true_log'],np.abs(df['s_true_log']),'r.',alpha=0.5)
'''


plt.show()
