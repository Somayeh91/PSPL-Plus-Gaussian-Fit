import glob,os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from numpy import *
import re
import scipy.stats as st
from os.path import expanduser
import cmath
import scipy.optimize as op
import time
import gzip
import pandas as pd

# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "serif"

def med_med (true,fitted):
    temp = fitted - true
    return (np.median(np.abs(temp-np.median(temp))))


start = time.time()

home = os.path.expanduser("~")

direc = os.listdir(".")
'''
name = '440_.CSV'  
tempdata = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/allc/'+str(name)
'''

name = 'alllc_result_v2.CSV'  
#Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/result_file/
tempdata = home+'/Desktop/trial_runs/'+str(name)
#&((df['f_s_true']*( (2 + df['u0_true']**2) / (df['u0_true']*np.sqrt(4 + df['u0_true']**2)) ) + (1-df['f_s_true']) )>1.5) & (df['u0_true']>0.1)
df1 = pd.read_csv(tempdata)

df1['u0_true'] = np.abs(df1['u0_true'])
df1['u0_fitted'] = np.abs(df1['u0_fitted'])

name = 'alllc_fails_v2.CSV'  
#Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/result_file/
tempdata2 = home+'/Desktop/trial_runs/'+str(name)
#&((df['f_s_true']*( (2 + df['u0_true']**2) / (df['u0_true']*np.sqrt(4 + df['u0_true']**2)) ) + (1-df['f_s_true']) )>1.5) & (df['u0_true']>0.1)
df2 = pd.read_csv(tempdata2)

df2['u0_true'] = np.abs(df2['u0_true'])
df2['u0_fitted'] = np.abs(df2['u0_fitted'])

#df = df[ (df['chi_2_2']>-25000) & (df['s_fitted']<5) & ((df['f_s_true']*( (2 + df['u0_true']**2) / (df['u0_true']*np.sqrt(4 + df['u0_true']**2)) ) + (1-df['f_s_true']) )>1.1)]
#df = df[(np.abs(df['t0_fitted']-df['tp_fitted'])>1) & ((-2*df['chi_2_2']/41039)>1.003) ]

#df['q_fitted_log'][df['q_fitted_log']>0] = np.log10( 1/df['q_fitted'][df['q_fitted_log']>0] )
#df['q_fitted'][df['q_fitted_log']>0] = ( 1/df['q_fitted'][df['q_fitted_log']>0] )

#df = df[df['s_fitted']<5][(np.abs(df['t0_fitted']-df['tp_fitted'])>1) & ((-2*df['chi_2_2']/41039)>1.003) ]
#df = df[df['s_fitted']<10][np.abs(df['u0_fitted'])>0.045]


frames = [df1, df2]
df = pd.concat(frames)

print(len(df))

color = ['#e41a1c', '#f781bf', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628','#e7298a', '#e6ab02']
temp_time = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/allc/time_stamps.CSV'
time_ = pd.read_csv(temp_time)

q_range =  np.array([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,1])
s_range =  np.array([0,0.4,0.8, 0.9 , 1, 1.1, 1.2 , 1.5, 2.5, 3.5, 5])
p_range =  [0.5,5,20,30,40,50,500]


plt.close()

'''
plt.figure(1)
plt.title('Mass Ratio',size=20)
plt.xlabel('True Mass Ratio - log',size=20)
plt.ylabel('Fitted Mass Ratio - log',size=20)
plt.plot(df['q_true_log'],df['q_true_log'],'g-',label='_nolegend_')
for i in range(len(s_range)-1):
    plt.plot(df['q_true_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],df['q_fitted_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],'o',color=color[i],alpha=0.8,label= '('+str((s_range[i]))+')-('+str((s_range[i+1]))+')')
    plt.axis([min(df['q_fitted_log']),max(df['q_fitted_log']),min(df['q_fitted_log']),max(df['q_fitted_log'])])
plt.legend(title='Separation')   
plt.figure(2)
plt.title('Projected Separation',size=20)
plt.xlabel('True Separation ',size=20)
plt.ylabel('Fitted Separation ',size=20) 
plt.plot(df['s_true'],df['s_true'],'g-',label='_nolegend_')
for i in range(len(q_range)-1):
       plt.plot(df['s_true'][(df['q_true']>=q_range[i])&(df['q_true']<q_range[i+1])],df['s_fitted'][(df['q_true']>=q_range[i])&(df['q_true']<q_range[i+1])],'o',color=color[i],alpha=0.8,label= '('+str((q_range[i]))+')-('+str((q_range[i+1]))+')')
plt.axis([min(df['s_true']),max(df['s_true']),min(df['s_true']),max(df['s_true'])])
plt.legend(title='Mass Ratio')   
   
plt.legend(title='Mass Ratio') 


i = 2
j= 2

plt.figure(1)
plt.title('Mass Ratio',size=25)
plt.xlabel('Log (True Mass Ratio)',size=25)
plt.ylabel('Log (Fitted Mass Ratio)',size=25)
plt.plot(df['q_true_log'],df['q_true_log'],'g-',label='_nolegend_')
plt.plot(df['q_true_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],df['q_fitted_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],'o',color=color[i],alpha=0.8,label= '('+str((s_range[i]))+')-('+str((s_range[i+1]))+')')
#plt.axis([min(df['q_fitted_log']),max(df['q_fitted_log']),min(df['q_fitted_log']),max(df['q_fitted_log'])])
plt.axis([-8,0,-8,0])
plt.text(-7.5,-0.5, 'Median Absolute Deviation = '+str(round(med_med(df['q_true_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])],df['q_fitted_log'][(df['s_true']>s_range[i])&(df['s_true']<s_range[i+1])]),6)))
plt.legend(title='Separation')   
plt.grid()
#fig = plt.gcf()
#fig.savefig(home+'/Desktop/'+'('+str((s_range[i]))+')-('+str((s_range[i+1]))+').png')

plt.figure(2)
plt.title('Projected Separation',size=25)
plt.xlabel('True Separation ',size=25)
plt.ylabel('Fitted Separation ',size=25) 
plt.plot(df['s_true'],df['s_true'],'g-',label='_nolegend_')
plt.plot(df['s_true'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])],df['s_fitted'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])],'o',color=color[j],alpha=0.8,label= '('+str((q_range[j]))+')-('+str((q_range[j+1]))+')')
#plt.axis([min(df['s_true']),max(df['s_true']),min(df['s_true']),max(df['s_true'])])
plt.axis([0,5,0,5])
plt.text(4,0.5, 'Median Absolute Deviation = '+str(round(med_med(df['s_true'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])],df['s_fitted'][(df['q_true']>=q_range[j])&(df['q_true']<q_range[j+1])]),3)))
plt.legend(title='Mass Ratio')   
plt.grid()   
#fig = plt.gcf()
#fig.savefig(home+'/Desktop/'+'('+str((q_range[j]))+')-('+str((q_range[j+1]))+').png')


      
plt.figure(6)
plt.plot(df_new_del_t['s_true'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5],df_new_del_t['u0_fitted'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5]/df_new_del_t['u0_true'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5],'b.',label='u0')
plt.plot(df_new_del_t['s_true'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5],df_new_del_t['t0_fitted'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5]/df_new_del_t['t0_true'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5],'k.',label='t0')
plt.plot(df_new_del_t['s_true'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5],df_new_del_t['tE_fitted'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5]/df_new_del_t['tE_true'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5],'g.',label='tE')
plt.plot(df_new_del_t['s_true'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5],np.abs(df_new_del_t['tp_fitted'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5]-df_new_del_t['t0_fitted'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5])/df_new_del_t['del_t'][df_new_del_t['s_fitted']<df_new_del_t['s_true']-0.5],'r.',label='del_t')
plt.title('df where fitted_s < true_s - 0.5')
plt.legend()

'''
red_chi_2_1 = np.abs(2*df['chi_2_1'])/38316
red_chi_2_2 = np.abs(2*df['chi_2_2'])/38316

red_chi_2_1_f = np.abs(2*df2['chi_2_1'])/38316
red_chi_2_2_f = np.abs(2*df2['chi_2_2'])/38316

bin_num = 1000

plt.close()
plt.figure(1)

plt.hist(np.abs(red_chi_2_1-red_chi_2_2),bins=np.logspace(np.log10(0.0001),np.log10(45.0),50)) 
plt.axvline(x=0.001 , color='black', linestyle = '--')

#plt.hist(np.abs(red_chi_2_1-red_chi_2_2)[np.abs(red_chi_2_1-red_chi_2_2)<0.001],bins=np.logspace(np.log10(0.0001),np.log10(45.0),50),color='#ff7f00') 
#plt.hist(np.abs(red_chi_2_1-red_chi_2_2)[red_chi_2_2>1.25],bins=np.logspace(np.log10(0.0001),np.log10(45.0),50),alpha=0.5) 
#plt.hist(np.abs(red_chi_2_2/red_chi_2_1)[np.abs(df['u0_true'])<0.05],bins=np.logspace(np.log10(0.0001),np.log10(45.0),50),alpha=0.5) 


#plt.title('Histogram of '+'${{\chi}_1}^2 - {{\chi}_2}^2$')
plt.xscale('log')
plt.yscale('log')
plt.tick_params(axis='y',labelsize=20)
plt.tick_params(axis='x',labelsize=20)

plt.ylabel('Counts',size=25)
plt.xlabel('Reduced ('+'${{\chi}_1}^2 - {{\chi}_2}^2$'+')',size=25)

plt.figure(2)
plt.hist(np.abs(red_chi_2_2/red_chi_2_1),bins=np.logspace(np.log10(0.01),np.log10(1.01),50))
plt.hist(np.abs(red_chi_2_2/red_chi_2_1)[np.abs(red_chi_2_1-red_chi_2_2)<0.005],bins=np.logspace(np.log10(0.0001),np.log10(45.0),50),alpha=0.5) 
#plt.hist(np.abs(red_chi_2_2/red_chi_2_1)[np.abs(df['u0_true'])<0.05],bins=np.logspace(np.log10(0.0001),np.log10(45.0),50),alpha=0.5) 

#plt.hist(np.abs(red_chi_2_2_f/red_chi_2_1_f),bins=np.logspace(np.log10(0.1),np.log10(1.0),50),alpha=0.5)
plt.title('Histogram of '+'${{\chi}_2}^2 / {{\chi}_1}^2$')  
plt.xscale('log')
plt.yscale('log')

#plt.hist(np.abs(red_chi_2_1-red_chi_2_2),bins=1000)  
#plt.hist(red_chi_2_2,bins=np.linspace(0,22,bin_num),alpha=0.5)  
plt.ylabel('Counts')
#plt.title('Histogram of '+'$|{{\chi}_2}^2 - {{\chi}_1}^2|$')

plt.figure(3)

plt.hist((red_chi_2_1)[df['u0_true']<0.05],bins=np.logspace(np.log10(0.8),np.log10(31300),50)) 
plt.hist((red_chi_2_1),bins=np.logspace(np.log10(0.8),np.log10(31300),50)) 

#plt.hist(np.abs(red_chi_2_1_f-red_chi_2_2_f),bins=np.logspace(np.log10(0.0001),np.log10(45.0),50),alpha=0.5) 

plt.title('Histogram of reduced '+'${{\chi}_1}^2 $')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Counts')

plt.figure(4)

plt.hist((red_chi_2_2),bins=np.logspace(np.log10(0.8),np.log10(21800),50)) 
plt.hist((red_chi_2_2)[df['u0_true']<0.05],bins=np.logspace(np.log10(0.8),np.log10(21800),50)) 

#plt.hist(np.abs(red_chi_2_1_f-red_chi_2_2_f),bins=np.logspace(np.log10(0.0001),np.log10(45.0),50),alpha=0.5) 

plt.title('Histogram of reduced '+'${{\chi}_2}^2 $')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Counts')

'''
plt.figure(2)
plt.plot(red_chi_2_1, red_chi_2_2, '.',red_chi_2_1, red_chi_2_1, 'g-')
plt.xlabel('${{\chi}_1}^2$')
plt.ylabel('${{\chi}_2}^2$')
'''

end = time.time()
print end-start
plt.show()