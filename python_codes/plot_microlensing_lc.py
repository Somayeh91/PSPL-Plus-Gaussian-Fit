# -*- coding: utf-8 -*-
import glob,os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
from numpy import *
import re
import scipy.stats as st
from os.path import expanduser
import cmath
import emcee
import scipy.optimize as op
import corner
import numpy.ma as ma
import pandas as pd
import gzip

def fun (t0,u0,tE,f_s):
    u = np.sqrt(u0**2+((df['t']-t0)/tE)**2)
    A = ((u**2)+2)/(u*np.sqrt(u**2+4))
    return (f_s * (A-1)) +1


home = os.path.expanduser("~")

#Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/lcsample2

name = 'cassan_32_136_3046.det.lc.gz'
tempdata = home+'/Library/Mobile Documents/com~apple~CloudDocs/project_microlensing/shared_lc/'+str(name)

t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)
df = pd.DataFrame({'t':t , 'f':f , 'f_err' : f_err , 'f_true': f_true, 'code':code})
dff = df
df = df[df['code']==4]
df = df.reset_index(drop=True)
        
#df['f'] = df['f_true']
        
fname = gzip.open(tempdata, 'rb')
x_0 = fname.readlines()[0:7]
f_s = float(x_0[0].split(' ')[4])
q = float(x_0[5].split(' ')[5])
s = x_0[5].split(' ')[6]
tE = float(x_0[6].split(' ')[4])
t0 = float(x_0[6].split(' ')[3])
u0 = float(x_0[6].split(' ')[1])
del_t = (np.sqrt(np.abs((float(s)-(1/float(s)))**2 - (float(u0))**2))) * float(tE)
tEp = np.sqrt(q) * tE
F_t = fun(t0,u0,tE,f_s)    

f_ris_true = df['f']-F_t   

if np.abs(max(f_ris_true))> np.abs(min(f_ris_true)):
    tpp = df['t'][f_ris_true.argmax()]
else:
    tpp = df['t'][f_ris_true.argmin()]  
tp_ = np.array([t0-del_t,t0+del_t])
tp = tp_[(tp_ - tpp).argmin()] 
     
plt.close()
plt.figure(11)
plt.rc('axes',edgecolor='black')
plt.rcParams['axes.facecolor'] = 'white'
plt.title('Plot of '+str(name),size= 20)
plt.xlabel('Time',size=18)
plt.tick_params(axis='x', labelsize=13)
plt.tick_params(axis='y', labelsize=13)
#plt.xlim(325,350)
#fig = plt.axes()
#fig.xaxis.set_major_formatter(ScalarFormatter(useMathText=False,useOffset=False))
#plt.axvline(x=t0 , color='black', linestyle = '--')
#plt.axvline(x=1592.34 , color='black', linestyle = '--')
#plt.axvline(x=1592.08 , color='black', linestyle = '--')
#plt.axvline(x=t0-tE , color='black', linestyle = '--')
#plt.axhline(y=1.34 , color='black', linestyle = '--')
#plt.axhline(y=max((df['f_true']-(1-f_s))/f_s) , color='black', linestyle = '--')



plt.ylabel('Magnification',size=18)
plt.plot(df['t'],(df['f_true']),'k.',linewidth=6)
#plt.axis([1586,1594, 2.8, 4.94])
#plt.grid()

plt.show()
#plt.plot(t,f_model)
'''
fig = plt.gcf()
fig.set_size_inches(11,8.5)
fig.savefig(home+'/Desktop/pres3.png',facecolor='white')
'''