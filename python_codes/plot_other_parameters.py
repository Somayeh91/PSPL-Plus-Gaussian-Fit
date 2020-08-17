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

name = 'alllc_result.CSV'  
#Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/result_file/
tempdata = home+'/Desktop/trial_runs/'+str(name)

df = pd.read_csv(tempdata)
df['u0_true'] = np.abs(df['u0_true'])
df['u0_fitted'] = np.abs(df['u0_fitted'])



df['q_fitted_log'][df['q_fitted_log']>0] = np.log10( 1/df['q_fitted'][df['q_fitted_log']>0] )
df['q_fitted'][df['q_fitted_log']>0] = ( 1/df['q_fitted'][df['q_fitted_log']>0] )

u = df['u0_true']
df['max_A'] = (2+u**2)/(u*np.sqrt(4+u**2))

df = df[df['s_fitted']<10][np.abs(df['u0_fitted'])>0.045] #02/20/2018


x1 = df['t0_fitted']
y1 = df['u0_fitted']
x1_name =' Ratio of the Fitted Separation Over True Separation' 
y1_name = 'Amplitude of the Planetary Event'
title = 'Amplitude of the Planetary Event Versus Ratio of the Fitted Separation Over True Separation'

x2 = df['t0_true']
y2 = df['u0_true']

plt.close('all')
plt.figure(1)
plt.rc('axes',edgecolor='black')
plt.rcParams['axes.facecolor'] = 'white'

plt.title(str(title),size=20)

plt.plot (x1,y1,'b.',label='_nolegend_',markersize=8 , alpha = 0.4)
plt.plot (x2,y2,'g.',label='_nolegend_',markersize=8 , alpha = 0.4)

#plt.axis([0,200,0,50])

#plt.legend()
plt.xlabel(str(x1_name),size=15)
plt.ylabel(str(y1_name),size=15)
plt.grid()

f = plt.gcf()

f.set_size_inches(15.0,7)
f.savefig(home+'/Desktop/'+ str(title)+'.png')
plt.show()
