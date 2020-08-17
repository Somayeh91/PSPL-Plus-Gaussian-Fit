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

start = time.time()

home = os.path.expanduser("~")

direc = os.listdir(".")

name1 = 'allc1_info.CSV'  
tempdata1 = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/allc/'+str(name1)
df1 = pd.read_csv(tempdata1)
temp_time = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/allc/time_stamps.CSV'

'''
name2 = 'allc2_info.CSV'  
tempdata2 = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/allc/'+str(name1)




df2 = pd.read_csv(tempdata2)



alllc_ = [df1,df2]
alllc = pd.concat(alllc_)
alllc.to_csv(home+'/Desktop/alllc_info.CSV')
'''

time = pd.read_csv(temp_time)

A_cut = 1

t = time['t'] 
x1= []
x_p_1 , x_p_2 = [], []

df1['del_t'] = (np.sqrt(np.abs((df1['s']-(1/df1['s']))**2 - (df1['u0'])**2))) * df1['tE']
df1['p_range_left'] = df1['t0']-df1['del_t']
df1['p_range_right'] = df1['t0']+df1['del_t']
'''
for i in range (len(df1['t0'])):
    x1.append (min(np.abs(df1['t0'][i]-t)))
    x_p_1.append (min(np.abs(df1['p_range_left'][i]-t)))
    x_p_2.append (min(np.abs(df1['p_range_right'][i]-t)))
    
  
x1 = asarray(x1)

A = df1['u0']**2 + 2 / ( np.abs(df1['u0']) * np.sqrt(df1['u0']**2 + 4))

df1['A_'] = df1['f_s']*(A) + (1-df1['f_s'])
df1['coverage'] =   df1['tE'] - x1
df1['coverage_p_left'] =   (np.sqrt(df1['q'])*df1['tE']) - x_p_1
df1['coverage_p_right'] =   (np.sqrt(df1['q'])*df1['tE']) - x_p_2



df1_new = df1[( (df1['s']<5)) & ( (df1['coverage_p_left']>0) & (df1['coverage_p_right']>0) )].reset_index(drop=True)
 
df1_excluded = []

for l in range(len(df1)):
    cc = 0
    for ll in range(len(df1_new)):
        if df1['name'][l] == df1_new['name'][ll]:
            cc = 1
    if cc == 0:
        df1_excluded.append(df1['name'][l])

#df1_new.to_csv(home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/allc/allc2_A_larger_than_'+str(A_cut)+'.CSV')        
df1_excluded = np.asarray(df1_excluded)
'''

print (len(df1))
#print (len(df1_new))


#plt.close()
#plt.figure(10)
#plt.plot(df1['s'],df1['A_'],'r.',df1_new['s'],df1_new['A_'],'k.')
#plt.xlabel('q')
#plt.ylabel('A')
#plt.plot(df1['s'][df1['A_']>2],df1['A_'][df1['A_']>2],'b.',df1['s'][df1['A_']<2],df1['A_'][df1['A_']<2],'r.')
#plt.plot(df2['s'][df2['A_']>2],df2['A_'][df2['A_']>2],'b.',df2['s'][df2['A_']<2],df2['A_'][df2['A_']<2],'r.')
#plt.show()




