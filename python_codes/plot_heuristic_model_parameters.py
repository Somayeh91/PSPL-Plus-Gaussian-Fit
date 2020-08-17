# -*- coding: utf-8 -*-
import glob,os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
from os.path import expanduser
import numpy.ma as ma
import pandas as pd
import gzip

def fun (t0,u0,tE,f_s):
    u = np.sqrt(u0**2+((df['t']-t0)/tE)**2)
    A = ((u**2)+2)/(u*np.sqrt(u**2+4))
    return (f_s * (A-1)) +1


home = os.path.expanduser("~")

#Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/lcsample2

name = 'cassan_18_123_294.det.lc.gz'
tempdata = home+'/Desktop/selected plots for paper/'+str(name)

t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)
df = pd.DataFrame({'t':t , 'f':f , 'f_err' : f_err , 'f_true': f_true, 'code':code})
df = df[df['code']==0]
        
#df['f'] = df['f_true']
        
fname = gzip.open(tempdata, 'rb')
x_0 = fname.readlines()[0:7]
f_s = float(x_0[0].split(' ')[1])
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
#plt.title('Plot of '+str(name),size= 20)
plt.xlabel('Time',size=23)
#plt.tick_params(axis='x', labelsize=20)
#plt.tick_params(axis='y', labelsize=20)
plt.tick_params(axis='x', which='both', bottom=False, top=False,labelbottom=False)
plt.tick_params(axis='y', which='both', left=False, right=False,labelleft=False)

plt.axvline(x=349.51,color='black', linestyle='dashed' , ymax = 0.958 )
plt.axvline(x=352.89,color='black', linestyle='dashed' , ymax = 0.75 )
plt.axhline(y=3.89,color='black', linestyle='dashed', xmin = 0.1128, xmax= 0.5  )
plt.axhline(y=13.8527,color='black', linestyle='dashed', xmin = 0.13, xmax= 0.5  )
plt.axhline(y=10.3303,color='black', linestyle='dashed', xmin = 0.6026, xmax= 0.616827  )

plt.text(341.081,4.10473,'$t_E$',size=20)
plt.text(335.353,13.702, '$A_{Max}$',size=20)
plt.text(349.675,2.54145, '$t_0$',size=20)
plt.text(353.125,2.54145, '$t_p$',size=20)
plt.text(354.119,10.5441, '$a$',size=20)
plt.text(350.494,9.4456, '$t_{Ep}$',size=20)

plt.annotate(s='', xy=(353.768,9.8536), xytext=(353.768,11.1403),arrowprops=dict(arrowstyle='<->'))
plt.annotate(s='', xy=(352.599,10.1989), xytext=(351.722,9.7594),arrowprops=dict(arrowstyle='->'))

plt.ylabel('Magnification',size=23)
plt.plot(df['t'],(df['f_true']-(1-f_s))/f_s,'k-',linewidth = 2.0, color = 'blue')
plt.axis([335,364, 2.259, 14.31])
#plt.grid()

plt.show()
#plt.plot(t,f_model)

fig = plt.gcf()
fig.set_size_inches(10,8)
fig.savefig(home+'/Desktop/pres3.png',facecolor='white')
