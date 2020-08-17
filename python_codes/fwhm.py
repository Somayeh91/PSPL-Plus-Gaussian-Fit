import glob,os,sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import *
from os.path import expanduser

home = os.path.expanduser("~")

### approach 1:

def fwhm3(valuelist, peakpos,base):
    """calculates the full width at half maximum (fwhm) of some curve.
    the function will return the fwhm with sub-pixel interpolation. It will start at the maximum position and 'walk' left and right until it approaches the half values.
    INPUT: 
    - valuelist: e.g. the list containing the temporal shape of a pulse 
    OPTIONAL INPUT: 
    -peakpos: position of the peak to examine (list index)
    the global maximum will be used if omitted.
    OUTPUT:
    -fwhm (value)
    """
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
    
tempdata = home+'/Library/Mobile Documents/com~apple~CloudDocs/OSU trip/Matt/cassan_32_4_1923.det.lc'

t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)

t = ma.masked_where(code != 0, t)
f = ma.masked_where(code != 0  , f)
f_err = ma.masked_where(code != 0 , f_err)
print min(f)
ind1, ind2 = fwhm3(f,f.argmax(),min(f))
print t[ind2]-t[ind1]



plt.close()
plt.plot(t,f)
plt.show()

### Approach 2:

home = os.path.expanduser("~")


name = 'cassan_0_30_255.det.lc.gz'  
tempdata = home+'/Desktop/alllc1/'+str(name)

t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)
df = pd.DataFrame({'t':t , 'f':f , 'f_err' : f_err , 'f_true': f_true, 'code':code})
df = df[df['code']==4]
df = df.reset_index(drop=True)

max_index = df['f'].argmax()
t_0 = df['t'][max_index]
f_max = max(df['f'])
base = min(df['f'])
phalf = (f_max - base)/2

for i in range (max_index,1,-1):
    #print i
    if df['f'][i] < phalf:
        ind1 = i+1
        break
for j in range (max_index,len(df['f'])-1,1):
    if df['f'][j] < phalf:
        ind2 = j-1 
        break
        
print df['t'][ind1], df['t'][ind2]       