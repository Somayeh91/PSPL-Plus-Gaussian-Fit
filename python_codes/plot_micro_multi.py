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


home = os.path.expanduser("~")

direc = os.listdir(".")

for i in direc:
    tempdata = i
    if tempdata.endswith('.gz'):
        t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)

        t = ma.masked_where(code != 0, t)
        f = ma.masked_where(code != 0  , f)
        f_err = ma.masked_where(code != 0 , f_err)
        f_true = ma.masked_where(code != 0 , f_true)
        fname = gzip.open(tempdata, 'rb')
        x_0 = fname.readlines()[0:7]
        t_0 = x_0[6].split(' ')[3]
        t_p = x_0[0].split(' ')[1]
        
        tE_p = x_0[5].split(' ')[6]
        tE = x_0[6].split(' ')[4]


#u = np.sqrt(0.15417495486**2+((t-538.763491798)/276.540512274)**2)
#A = ((u**2)+2)/(u*np.sqrt(u**2+4))
#f_model = (0.987150671112 * (A-1)) +1
#print t
        fig = plt.gcf()        
        plt.close()
        plt.plot(t,f_true,'b')
#plt.plot(t,f_model)
        fig.savefig(home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/plots/'+str(i)+'.png')