import glob,os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
from numpy import *
import re
import scipy.stats as st
from os.path import expanduser
import cmath
import emcee
import scipy.optimize as op
import corner
import time
import gzip

start = time.time()

home = os.path.expanduser("~")

direc = os.listdir(".")


for i in direc:

    name = str(i)  
    #tempdata = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/lcsample/'+str(name)
    tempdata = home+'/Desktop/alllc1/'+str(name)
    
    if tempdata.endswith('.gz'):
        t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)

        t = ma.masked_where(code != 0, t)
        f = ma.masked_where(code != 0  , f)
        f_true = ma.masked_where(code != 0 , f_true)
        fname = gzip.open(str(tempdata), 'rb')
        f_s_true = fname.readline().split(' ')[1]
        x = fname.readlines()[4].split(' ')
        y = gzip.open(str(tempdata), 'rb').readlines()[6].split(' ')
        q = x[5]
        s = x[6]
        
        t0 = y[3]
        tE = y[4]

#print "mass_ratio_true= "+q,"separation_true= "+s
#print "t0= "+t0, "tE= "+tE
#print "theo t0-tp = "+str(float(s)*float(tE)),"theo tEp = "+str(np.sqrt(float(q))*float(tE))


        plt.close()
        plt.title(str(name))
        plt.plot(t,f_true,'k')
        plt.axvline(x=float(t0)+(float(s)*float(tE)), color='red',linewidth=2)
        plt.axvline(x=float(t0)-(float(s)*float(tE)), color='red',linewidth=2)
        plt.axis([float(t0)-(float(s)*float(tE))-5,float(t0)+(float(s)*float(tE))+5, min(f_true), max(f_true)])
        fig = plt.gcf()
        fig.set_size_inches(12.0,6.0)
        fig.savefig('./../lc_check/'+str(i.split('.det')[0])+'.png')
        #plt.plot(t,f_true-fun(t0_ml, u0_ml, tE_ml,f_s_ml),'b')
        end = time.time()
        #print end-start
