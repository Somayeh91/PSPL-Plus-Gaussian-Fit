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


home = os.path.expanduser("~")

direc = os.listdir(".")
df1 = pd.DataFrame({'f_s' : [], 'name' : [], 'q' : [], 's' : [], 't0':[], 'tE' : [], 'u0' : []})
for i in range(len(direc)):
    
    tempdata = direc[i]
    #print tempdata
    
    if tempdata.endswith('.gz'):
        fname = gzip.open(tempdata, 'rb')
        x_0 = fname.readlines()[0:7]
        f_s_true = x_0[0].split(' ')[1]
        q = x_0[5].split(' ')[5]
        s = x_0[5].split(' ')[6]
        tE_theo = x_0[6].split(' ')[4]
        t0_theo = x_0[6].split(' ')[3]
        u0_theo = x_0[6].split(' ')[1]
        
        
        par = [f_s_true,tempdata,q,s,t0_theo,tE_theo,u0_theo]
        #print par
        df1.loc[i]=[par[n] for n in range(7)]
        #print df1

df1.to_csv('440_.CSV')        
    




 