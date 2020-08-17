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
    file_content = []
    if i.endswith('.gz'):
        fname = gzip.open(str(i), 'rb')
        print fname
        f_s_true = fname.readline().split(' ')[1]
        t,f,f_err,f_true,code = np.loadtxt(str(i),usecols=(0,1,2,3,5),unpack=True)
        print t