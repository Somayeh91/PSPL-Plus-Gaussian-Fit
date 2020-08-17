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

def empty(df):
    return len(df.index) == 0
    
def fwhm(valuelist, peakpos,base):
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

def tE_finder (t,f,f_s):
    df = pd.DataFrame({'t' : t, 'f' : f})
    
    A_base = (float(f_s)*0.34)+1
    
    t_max = df['t'][df['f'].argmax()]
    t_right = df['t'][df['t']>t_max]
    t_left = df['t'][df['t']<t_max]
    
    if empty(t_right) == 'False':
        tE_right = np.abs(t_right[(np.abs(df['f'][df['t']>t_max]-(A_base))).argmin()]-t_max)
    else:
        tE_right = 1
    if empty(t_left) == 'False':
        tE_left = np.abs(t_left[(np.abs(df['f'][df['t']<t_max]-(A_base))).argmin()]-t_max)
    else:
        tE_left = 1
    
    return min([tE_right,tE_left])

def fun (t0,u0,tE,f_s):
    u = np.sqrt(u0**2+((df['t']-t0)/tE)**2)
    A = ((u**2)+2)/(u*np.sqrt(u**2+4))
    return (f_s * (A-1)) +1
        
def fun2 (mean, sigma,amp, t0,u0,tE,f_s):
    u = np.sqrt(u0**2+((df['t']-t0)/tE)**2)
    A = (((amp/np.sqrt(2*pi*(sigma**2)))*np.exp(-((df['t']-mean)**2)/(2*(sigma**2)))))+((u**2)+2)/(u*np.sqrt((u**2)+4))
    return (f_s * (A-1)) +1

def lnlike(theta, t, f, f_err):
    t0, u0, tE,f_s = theta
    model = fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    return -0.5*(np.sum((f-model)**2*inv_sigma2))

def lnlike2(theta, t, f, f_err):
    mean, sigma,amp, t0,u0,tE,f_s = theta
    model = fun2(mean, sigma,amp, t0,u0,tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    return -0.5*(np.sum((f-model)**2*inv_sigma2))

def chisq(theta, t, f, f_err):
    t0, u0, tE,f_s = theta
    model = fun(t0, u0, tE,f_s)
    inv_sigma2 = 1.0/(f_err**2)
    chii = ((f-model)**2*inv_sigma2)
    
    c = (np.cumsum(chii))
    return c
def subtract(chi):
    c=0
    r = 0
    for i in range (len(chi)-1):
        if ((chi[i+1]-chi[i])) > c:
            c = chi[i+1]-chi[i]
            r = i
        
    return chi[r+1]
    
def slope(n):
    return float(n[3]-n[1])/float(n[2]-n[0]), n[1]-(float(n[3]-n[1])/float(n[2]-n[0]))*n[0]
    
def lin_fit (line_,y_):
    return (y_ - line_[1])/line_[0]
    
def x_y_finder (t,f,y):

    t_chi = pd.DataFrame({'t' : t, 'chi2' : f})
    
    x_m = (t_chi[t_chi['chi2']>y]).reset_index()
    if x_m.empty :
        x2,y2 = max(t_chi['t']),1
    else:
        x2 , y2 =   x_m['t'][0], x_m ['chi2'][0]
        
    x_l = (t_chi[(t_chi['chi2']<y)&(t_chi['t']<x2)]).reset_index()
    if x_l.empty :
        x1,y1 = min(t_chi['t']),0
    else:
        x1 , y1 =   x_l ['t'][max(x_l.index)], x_l ['chi2'][max(x_l.index)]

    return x1,y1,x2,y2
        
#s_ = [[],[]]
#q_ = [[],[]]

name = 'cassan_0_115_2759.det.lc.gz'  
tempdata = home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/lcsample2/'+str(name)

t,f,f_err,f_true,code = np.loadtxt(tempdata,usecols=(0,1,2,3,5),unpack=True)
df = pd.DataFrame({'t':t , 'f':f , 'f_err' : f_err , 'f_true': f_true, 'code':code})
df = df[df['code']==0]
fname = gzip.open(tempdata, 'rb')
x_0 = fname.readlines()[0:7]
f_s_true = x_0[0].split(' ')[1]
q = x_0[5].split(' ')[5]
s = x_0[5].split(' ')[6]
tE_theo = x_0[6].split(' ')[4]
t0_theo = x_0[6].split(' ')[3]
u0_theo = x_0[6].split(' ')[1]
        
    
u0_true = float(f_s_true)/(max(df['f'])-1+float(f_s_true))
t0_true =  df['t'][df['f'].argmax()] #float(t0_theo)
ind1, ind2 = fwhm(df['f'],df['f'].argmax(),min(df['f']))
#tE_true = t[ind2]-t[ind1]
tE_true = [tE_finder (df['t'],df[ 'f'], f_s_true),t[ind2]-t[ind1]]
#print 'tE_true = '+ str( tE_true)

tE_ = [[],[]]    
for i in tE_true:
               
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [t0_true, u0_true, i,f_s_true], args=(df['t'],df[ 'f'], df['f_err']),method = 'Nelder-Mead')
    t0_ml, u0_ml, tE_ml,f_s_ml = result['x']
    tE_[0].append(lnlike([t0_ml, u0_ml, tE_ml,f_s_ml],df['t'],df[ 'f'], df['f_err']))
    tE_[1].append([t0_ml, u0_ml, tE_ml,f_s_ml])
    
#print tE_[1]
mm = np.asarray( tE_[0])
tE__ = tE_[1][mm.argmax()]

t0_ml, u0_ml, tE_ml,f_s_ml = tE__[0],tE__[1],tE__[2],tE__[3]

f_ris = df['f']-fun(t0_ml, u0_ml, tE_ml,f_s_ml)
f_ris_true = df['f_true']-fun(t0_ml, u0_ml, tE_ml,f_s_ml)

u0_true = float(f_s_true)/(max(df['f']-f_ris)-1+float(f_s_true))
t0_true = df['t'][(df['f']-f_ris).argmax()]

#print 't0_true = ' + str(t0_true)
        

        
N = np.arange(1,len(df['t'])+1)
chi_2 = chisq([t0_ml, u0_ml, tE_ml,f_s_ml],df['t'],df[ 'f'], df['f_err'])
t_chi = pd.DataFrame({'t' : df['t'],'chi2' : np.abs((chi_2)-N)/max(np.abs((chi_2)-N))})
        
half_chi_ = [t_chi['chi2'][(f_ris).argmax()],t_chi['chi2'][(f_ris).argmin()]]
t_chi_ =  [t_chi['t'][(f_ris).argmax()],t_chi['t'][(f_ris).argmin()]] 
f_ris__ = [(f_ris).max(),(f_ris).min()]
#half_chi = half_chi_[np.asarray([np.abs(l - 0.5) for l in half_chi_]).argmin()]
#t_half_chi = t_chi_[np.asarray([np.abs(l - 0.5) for l in half_chi_]).argmin()]
#f_ris_max = f_ris__[np.asarray([np.abs(l - 0.5) for l in half_chi_]).argmin()]
t_half_50 = lin_fit(slope(x_y_finder(t_chi['t'],t_chi['chi2'],0.5)),0.5)
jj = [0.05,0.2,0.3,0.4] #np.linspace(0.15,1-half_chi,4)
#print jj

#ttt =t_half_chi #(lin_fit(slope(x_y_finder(t_chi,half_chi)),half_chi))    
#duration = [min([np.abs((lin_fit(slope(x_y_finder(t_chi['t'],t_chi['chi2'],half_chi+j)),half_chi+j))-ttt ), np.abs(ttt-(lin_fit(slope(x_y_finder(t_chi['t'],t_chi['chi2'],half_chi-j)),half_chi-j)))]) for j in jj]
tt = t_half_50

duration = [min([np.abs((lin_fit(slope(x_y_finder(t_chi['t'],t_chi['chi2'],0.5+j)),0.5+j))-tt ), np.abs(tt-(lin_fit(slope(x_y_finder(t_chi['t'],t_chi['chi2'],0.5-j)),0.5-j)))]) for j in jj]
dur_min = min(duration)
duration.append(dur_min/10)
duration.append(dur_min*10)
duration.append(dur_min*100)


min_model__ = [[],[]]
min_model_ = [[],[]]
    

        

# making sure we are finding time of maximum absolute deviation:
                
if np.asarray([np.abs(l - 0.5) for l in half_chi_]).argmin() == 1:
    amp = 'neg'
else:
    amp = 'pos'
                        

for sigma in duration:
            
    for a in range(0,2):        
        amp_ = f_ris__[a]
        #print amp_
        t_mean_ = t_chi_[a]
        amp_ = amp_ * sigma * np.sqrt(2*pi)
        #print 'sigma = '+str(sigma) + ', amp = '+ str(amp_)
    
        nll = lambda *args: -lnlike2(*args)
        result = op.minimize(nll, [t_mean_,sigma,amp_,t0_ml, u0_ml, tE_ml,f_s_ml], args=(df['t'],df[ 'f'], df['f_err']),method = 'Nelder-Mead')
        mean_mll, sigma_mll,amp_mll,t0_mll, u0_mll, tE_mll,f_s_mll = result['x']
        #print result['x']
        min_model_[0].append(lnlike2([mean_mll, sigma_mll,amp_mll,t0_mll, u0_mll, tE_mll,f_s_mll],df['t'],df[ 'f'], df['f_err']))
        min_model_[1].append([mean_mll, sigma_mll,amp_mll,t0_mll, u0_mll, tE_mll,f_s_mll])

        

mmm_ = np.asarray( min_model_[0])



final_param = min_model_[1][mmm_.argmax()]
    
#print 'second fit =' + str(final_param)
    

    

x_c = np.sqrt((final_param[4]**2)+((final_param[3]-final_param[0])/final_param[5])**2 )

if amp == 'pos':
    s_exp = (x_c+np.sqrt((x_c**2)+4))/2   
elif amp == 'neg':
    s_exp = np.abs((x_c-np.sqrt((x_c**2)+4))/2)
#print 'theo q = ' +q +',theo s = ' +s
#print 'exp q = '+str((final_param[1]/final_param[5])**2)+ ', exp s = ' +str(s_exp)
    
# theoretical values of mass ratio
#q_[1].append(np.log10(float(q)))
# experimental values of mass ratio
#q_[0].append(np.log10((final_param[1]/final_param[5])**2))
# theoretical values of separation
#s_[1].append(np.log10(float(s)))
# experimental values of separation
#q_[0].append(np.log10(s_exp))

del_t = round((np.sqrt((float(s)-(1/float(s)))**2 - (float(u0_theo))**2)),2)
cell_text = [[round(final_param[3],2),round(float(t0_theo),2)],[round(final_param[4],2),round(float(u0_theo),2)],[round(final_param[5],3),round(float(tE_theo),3)],[round(final_param[6],3),round(float(f_s_true),3)],[round(final_param[0],2),str(round(float(t0_theo),2)-del_t) + ' or '+str(round(float(t0_theo),2)+del_t)],[round(final_param[1],3),round(np.sqrt(float(q))*float(tE_theo),3)],[round(final_param[2],2),0]]
rows = ('t0','u0','tE','f_s','tp','tEp','amp')   
columns = ('Fitted values','True values')

           
f_res_final = df['f']-fun2(final_param[0],final_param[1],final_param[2],final_param[3],final_param[4],final_param[5],final_param[6])                        
                                 
# definitions for the axes
left, width = 0.1, 0.8
bottom, height_1,height_2, height_3 = 0.05,0.1, 0.53 ,0.2
bottom_2 = bottom+ height_1 +0.02
bottom_3 = bottom_2 + height_2 +0.07

rect_1 = [left, bottom, width, height_1]
rect_2 = [left, bottom_2, width, height_2]
rect_3 = [left, bottom_3, width, height_3]                                            
                                                       
plt.close()
plt.figure(1, figsize=(15,12)) 
                                                          
residual = plt.axes(rect_1)
main = plt.axes(rect_2)
tbl = plt.axes(rect_3)                                                                            
                                                                                        
tbl.axis('tight')
tbl.axis('off')
the_table = tbl.table(cellText=cell_text,rowLabels=rows,colLabels=columns,loc='center',cellLoc='center')
prop = the_table.properties()
cells = prop['child_artists']
for cell in cells:
    cell.set_height(0.15)
main.set_title('Simulated WFIRST lc plus the functional modeling for '+'"'+str(name)+'"')
main.plot(df['t'],df['f'],'b.',df['t'],fun2(final_param[0],final_param[1],final_param[2],final_param[3],final_param[4],final_param[5],final_param[6]),'r')#plt.text(485,25, str(info), style='italic',bbox={'facecolor':'none','alpha':0.5, 'pad':10})
#main.axis([float(t0_theo)-(float(tE_theo)*1.5),float(t0_theo)+(float(tE_theo)*1.5),min(df['f'])-0.5,max(df['f'])+0.5])
main.axes.get_xaxis().set_ticks([])

residual.plot(df['t'],f_res_final,'r-')
#residual.axis([float(t0_theo)-(float(tE_theo)*1.5),float(t0_theo)+(float(tE_theo)*1.5),min(f_ris),max(f_ris)])
fig = plt.gcf()
#fig.savefig(home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/plots_single/'+str(name)+'_lc.png')

#plt.close()
plt.figure(2)
fig = plt.gcf()
fig.set_size_inches(12.0,6.0)
plt.title('Plot of normalized (cumulative chi-squared - N) VS time for '+'"'+str(name)+'"')
plt.plot(t_chi['t'],t_chi['chi2'],'k.')
#fig.savefig(home+'/Library/Mobile Documents/com~apple~CloudDocs/Microlensing/OSU trip/Matt/plots_single/'+str(name)+'_chi.png')
plt.show()


end = time.time()
print (end-start)


 