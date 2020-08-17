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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import (mark_inset,inset_axes,InsetPosition)
from tqdm import tqdm

start = time.time()

home = os.path.expanduser("~")

name = 'rerun_chi_squared.csv' 
tempdata = 'data/'+str(name) #home+'/Desktop/trial runs/'+str(name)
dff = pd.read_csv(tempdata)

direc = os.listdir(".")

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
    
    if len(t_right) != 0:
        tE_right = np.abs(t_right[(np.abs(df['f'][df['t']>t_max]-(A_base))).argmin()]-t_max)
    else:
        tE_right = 1
    if len(t_left) != 0:
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

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
        
s_ = [[],[]]
q_ = [[],[]]
tE_list =[[],[]]
t0_ =[[],[]]
u0_ =[[],[]]
f_s_ =[[],[]]
tp_ =[]
tEp_ =[[],[]]
ampli_ =[]
name_lc = []
chi_2_1 = []
chi_2_2 = []
del_t_ = []
x_c_ = []
s1 , s2 = [], []

s_f = [[],[]]
q_f = [[],[]]
tE_list_f =[[],[]]
t0_f =[[],[]]
u0_f =[[],[]]
f_s_f =[[],[]]
tp_f =[]
tEp_f =[[],[]]
ampli_f =[]
name_lc_f = []
chi_2_1_f = []
chi_2_2_f = []
del_t_f = []
x_c_f = []
s1_f , s2_f = [], []
length, length_f = [], []
alpha_ = []
rho_ = []
m_l_ = []
D_l_ = []
D_s_ = []
mu_s_l_ = []
mu_s_b_ = []
mu_l_l_ = []
mu_l_b_ = []
alpha_f = []
rho_f = []
m_l_f = []
D_l_f = []
D_s_f = []
mu_s_l_f = []
mu_s_b_f = []
mu_l_l_f = []
mu_l_b_f = []


# definitions for the axes
left, width = 0.1, 0.8
bottom, height_1,height_2, height_3 = 0.05,0.1, 0.53 ,0.2
bottom_2 = bottom+ height_1 +0.02
bottom_3 = bottom_2 + height_2 +0.07

rect_1 = [left, bottom, width, height_1]
rect_2 = [left, bottom_2, width, height_2]
rect_3 = [left, bottom_3, width, height_3]

alllc1 = pd.read_csv(home+'/Desktop/alllc1/alllc1_info.CSV')
alllc2 = pd.read_csv(home+'/Desktop/alllc2/alllc2_info.CSV')


for i in tqdm( range(len(dff['name']))): #len(dff['name'])
    
    if dff['name'][i] in np.asarray(alllc1['name']):
	tempdata1 = home+'/Desktop/alllc1/'+dff['name'][i]
    if dff['name'][i] in np.asarray(alllc2['name']):
	tempdata1 = home+'/Desktop/alllc2/'+dff['name'][i]
    #print dff['name'][i]

    if dff['name'][i].endswith('.gz'):
        t,f,f_err,f_true,code = np.loadtxt(tempdata1,usecols=(0,1,2,3,5),unpack=True)
        df = pd.DataFrame({'t':t , 'f':f , 'f_err' : f_err , 'f_true': f_true, 'code':code})
        df = df[df['code']==4] # WFIRST filter W149
        df = df.reset_index(drop=True) #reset indeces

        #df['f'] = df['f_true']
        
        fname = gzip.open(tempdata1, 'rb')
        x_0 = fname.readlines()[0:7]
        f_s_true = x_0[0].split(' ')[4] # WFIRST filter W149
        q = x_0[5].split(' ')[5]
        s = x_0[5].split(' ')[6]
        tE_theo = x_0[6].split(' ')[4]
        t0_theo = x_0[6].split(' ')[3]
        u0_theo = x_0[6].split(' ')[1]
        alpha = (x_0[6].split(' ')[2])
        rho = (x_0[6].split(' ')[5])
        m_l = (x_0[4].split(' ')[3])
        D_l = (x_0[4].split(' ')[2])
        D_s = (x_0[3].split(' ')[2])
        mu_s_l = (x_0[3].split(' ')[4])
        mu_s_b = (x_0[3].split(' ')[5])
        mu_l_l = (x_0[4].split(' ')[4])
        mu_l_b = (x_0[4].split(' ')[5])
        
        t_len = len(df['t'])
        
    
        A_max = 1.0/(float(f_s_true)/(max(df['f'])-1+float(f_s_true)))
        u0_true = np.sqrt( ( (1+np.sqrt(1+16*(A_max**2)))/(2* A_max) ) - 2 )
        t0_true = df['t'][df['f'].argmax()]
        ind1, ind2 = fwhm(df['f'],df['f'].argmax(),min(df['f']))
        #tE_true = t[ind2]-t[ind1]
        tE_true = [tE_finder (df['t'],df[ 'f'], f_s_true),t[ind2]-t[ind1]]
        #print 'tE_true = '+ str( tE_true)

        tE_ = [[],[]]    
        for j in tE_true:
               
            nll = lambda *args: -lnlike(*args)
            result = op.minimize(nll, [t0_true, u0_true, j,f_s_true], args=(df['t'],df[ 'f'], df['f_err']),method = 'Nelder-Mead')
            t0_ml, u0_ml, tE_ml,f_s_ml = result['x']
            tE_[0].append(lnlike([t0_ml, u0_ml, tE_ml,f_s_ml],df['t'],df[ 'f'], df['f_err']))
            tE_[1].append([t0_ml, u0_ml, tE_ml,f_s_ml])
    
        #print tE_[1]
        mm = np.asarray( tE_[0])
        tE__ = tE_[1][mm.argmax()]

        t0_ml, u0_ml, tE_ml,f_s_ml = tE__[0],tE__[1],tE__[2],tE__[3]

        f_ris = df['f']-fun(t0_ml, u0_ml, tE_ml,f_s_ml)
        f_ris_true = df['f_true']-fun(t0_ml, u0_ml, tE_ml,f_s_ml)
        




        #u0_true = float(f_s_true)/(max(df['f']-f_ris)-1+float(f_s_true))
        A_max = 1.0/(float(f_s_true)/(max(df['f']-f_ris)-1+float(f_s_true)))
        u0_true = np.sqrt( ( (1+np.sqrt(1+16*(A_max**2)))/(2* A_max) ) - 2 )
        t0_true = df['t'][(df['f']-f_ris).argmax()]

        #print 't0_true = ' + str(t0_true)
        
        duration = [0.01, 0.1,1]  

     
        cut = 100     
        #f_ris__ = [(f_ris[(df['t']>(float(t0_true)-cut))&(df['t']<(float(t0_true)+cut))]).max(),(f_ris[(df['t']>(float(t0_true)-cut))&(df['t']<(float(t0_true)+cut))]).min()]
        #t_ris__ =  [df['t'][(f_ris[(df['t']>(float(t0_true)-cut))&(df['t']<(float(t0_true)+cut))]).argmax()], df['t'][(f_ris[(df['t']>(float(t0_true)-cut))&(df['t']<(float(t0_true)+cut))]).argmin()]]  
        

        f_ris = smooth(f_ris,10)
        f_ris_true = smooth(f_ris_true,10)

        if u0_ml>0.7:
            f_ris = smooth(f_ris,19)
            f_ris_true = smooth(f_ris_true,19)

        f_ris__ = [f_ris.max(),f_ris.min()]
        t_ris__ = [df['t'][f_ris.argmax()],df['t'][f_ris.argmin()]]

        
        min_model_ = [[],[],[]]
                
        # making sure we are finding time of maximum absolute deviation:
                  

        for sigma in duration:
            
            
            for a in range(0,2):
                
                amp_ = f_ris__[a]
                t_mean_ = t_ris__[a]
                amp_ = amp_ * sigma * np.sqrt(2*np.pi)
                #print 'sigma = '+str(sigma) + ', amp = '+ str(amp_)
            
                nll = lambda *args: -lnlike2(*args)
                result = op.minimize(nll, [t_mean_,sigma,amp_,t0_ml, u0_ml, tE_ml,f_s_ml], args=(df['t'],df[ 'f'], df['f_err']),method = 'Nelder-Mead')
                mean_mll, sigma_mll,amp_mll,t0_mll, u0_mll, tE_mll,f_s_mll = result['x']
                #print result['x']
            
                min_model_[0].append(lnlike2([mean_mll, sigma_mll,amp_mll,t0_mll, u0_mll, tE_mll,f_s_mll],df['t'],df[ 'f'], df['f_err']))
                min_model_[1].append([mean_mll, sigma_mll,amp_mll,t0_mll, u0_mll, tE_mll,f_s_mll])

        

        mmm_ = np.asarray( min_model_[0])


        print 2*(np.abs(mm.max())-np.abs(mmm_.max()))        
        
        final_param = min_model_[1][mmm_.argmax()]
        #print 'second fit =' + str(final_param)
    
        
    

        x_c = np.sqrt((final_param[4]**2)+((final_param[3]-final_param[0])/final_param[5])**2 )

        if final_param[2] >0:
            s_exp = (x_c+np.sqrt((x_c**2)+4))/2   
        elif final_param[2] <0:
            s_exp = np.abs(((x_c)-np.sqrt((x_c**2)+4))/2)

        del_t = round( (np.sqrt( np.abs((float(s)-(1/float(s)))**2  - (float(u0_theo) )**2) )) ,3)* float(tE_theo)
        
        #(np.abs(final_param[0] - final_param[3])/np.abs(del_t) <10) or ( np.abs(mmm_.max())<22000 )
        
        if  (np.abs(final_param[0] - final_param[3])<50) and  ( np.abs(mmm_.max() - mm.max()) >25) and ( np.abs(mmm_.max())<26000 ) :
        
        

            f_res_final = df['f']-fun2(final_param[0],final_param[1],final_param[2],final_param[3],final_param[4],final_param[5],final_param[6]) 

                   
        
            cell_text = [[round(final_param[3],2),round(float(t0_theo),2)],[round(final_param[4],2),round(float(u0_theo),2)],[round(final_param[5],3),round(float(tE_theo),3)],
            [round(final_param[6],3),round(float(f_s_true),3)],[round(final_param[0],2),str(round(float(t0_theo)-del_t,2)) + ' or '+str(round(float(t0_theo)+del_t,2))],
            [round(final_param[1],4),round(np.sqrt(float(q))*float(tE_theo),4)],[round(final_param[2],4),0]]
        
        
            rows = ('t0','u0','tE','f_s','tp','tEp','amp')   
            columns = ('Fitted values','True values')
                
                                               
                                                            
            plt.close()
            plt.figure(1, figsize=(15,12)) 
                                                            
            residual = plt.axes(rect_1)
            main = plt.axes(rect_2)
            tbl = plt.axes(rect_3)        
        
            tp_range = np.array([float(t0_theo)-del_t,float(t0_theo)+del_t])
                                                                    
         
            tp_theo = tp_range[np.abs(tp_range-final_param[0]).argmin()]

            range_p = [final_param[0]-((np.sqrt(float(q))*float(tE_theo))*3),     
                    final_param[0]+((np.sqrt(float(q))*float(tE_theo))*3)]

            if len(df['f'][(df['t']>range_p[0]) & (df['t']<range_p[1])]) == 0:
                range_p = [float(tp_theo)-((np.sqrt(float(q))*float(tE_theo))*3),     
                            float(tp_theo)+((np.sqrt(float(q))*float(tE_theo))*3)]
               
            f_range_p = df['f'][(df['t']>range_p[0]) & (df['t']<range_p[1])]
            t_range_p = df['t'][(df['t']>range_p[0]) & (df['t']<range_p[1])]
         
            khengool = float(tE_theo)*1.5
            if del_t> np.abs(final_param[0]-float(t0_theo)):
                while (del_t> khengool ):
                    khengool = khengool + 0.1 * khengool
            elif del_t< np.abs(final_param[0]-float(t0_theo)):
                while (np.abs(final_param[0]-float(t0_theo)) > khengool ):
                    khengool = khengool + 0.1 * khengool
        
            main_t_range =  [ float(t0_theo)-khengool,float(t0_theo)+khengool]                                                                                                                                                         
            if len(t_range_p)==0:
                print dff['name'][i], ' : light curve not complete'
            if not len(t_range_p)==0:
                
                # theoretical values of mass ratio
                q_[0].append(np.log10(float(q)))
                # experimental values of mass ratio
                q_[1].append(np.log10((final_param[1]/final_param[5])**2))
                # theoretical values of separation
                s_[0].append(np.log10(float(s)))
                # experimental values of separation
                s_[1].append(np.log10(s_exp))
                # name
                name_lc.append(dff['name'][i])
                # fitted tp
                tp_.append(final_param[0])
                #theoretical tEp
                tEp_[0].append(round(np.sqrt(float(q))*float(tE_theo),4))
                # fitted tEp
                tEp_[1].append(final_param[1])
                # fitted amplitude
                ampli_.append(final_param[2])  
                #theoretical t0
                t0_[0].append(float(t0_theo)) 
                # fitted t0
                t0_[1].append(final_param[3])
                # Theoretical tE
                tE_list[0].append(float(tE_theo))
                # fitted tE
                tE_list[1].append(float(final_param[5]))
                # Theoretical u0
                u0_[0].append(float(u0_theo))     
                # fitted u0
                u0_[1].append(final_param[4])
                # Theoretical f_s
                f_s_[0].append(float(f_s_true))
                # fitted f_s
                f_s_[1].append(final_param[6])
                # minimum chi squared 2
                chi_2_2.append(mmm_.max())
                # minimum chi squared 1
                chi_2_1.append(mm.max())
                # Save |tp-t0| true
                del_t_.append(del_t)
                # calculate x_c
                x_c_.append(x_c)
                #save degenarate values of s for one value of X_C
                s1.append((x_c+np.sqrt((x_c**2)+4))/2)
                s2.append ( np.abs(((x_c)-np.sqrt((x_c**2)+4))/2 ))
                #save number of data points
                length.append(t_len)
                alpha_.append(float(alpha))
                rho_.append(float(rho))
                m_l_ .append(float(m_l))
                D_l_.append(float(D_l))
                D_s_.append(float(D_s))
                mu_s_l_.append(float(mu_s_l))
                mu_s_b_.append(float(mu_s_b))
                mu_l_l_.append(float(mu_l_l))
                mu_l_b_.append(float(mu_l_b))
                
            
                                                                        
                                                                                                                                                                                                                                                                    
                tbl.axis('off')
                the_table = tbl.table(cellText=cell_text,rowLabels=rows,colLabels=columns,loc='center',cellLoc='center')
                prop = the_table.properties()
                cells = prop['child_artists']
                for cell in cells:
                    cell.set_height(0.15)
                main.set_title('Simulated WFIRST lc plus the functional modeling for '+'"'+str(dff['name'][i])+'"')
                main.plot(df['t'],df['f'],'b.',df['t'],fun2(final_param[0],final_param[1],final_param[2],final_param[3],final_param[4],final_param[5],final_param[6]),'r')
                main.axis([main_t_range[0],main_t_range[1],min(df['f']),max(df['f'])+(0.05*(max(df['f'])-min(df['f'])))])

                main.axes.get_xaxis().set_ticks([])

                if df['t'][df['f'].argmax()] >856 :
                    inset = plt.axes([left+0.01, bottom_3-0.07-0.2-0.05, 0.2, 0.2])
                else:
                    inset = plt.axes([0.9-0.2-0.01, bottom_3-0.07-0.2-0.05, 0.2, 0.2])

                inset.set_title('Planetary Event')
                inset.plot(df['t'],df['f'],'b.',df['t'],fun2(final_param[0],final_param[1],final_param[2],final_param[3],final_param[4],final_param[5],final_param[6]),'r')
                x1, x2, y1, y2 = min(t_range_p),max(t_range_p),min(f_range_p),max(f_range_p)
                inset.set_xlim(x1, x2)
                inset.set_ylim(y1-(0.05*(y2-y1)), y2+(0.05*(y2-y1)))


                inset.axes.get_xaxis().set_ticks([])
                inset.axes.get_yaxis().set_ticks([])
                mark_inset(main, inset, loc1=2, loc2=4, fc="none", ec="0.5")


                residual.plot(df['t'],f_res_final,'r-')
                residual.axis([main_t_range[0],main_t_range[1],min(f_res_final),max(f_res_final)+(0.05*(max(f_res_final)-min(f_res_final)))])
                fig = plt.gcf()
                fig.savefig(home+'/Desktop/trial_runs/new_tests/chi_2_test/'+str(dff['name'][i])+'.png')
                
        else:
            
            
            
            
            
            cell_text = [[round(final_param[3],2),round(float(t0_theo),2)],[round(final_param[4],2),round(float(u0_theo),2)],[round(final_param[5],3),round(float(tE_theo),3)],
            [round(final_param[6],3),round(float(f_s_true),3)],[round(final_param[0],2),str(round(float(t0_theo)-del_t,2)) + ' or '+str(round(float(t0_theo)+del_t,2))],
            [round(final_param[1],4),round(np.sqrt(float(q))*float(tE_theo),4)],[round(final_param[2],4),0]]
        
        
            rows = ('t0','u0','tE','f_s','tp','tEp','amp')   
            columns = ('Fitted values','True values')
                
            f_res_final = df['f']-fun2(final_param[0],final_param[1],final_param[2],final_param[3],final_param[4],final_param[5],final_param[6])                                    
                                                            
            plt.close()
            plt.figure(1, figsize=(15,12)) 
                                                            
            residual = plt.axes(rect_1)
            main = plt.axes(rect_2)
            tbl = plt.axes(rect_3)        
        
            tp_range = array([float(t0_theo)-del_t,float(t0_theo)+del_t])
                                                                    
         
            tp_theo = tp_range[np.abs(tp_range-final_param[0]).argmin()]

            range_p = [final_param[0]-((np.sqrt(float(q))*float(tE_theo))*3),     
                    final_param[0]+((np.sqrt(float(q))*float(tE_theo))*3)]

            if len(df['f'][(df['t']>range_p[0]) & (df['t']<range_p[1])]) == 0:
                range_p = [float(tp_theo)-((np.sqrt(float(q))*float(tE_theo))*3),     
                            float(tp_theo)+((np.sqrt(float(q))*float(tE_theo))*3)]
               
            f_range_p = df['f'][(df['t']>range_p[0]) & (df['t']<range_p[1])]
            t_range_p = df['t'][(df['t']>range_p[0]) & (df['t']<range_p[1])]
         
            khengool = float(tE_theo)*1.5
            if del_t> np.abs(final_param[0]-float(t0_theo)):
                while (del_t> khengool ):
                    khengool = khengool + 0.1 * khengool
            elif del_t< np.abs(final_param[0]-float(t0_theo)):
                while (np.abs(final_param[0]-float(t0_theo)) > khengool ):
                    khengool = khengool + 0.1 * khengool
        
            main_t_range =  [ float(t0_theo)-khengool,float(t0_theo)+khengool]                                                                                                                                                         
            if len(t_range_p)==0:
                print dff['name'], ' : light curve not complete'
            if not len(t_range_p)==0:
                                
                # theoretical values of mass ratio
                q_f[0].append(np.log10(float(q)))
                # experimental values of mass ratio
                q_f[1].append(np.log10((final_param[1]/final_param[5])**2))
                # theoretical values of separation
                s_f[0].append(np.log10(float(s)))
                # experimental values of separation
                s_f[1].append(np.log10(s_exp))
                # name
                name_lc_f.append(dff['name'][i])
                # fitted tp
                tp_f.append(final_param[0])
                #theoretical tEp
                tEp_f[0].append(round(np.sqrt(float(q))*float(tE_theo),4))
                # fitted tEp
                tEp_f[1].append(final_param[1])
                # fitted amplitude
                ampli_f.append(final_param[2])  
                #theoretical t0
                t0_f[0].append(float(t0_theo)) 
                # fitted t0
                t0_f[1].append(final_param[3])
                # Theoretical tE
                tE_list_f[0].append(float(tE_theo))
                # fitted tE
                tE_list_f[1].append(float(final_param[5]))
                # Theoretical u0
                u0_f[0].append(float(u0_theo))     
                # fitted u0
                u0_f[1].append(final_param[4])
                # Theoretical f_s
                f_s_f[0].append(float(f_s_true))
                # fitted f_s
                f_s_f[1].append(final_param[6])
                # minimum chi squared 2
                chi_2_2_f.append(mmm_.max())
                # minimum chi squared 1
                chi_2_1_f.append(mm.max())
                # Save |tp-t0| true
                del_t_f.append(del_t)
                # calculate x_c
                x_c_f.append(x_c)
                #save degenarate values of s for one value of X_C
                s1_f.append((x_c+np.sqrt((x_c**2)+4))/2)
                s2_f.append ( np.abs(((x_c)-np.sqrt((x_c**2)+4))/2 ))
                #save number of data points
                length_f.append(t_len)
                alpha_f.append(float(alpha))
                rho_f.append(float(rho))
                m_l_f .append(float(m_l))
                D_l_f.append(float(D_l))
                D_s_f.append(float(D_s))
                mu_s_l_f.append(float(mu_s_l))
                mu_s_b_f.append(float(mu_s_b))
                mu_l_l_f.append(float(mu_l_l))
                mu_l_b_f.append(float(mu_l_b))
                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                tbl.axis('tight')
                tbl.axis('off')
                the_table = tbl.table(cellText=cell_text,rowLabels=rows,colLabels=columns,loc='center',cellLoc='center')
                prop = the_table.properties()
                cells = prop['child_artists']
                for cell in cells:
                    cell.set_height(0.15)
                main.set_title('Simulated WFIRST lc plus the functional modeling for '+'"'+str(dff['name'][i])+'"')
                main.plot(df['t'],df['f'],'b.',df['t'],fun2(final_param[0],final_param[1],final_param[2],final_param[3],final_param[4],final_param[5],final_param[6]),'r')
                main.axis([main_t_range[0],main_t_range[1],min(df['f']),max(df['f'])+(0.05*(max(df['f'])-min(df['f'])))])
                main.axes.get_xaxis().set_ticks([])

                if df['t'][df['f'].argmax()] >856 :
                    inset = plt.axes([left+0.01, bottom_3-0.07-0.2-0.05, 0.2, 0.2])
                else:
                    inset = plt.axes([0.9-0.2-0.01, bottom_3-0.07-0.2-0.05, 0.2, 0.2])

                inset.set_title('Planetary Event')
                inset.plot(df['t'],df['f'],'b.',df['t'],fun2(final_param[0],final_param[1],final_param[2],final_param[3],final_param[4],final_param[5],final_param[6]),'r')
                x1, x2, y1, y2 = min(t_range_p),max(t_range_p),min(f_range_p),max(f_range_p)
                inset.set_xlim(x1, x2)
                inset.set_ylim(y1-(0.05*(y2-y1)), y2+(0.05*(y2-y1)))
    
                inset.axes.get_xaxis().set_ticks([])
                inset.axes.get_yaxis().set_ticks([])
                mark_inset(main, inset, loc1=2, loc2=4, fc="none", ec="0.5")


                residual.plot(df['t'],f_res_final,'r-')
                residual.axis([main_t_range[0],main_t_range[1],min(f_res_final),max(f_res_final)+(0.05*(max(f_res_final)-min(f_res_final)))])
                fig = plt.gcf()
                fig.savefig(home+'/Desktop/trial_runs/new_tests/chi_2_test/'+'failure?_'+str(dff['name'][i])+'.png')
                
            
            
result = pd.DataFrame({'name': name_lc ,'q_true_log': q_[0],'q_fitted_log': q_[1],'s_true_log': s_[0],'s_fitted_log': s_[1],
                    'q_true' : 10**np.asarray((q_[0])) , 'q_fitted' : 10**np.asarray((q_[1])), 's_true' : 10**np.asarray(s_[0]), 's_fitted' : 10**np.asarray(s_[1]),
                    'tp_fitted' : tp_ , 'tEp_true' : tEp_[0], 'tEp_fitted' : tEp_[1], 'ampl_fitted' : ampli_, 'del_t' : del_t_, 's_>_1' : s1,
                    'tE_true' : tE_list[0], 'tE_fitted' : tE_list[1], 't0_true' : t0_[0], 't0_fitted' : t0_[1], 'x_c' : x_c_, 's_<_1' : s2,
                    'u0_true' : u0_[0], 'u0_fitted' : u0_[1], 'f_s_true' : f_s_[0], 'f_s_fitted': f_s_[1], 'chi_2_2' : chi_2_2, 'chi_2_1' : chi_2_1,
                    'length': length, 'alpha':alpha_, 'rho': rho_, 'm_l': m_l_, 'D_l': D_l_, 'D_s':D_s_ , 'mu_s_l':mu_s_l_, 'mu_s_b':mu_s_b_,
                    'mu_l_l': mu_l_l_, 'mu_l_b': mu_l_b_})
                    
failure_list = pd.DataFrame({'name': name_lc_f ,'q_true_log': q_f[0],'q_fitted_log': q_f[1],'s_true_log': s_f[0],'s_fitted_log': s_f[1],
                    'q_true' : 10**np.asarray((q_f[0])) , 'q_fitted' : 10**np.asarray((q_f[1])), 's_true' : 10**np.asarray(s_f[0]), 's_fitted' : 10**np.asarray(s_f[1]),
                    'tp_fitted' : tp_f , 'tEp_true' : tEp_f[0], 'tEp_fitted' : tEp_f[1], 'ampl_fitted' : ampli_f, 'del_t' : del_t_f, 's_>_1' : s1_f,
                    'tE_true' : tE_list_f[0], 'tE_fitted' : tE_list_f[1], 't0_true' : t0_f[0], 't0_fitted' : t0_f[1], 'x_c' : x_c_f, 's_<_1' : s2_f,
                    'u0_true' : u0_f[0], 'u0_fitted' : u0_f[1], 'f_s_true' : f_s_f[0], 'f_s_fitted': f_s_f[1], 'chi_2_2' : chi_2_2_f, 'chi_2_1' : chi_2_1_f,
                    'length':length_f, 'alpha':alpha_f, 'rho': rho_f, 'm_l': m_l_f, 'D_l': D_l_f, 'D_s':D_s_f , 'mu_s_l':mu_s_l_f, 'mu_s_b':mu_s_b_f,
                    'mu_l_l': mu_l_l_f, 'mu_l_b': mu_l_b_f})                    

                    
result.to_csv(home+'/Desktop/trial_runs/new_tests/chi_2_test/alllc_w_c.CSV')
failure_list.to_csv(home+'/Desktop/trial_runs/new_tests/chi_2_test/alllc_w_c_failed.CSV')


plt.close()
plt.figure(3)
plt.title('Mass Ratio_Method1')
plt.plot(result['q_true_log'],result['q_true_log'],'g-',result['q_true_log'],result['q_fitted_log'],'b.')
plt.xlabel('True Mass Ratio - log',size=15)
plt.ylabel('Fitted Mass Ratio - log',size=15)
fig = plt.gcf()
fig.set_size_inches(12.0,6.0)
fig.savefig(home+'/Desktop/trial_runs/new_tests/chi_2_test/Mass Ratio.png')
plt.close()
plt.figure(4)
plt.title('Projected Separation_Method1')
plt.plot(result['s_true_log'],result['s_true_log'],'g-',result['s_true_log'],result['s_fitted_log'],'b.')
plt.xlabel('True Separation - log',size=15)
plt.ylabel('Fitted Separation - log',size=15)
fig = plt.gcf()
fig.set_size_inches(12.0,6.0)
fig.savefig(home+'/Desktop/trial_runs/new_tests/chi_2_test/Projected Separation.png')



end = time.time()
print (end-start)/(60*60)




 