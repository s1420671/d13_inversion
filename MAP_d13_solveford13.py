#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 09:09:23 2019

@author: s1420671
"""
import numpy as np
from netCDF4 import Dataset
import os
import matplotlib.pyplot as plt


dte = ['200401','200402','200403','200404','200405','200406','200407','200408',
        '200409','200410','200411','200412','200501','200502','200503','200504',
        '200505','200506','200507','200508','200509','200510','200511','200512',
        '200601','200602','200603','200604','200605','200606','200607','200608',
        '200609','200610','200611','200612','200701','200702','200703','200704',
        '200705','200706','200707','200708','200709','200710','200711','200712',
        '200801','200802','200803','200804','200805','200806','200807','200808',
        '200809','200810','200811','200812','200901','200902','200903','200904',
        '200905','200906','200907','200908','200909','200910','200911','200912',
        '201001','201002','201003','201004','201005','201006','201007','201008',
        '201009','201010','201011','201012','201101','201102','201103','201104',
        '201105','201106','201107','201108','201109','201110','201111','201112',
        '201201','201202','201203','201204','201205','201206','201207','201208',
        '201209','201210','201211','201212','201301','201302','201303','201304',
        '201305','201306','201307','201308','201309','201310','201311','201312',
        '201401','201402','201403','201404','201405','201406','201407','201408',
        '201409','201410','201411','201412','201501','201502','201503','201504',
        '201505','201506','201507','201508','201509','201510','201611','201512',
        '201601','201602','201603','201604','201605','201606','201607','201608',
        '201609','201610','201611','201612']

date = ['200401','200402','200403','200404','200405','200406','200407','200408',
        '200409','200410','200411','200412','200501','200502','200503','200504',
        '200505','200506','200507','200508','200509','200510','200511','200512',
        '200601','200602','200603','200604','200605','200606','200607','200608',
        '200609','200610','200611','200612','200701','200702','200703','200704',
        '200705','200706','200707','200708','200709','200710','200711','200712',
        '200801','200802','200803','200804','200805','200806','200807','200808',
        '200809','200810','200811','200812','200901','200902','200903','200904',
        '200905','200906','200907','200908','200909','200910','200911','200912',
        '201001','201002','201003','201004','201005','201006','201007','201008',
        '201009','201010','201011','201012','201101','201102','201103','201104',
        '201105','201106','201107','201108','201109','201110','201111','201112',
        '201201','201202','201203','201204','201205','201206','201207','201208',
        '201209','201210','201211','201212','201301','201302','201303','201304',
        '201305','201306','201307','201308','201309','201310','201311','201312',
        '201401','201402','201403','201404','201405','201406','201407','201408',
        '201409','201410','201411','201412','201501','201502','201503','201504',
        '201505','201506','201507','201508','201509','201510','201611','201512',
        '201601','201602','201603','201604','201605','201606','201607','201608',
        '201609','201610','201611','201612','201701','201702','201703']


site = ['alt','asc','brw','cgo',
        'kum','mhd','mlo',
        'nwr','smo','tap','wlg']
nstes=len(site)
LTI=[43,21,40,12,27,36,27,33,19,32,32]
LNI=[23,33, 5,65, 5,34, 5,15, 2,61,56]

sites = ['alt','asc','brw','cgo',
        'kum','mhd','mlo','und',
        'nwr','smo','tap','wlg']
nsites=len(sites)

lati=[15,34,39,28,26,40,30,15,29,22,17,25]
loni=[62,42,51,53,41, 5,20,24, 8,24,39,40]

region=['aus','euro','eubor','eutemp','nafrica', 'nabor','natemp',
        'satemp','oceans','satemp','safrica','tropicalasia']

nregions = len(region)

mix13 = []

os.chdir('/geos/d21/s1420671/13CH4/forward/tries/forward')

for n in np.arange(0,len(dte),1):
    os.chdir('/geos/d21/s1420671/13CH4/forward/tries/forward')
    #print(date[n])
    methane_12c= 'out.{}01.nc'.format((str(dte[n])))
    fh = Dataset(methane_12c, mode='r')  
    CH4 = fh.variables['IJ-AVG-S__CH4'][:]
    fh.close()
    MixR =CH4[0,:,:]
    
    for r in np.arange(nregions):
        reg=region[r]
        os.chdir('/geos/u73/s1420671/code/maps/transcom_regions')
        msk = '{}.nc'.format((str(reg)))
        fh = Dataset(msk, mode='r')  
        mask = fh.variables['map'][:]
        
        avd = np.mean(MixR * mask)
        mix13 = np.append(mix13, avd)
    
    
mix12 = []



for n in np.arange(0,len(dte),1):
    os.chdir('/geos/d21/s1420671/12CH4/prior')
    methane_12c= 'out.{}01.nc'.format((str(dte[n])))
    fh = Dataset(methane_12c, mode='r')  
    CH4 = fh.variables['IJ-AVG-S__CH4'][:]
    fh.close()
    MixR =CH4[0,:,:]
    
    for r in np.arange(nregions):
        reg=region[r]
        os.chdir('/geos/u73/s1420671/code/maps/transcom_regions')
        msk = '{}.nc'.format((str(reg)))
        fh = Dataset(msk, mode='r') 
        mask = fh.variables['map'][:]
        
        avd = np.mean(MixR * mask)
        mix12 = np.append(mix12, avd)


mix13=mix13*1.035
xa = (((mix13/mix12)/0.0112372)-1)*1000

  




mix_13 = []

os.chdir('/geos/d21/s1420671/13CH4/forward/tries/forward')

for lcount in np.arange(nstes):
    lat = lti[lcount]
    lon = lni[lcount]
    
    for n in np.arange(0,len(date),1):
        methane_12c= 'out.{}01.nc'.format((str(date[n])))
        fh = Dataset(methane_12c, mode='r')  
        CH4 = fh.variables['IJ-AVG-S__CH4'][:]
        fh.close()
        MixR =CH4[0,lat,lon]
        mix_13 = np.append(mix_13, MixR)
    
    
mix_12 = []

os.chdir('/geos/d21/s1420671/12CH4/prior')

for lcount in np.arange(nstes):
    lat = lti[lcount]
    lon = lni[lcount]
    
    for n in np.arange(0,len(date),1):
        methane_12c= 'out.{}01.nc'.format((str(date[n])))
        fh = Dataset(methane_12c, mode='r')  
        CH4 = fh.variables['IJ-AVG-S__CH4'][:]
        fh.close()
        MixR =CH4[0,lat,lon]
        mix_12 = np.append(mix_12, MixR)
    


mix_13=mix_13*1.035
mix = (((mix_13/mix_12)/0.0112372)-1)*1000


mx=[]
y=[]
      
for siten in np.arange(len(site)):
    #print(len(mx))
    si = site[siten]
    os.chdir('/home/s1420671/Downloads/13CH4_NOAA_observations')
    f = open("ch4c13_"+si+"_surface-flask_1_sil_month.txt", "r")
    for line in f:
        if '2004' in line:
           mx = np.append(mx,line)
        if '2005' in line:
           mx = np.append(mx,line)
        if '2006' in line:
           mx = np.append(mx,line)
        if '2007 ' in line:
           mx = np.append(mx,line)
        if '2008' in line:
           mx = np.append(mx,line)
        if '2009' in line:
           mx = np.append(mx,line)
        if '2010' in line:
           mx = np.append(mx,line)
        if '2011' in line:
           mx = np.append(mx,line)
        if '2012' in line:
           mx = np.append(mx,line)
        if '2013' in line:
           mx = np.append(mx,line)
        if '2014' in line:
           mx = np.append(mx,line)
        if '2015' in line:
            mx=np.append(mx,line)
        if '2016' in line:
           mx=np.append(mx,line)
        if '2017  1' in line:
            mx=np.append(mx,line)
        if '2017  2' in line:
            mx=np.append(mx,line)
        if '2017  3' in line:
            mx=np.append(mx,line)
                   
    f.close


for line in mx:
    spl = line.split()
    y = np.append(y, float(spl[3]))


#here, you can change the error in Sa and Se 
#and construct the errors matrices   
#Sa has 50% error on the prior
#Se has 10ppb error
sa_err = np.zeros(1872)

sa_err[0:156] = np.square(xa[0:156]*0.2)
sa_err[156:312] = np.square(xa[156:312]*0.5)
sa_err[312:468] = np.square(xa[312:468]*0.3)
sa_err[468:624] = np.square(xa[468:624]*0.7)
sa_err[624:780] = np.square(xa[624:780]*0.8)
sa_err[780:936] = np.square(xa[780:936]*0.1)
sa_err[936:1092] = np.square(xa[936:1092]*0.3)
sa_err[1092:1248] = np.square(xa[1092:1248]*0.8)
sa_err[1248:1404] = np.square(xa[1248:1404]*0.15)
sa_err[1404:1560] = np.square(xa[1404:1560]*0.8)
sa_err[1560:1716] = np.square(xa[1560:1716]*0.8)
sa_err[1716:1872] = np.square(xa[1716:1872]*0.8)

se_err = np.square(0.25)


Sa = np.zeros([nmonths*nregions,nmonths*nregions])
for ii in np.arange(nmonths*nregions): Sa[ii,ii] = sa_err[ii]

Se = np.zeros([len(y),len(y)])
for ii in np.arange(len(y)): Se[ii,ii] = se_err


#-----------------------------
# MAP
#-----------------------------

from numpy.linalg import inv
from numpy.linalg import multi_dot

invse = inv(Se)
tranjac = np.transpose(jacobian)
invsa = inv(Sa)

xx = xa + multi_dot([inv(multi_dot([tranjac,invse,jacobian])+invsa),tranjac,invse,(y-mix)])
#xx=-xx

s = inv(multi_dot([jacobian.T,invse,jacobian])+invsa)


#xx = np.negative(xx)
xxsum = np.array2string(np.sum(xx))
xasum = np.array2string(np.sum(xa))


#set up errors from error covariance matrix for post (d)
#and prior error

err=[]

for errcount in np.arange(nmonths*nregions):
    err = np.append(err,np.sqrt(s[errcount,errcount]))

err_scale = np.true_divide(err,xa)
err[156:312]=err[156:312]*0.48 
err[468:780]=err[468:780]*0.28 
err[936:1248]=err[936:1248]*0.23
err[1404:1872]=err[1404:1872]*0.18 

pri_err =[]

for pri_count in np.arange(nmonths*nregions):
    pri_err = np.append(pri_err,np.sqrt(Sa[pri_count,pri_count]))
    
pri_err[156:312]=pri_err[156:312]*0.5
pri_err[468:780]=pri_err[468:780]*0.3 
pri_err[936:1248]=pri_err[936:1248]*0.25 
pri_err[1404:1872]=pri_err[1404:1872]*0.2 

#plot subplots of posterior estimate
#on same axes as prior
#plot errorbars as shaded area
#write out relative error on posterior
#store scaling factor to go from prior to posterior emissions
#to use in a forward run for the NOAA comparison

import pandas as pd
x = pd.date_range(start='1/1/2004', end='1/1/2017', periods=nmonths)

scale=[]
pl=0

fig=plt.figure()
fig.set_figheight(20)
fig.set_figwidth(35)

region=['australia','europe','eurasianboreal','eurasiantemperate','northafrica',

        'northamericanboreal','northamericantemperate','southamericantemperate',

        'oceans','southamericantropical','southernafrica','tropicalasia']

for plot in np.arange(nregions):

    rg =region[plot]
    plt.subplot(3,4,plot+1)  
    plt.rc('font', size=20)
    plt.plot(x,xx[pl:pl+nmonths], x,xa[pl:pl+nmonths])
    plt.fill_between(x,xx[pl:pl+nmonths]-err[pl:pl+nmonths],
                     xx[pl:pl+nmonths]+err[pl:pl+nmonths], alpha=0.4, facecolor='#089FFF')
    plt.fill_between(x,xa[pl:pl+nmonths]-pri_err[pl:pl+nmonths],
                     xa[pl:pl+nmonths]+pri_err[pl:pl+nmonths], alpha=0.4,facecolor='#FF9848')
    avgerr = np.around(np.mean(err_scale[pl:pl+nmonths]),decimals=5)
    avgerr = np.array2string(avgerr)
    #plt.annotate("Avg post error = "+ avgerr, xy=(0.05,0.9), xycoords='axes fraction',fontsize=15)
    plt.title(rg,fontsize=20)
    plt.legend(('posterior', 'prior'),prop={'size':15})
    
    scale[pl:pl+nmonths] = np.around(np.divide(xx[pl:pl+nmonths],xa[pl:pl+nmonths]),decimals=5)
    
    if plot == 4:
        plt.ylabel('d13C (ppt)',fontsize=40)
    if plot == 8:
        plt.xlabel('Year',fontsize=20)
    if plot == 9:
        plt.xlabel('Year',fontsize=20)
    if plot == 10:
        plt.xlabel('Year',fontsize=20)
    if plot == 11:
        plt.xlabel('Year',fontsize=20)
    plt.ylim(-70,-30)

    pl=pl+nmonths
 