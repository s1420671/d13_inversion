#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:51:26 2019

@author: s1420671
"""
import numpy as np
from netCDF4 import Dataset
import os

site = ['alt','asc','brw','cgo',
        'kum','mhd','mlo',
        'nwr','smo','tap','wlg']

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
                
    f.close

#for l in np.arange(len(mx)-1):
#    if mx[l]==mx[l+1]:
#        np.remove(mx[l])

for line in mx:
    spl = line.split()
    y = np.append(y, float(spl[3]))

#y =list(dict.fromkeys(y))

MR   = []

lati=[43,21,40,12,27,36,27,33,19,32,32]
loni=[23,33, 5,65, 5,34, 5,15, 2,61,56]


nsites = len(lati)

#read in MR and emissions data, read through the GC output files for emissions
#roll through regions and take out MR and emission
#order: R1S1,R1S2..R2S1,R2S2...

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
        '201609','201610','201611','201612']

nmonths = len(date)


c=0               
d13_pri=[]
pri12=[]
pri13=[]
post13=[]

os.chdir('/geos/d21/s1420671/jacobian/12CH4/14yr/forward/goodsites')
for c in np.arange(nsites):
        lat = lati[c]
        lon = loni[c]
        for n in np.arange(0,len(date),1):
            methane= 'out.{}01.nc'.format(str(date[n]))
            fh = Dataset(methane, mode='r') 
            C12 = fh.variables['IJ-AVG-S__CH4'][:]
            fh.close()
            pri12 = np.append(pri12, C12[0,lat,lon])
            
os.chdir('/geos/d21/s1420671/13CH4/forward/fulllength')
for c in np.arange(nsites):
        lat = lati[c]
        lon = loni[c]
        for n in np.arange(0,nmonths,1):
            methane= 'out.{}01.nc'.format((str(date[n])))
            fh = Dataset(methane, mode='r') 
            C13 = fh.variables['IJ-AVG-S__CH4'][:]
            fh.close()
            pri13 = np.append(pri13, C13[0,lat,lon])
        
d13_pri = (((pri13/(pri12))/(0.01118))-1)*1000               

#mr=MR
#c=0
#d13_post=[]
#os.chdir('/geos/d21/s1420671/jacobian/13CH4/forward/invfirsttry')
#for c in np.arange(nsites):
#        lat = lati[c]
#        lon = loni[c]
#        for n in np.arange(1,13,1):
#            methane= 'out.2014{}01.000000.nc'.format((str(n).zfill(2)))
#            fh = Dataset(methane, mode='r') 
#            C13 = fh.variables['IJ-AVG-S__CH4'][:]
#            fh.close()
#            post13 = np.append(post13, C13[0,lat,lon])
#            
#            
#d13_post = (((post13/(pri12))/(0.01118))-1)*1000  

#respri =   (y*0.97) - (MR_pri*1)
#respost = (y*0.97) - (MR_post*1)
respri = y - d13_pri
#respost = y - d13_post


#print(np.mean(respri),np.mean(respost))
avrespri = np.array2string(np.mean(respri))
#avrespost = np.array2string(np.mean(respost))


import matplotlib.pyplot as plt

#colors = ['darkorange','dodgerblue']

#plt.hist([respri,respost],20,color=colors,rwidth=1)
plt.hist(respri,30,color='darkorange',rwidth=1,alpha=0.8)
plt.xlim(-10,60)
#plt.hist(respost,20,color='dodgerblue',rwidth=1,alpha=0.8)
plt.xlabel('Residual value')
plt.ylabel('Count')
plt.legend(('Prior','Posterior'))
plt.title('Mean prior residual = '+avrespri)#; Mean posterior residual = '+avrespost)
#
#os.chdir("/home/s1420671/python/Timeseries_Plots/MAP_setup/bestNOAA")
#plt.savefig("residuals_bestNOAAabs.jpg")
