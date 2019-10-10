import numpy as np
from netCDF4 import Dataset
import os
import matplotlib.pyplot as plt


def readfile(monthstring,n):
    methane_12c= monthstring+'{}.nc'.format((str(n).zfill(2)))
    fh = Dataset(methane_12c, mode='r')  
    CH4 = fh.variables['IJ-AVG-S__CH4'][:]
    emissions = fh.variables['CH4-EMIS__CH4-TOT'][:]
    fh.close()
    MR =CH4[0,lat,lon]
    emis = np.sum(emissions[0,:,:])
    return MR, emis

def jacob_element(MR,emis):
    outval = np.true_divide(np.array(MR),emis)
    return outval

#------------------------------------------------

# Define times

#------------------------------------------------

nwindow       = 4 # number of months for the MAP window

monthstr = ['jan04','feb04','mar04','apr04','may04','jun04','jul04','aug04','sep04','oct04','nov04','dec04',
            'jan05','feb05','mar05','apr05','may05','jun05','jul05','aug05','sep05','oct05','nov05','dec05',
            'jan06','feb06','mar06','apr06','may06','jun06','jul06','aug06','sep06','oct06','nov06','dec06',
            'jan07','feb07','mar07','apr07','may07','jun07','jul07','aug07','sep07','oct07','nov07','dec07',
            'jan08','feb08','mar08','apr08','may08','jun08','jul08','aug08','sep08','oct08','nov08','dec08',
            'jan09','feb09','mar09','apr09','may09','jun09','jul09','aug09','sep09','oct09','nov09','dec09',
            'jan10','feb10','mar10','apr10','may10','jun10','jul10','aug10','sep10','oct10','nov10','dec10',
            'jan11','feb11','mar11','apr11','may11','jun11','jul11','aug11','sep11','oct11','nov11','dec11',
            'jan12','feb12','mar12','apr12','may12','jun12','jul12','aug12','sep12','oct12','nov12','dec12',
            'jan13','feb13','mar13','apr13','may13','jun13','jul13','aug13','sep13','oct13','nov13','dec13',
            'jan14','feb14','mar14','apr14','may14','jun14','jul14','aug14','sep14','oct14','nov14','dec14',
            'jan15','feb15','mar15','apr15','may15','jun15','jul15','aug15','sep15','oct15','nov15','dec15',
            'jan16','feb16','mar16','apr16','may16','jun16','jul16','aug16','sep16','oct16','nov16','dec16']


nmonths = len(monthstr)

#------------------------------------------------

# Read model mole fraction data

#------------------------------------------------

#GEOS-Chem grid boxes that correspond with each NOAA site included 
site = ['alt','asc','brw','cgo',
        'kum','mhd','mlo',
        'nwr','smo','tap','wlg']


lati=[43,21,40,12,27,36,27,33,19,32,32]
loni=[23,33, 5,65, 5,34, 5,15, 2,61,56]



nsites = len(lati) 

region=['australia','europe','eurasianboreal','eurasiantemperate','northafrica',

        'northamericanboreal','northamericantemperate','southamericantemperate',

        'oceans','southamericantropical','southernafrica','tropicalasia']

nregions = len(region)


MR12 = []
MR13 = []
emis = []

# read in MR and emissions data using readfile
# order: SITE x REGION x MONTH X RUN
# size: nsites x nregion x nmonths x nwindow

for sitecount in np.arange(nsites):
    lat = lati[sitecount]
    lon = loni[sitecount]
    
    for rgcount in np.arange(nregions):
        reg  =region[rgcount]
        os.chdir('/geos/d21/s1420671/jacobian/12CH4/14yr/jacobian/regions/'+reg)

        for monthcount in np.arange(nmonths):
            monthstring = monthstr[monthcount]

            for windowcount in np.arange(nwindow):
                windowuse = windowcount + 1

                MRtmp, emistmp = readfile(monthstring,windowuse)

                MR12 =np.append(MR12,MRtmp)
                #emis=np.append(emis,emistmp)


for sitecount in np.arange(nsites):
    lat = lati[sitecount]
    lon = loni[sitecount]
    
    for rgcount in np.arange(nregions):
        reg  =region[rgcount]
        os.chdir('/geos/d21/s1420671/jacobian/13CH4/14yr/jacobian/regions/'+reg)

        for monthcount in np.arange(nmonths):
            monthstring = monthstr[monthcount]

            for windowcount in np.arange(nwindow):
                windowuse = windowcount + 1

                MRtmp, emistmp = readfile(monthstring,windowuse)

                MR13=np.append(MR13,MRtmp)
                emis=np.append(emis,emistmp)


MR13=MR13*1.035
MR = (((MR13/MR12)/0.0112372)-1)*1000

mr=MR

##set up jacobian values by dividing MR by emissions
#split by each sensitivity test
#MR at specific site and emis is total for region 
#size = REGION  x SITE x MONTH
#size = 12 x 23 x 12

values = np.zeros([nwindow,(len(emis)//nwindow)]) 

i=0

for n in np.arange((len(emis)//nwindow)):
    MRUse  = MR[i:i+nwindow]
    outval = jacob_element(MRUse,emis[i])
    values[:,n] = outval
    i=i+nwindow                


#put jacobian together
#read out first 12 (region1, site 1)
#Then skip down 15 rows
#and move onto next site

jacobian = np.zeros([(nmonths+nwindow-1)*nsites,nmonths*nregions])

vl=0
col=0
mo=0

for sitecount in np.arange(nsites):
    col = 0
    
    for rgn in np.arange(nregions):
        mo = sitecount*(nmonths+nwindow-1)
        
        for month in np.arange(nmonths):
            jacobian[mo:mo+nwindow,col] = values[:,vl]
            
            vl  += 1
            col  += 1
            mo += 1

#------------------------------------------------

# Read prior fluxes

#------------------------------------------------

#every 4th contains a value
#only need to read out first 144 (12 regions, 12 months)
#as for each site within a region, the emissions are the same

xa = []

emiscount=0

for loop in np.arange(nmonths*nregions):
    xa= np.append(xa,emis[emiscount])
    emiscount=emiscount+nwindow        
 
##Read out a 'normal' GEOS-Chem run for mixing ratios
#with no region doubled 
#used for the caculation 
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

mix13 = []

os.chdir('/geos/d21/s1420671/13CH4/forward/tries/forward')

for lcount in np.arange(nsites):
    lat = lati[lcount]
    lon = loni[lcount]
    
    for n in np.arange(0,len(date),1):
        methane_12c= 'out.{}01.nc'.format((str(date[n])))
        fh = Dataset(methane_12c, mode='r')  
        CH4 = fh.variables['IJ-AVG-S__CH4'][:]
        fh.close()
        MixR =CH4[0,lat,lon]
        mix13 = np.append(mix13, MixR)
    
    
mix12 = []

os.chdir('/geos/d21/s1420671/12CH4/prior')

for lcount in np.arange(nsites):
    lat = lati[lcount]
    lon = loni[lcount]
    
    for n in np.arange(0,len(date),1):
        methane_12c= 'out.{}01.nc'.format((str(date[n])))
        fh = Dataset(methane_12c, mode='r')  
        CH4 = fh.variables['IJ-AVG-S__CH4'][:]
        fh.close()
        MixR =CH4[0,lat,lon]
        mix12 = np.append(mix12, MixR)
    


mix13=mix13*1.035
mix = (((mix13/mix12)/0.0112372)-1)*1000

  
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

sa_err = np.square(xa*0.5)
se_err = np.square(1)


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

xx = mix + multi_dot([inv(multi_dot([tranjac,invse,jacobian])+invsa),tranjac,invse,(y-mix)])

s = inv(multi_dot([jacobian.T,invse,jacobian])+invsa)


xx = abs(xx)
xxsum = np.array2string(np.sum(xx))
xasum = np.array2string(np.sum(xa))


#set up errors from error covariance matrix for post (d)
#and prior error

err=[]

for errcount in np.arange(nmonths*nregions):
    err = np.append(err,np.sqrt(s[errcount,errcount]))

err_scale = np.true_divide(err,xa)
   

pri_err =[]

for pri_count in np.arange(nmonths*nregions):
    pri_err = np.append(pri_err,np.sqrt(Sa[pri_count,pri_count]))

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
    plt.annotate("Avg post error = "+ avgerr, xy=(0.05,0.9), xycoords='axes fraction',fontsize=15)
    plt.title(rg,fontsize=20)
    plt.legend(('posterior', 'prior'),prop={'size':15})
    
    scale[pl:pl+nmonths] = np.around(np.divide(xx[pl:pl+nmonths],xa[pl:pl+nmonths]),decimals=5)
    
    if plot == 4:
        plt.ylabel('CH4 emissions (kg)',fontsize=40)
    if plot == 8:
        plt.xlabel('Year',fontsize=20)
    if plot == 9:
        plt.xlabel('Year',fontsize=20)
    if plot == 10:
        plt.xlabel('Year',fontsize=20)
    if plot == 11:
        plt.xlabel('Year',fontsize=20)

    pl=pl+nmonths
 
    
#os.chdir("/home/s1420671/python/Timeseries_Plots/long_inversion/CH4/inv_goodsites")
#plt.savefig("CH4_inversion_goodsite.jpg")

