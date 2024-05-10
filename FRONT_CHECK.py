#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 10:06:06 2022

@author: jono
"""
# Will compare fronts to vertical velocity as well as current speed

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage.util import random_noise
from skimage import feature
import datetime
from netCDF4 import Dataset
import os
os.environ['PROJ_LIB'] = '/home/jono/anaconda3/share/proj/'
from mpl_toolkits.basemap import Basemap
import matplotlib
import cmocean
import cmocean.cm as cmo
import pickle
from geopy import distance

#%%

def extract_CROCO(file_name, file_base = '/media/data/CHPC_SBUS_3km/', var = 'temp', level = 59, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.):
    import numpy as np
    from netCDF4 import Dataset
    nc_file = file_base + file_name
    root = Dataset(nc_file)
    lat = root.variables['lat_rho'][:]
    lon = root.variables['lon_rho'][:]
    time_CROCO = root.variables['time'][:]
    
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx], idx
    
    lat_min_index = find_nearest(lat[:,0],lat_min)[1]
    lat_max_index = find_nearest(lat[:,0],lat_max)[1]
    lon_min_index = find_nearest(lon[0,:],lon_min)[1]
    lon_max_index = find_nearest(lon[0,:],lon_max)[1]
    
    lon_CROCO = lon[0,lon_min_index:lon_max_index]
    lat_CROCO = lat[lat_min_index:lat_max_index,0]
    
    # Quick check
    if len(root.variables[var].shape) == 4 and str(level).isnumeric() == True:
        var_CROCO = root.variables[var][:,level,lat_min_index:lat_max_index, lon_min_index:lon_max_index]
    elif len(root.variables[var].shape) == 4 and level == False:
        var_CROCO = root.variables[var][:,:,lat_min_index:lat_max_index, lon_min_index:lon_max_index]
    elif len(root.variables[var].shape) == 3 and level == False:
        var_CROCO = root.variables[var][:,lat_min_index:lat_max_index, lon_min_index:lon_max_index]
    else:
        pass
    
    var_CROCO = np.squeeze(var_CROCO)
    h_CROCO = root.variables['h'][lat_min_index:lat_max_index, lon_min_index:lon_max_index]
    mask_CROCO = root.variables['mask_rho'][lat_min_index:lat_max_index, lon_min_index:lon_max_index]
    pn = root.variables['pn'][lat_min_index:lat_max_index, lon_min_index:lon_max_index] 
    pm = root.variables['pm'][lat_min_index:lat_max_index, lon_min_index:lon_max_index] 
    f = root.variables['f'][lat_min_index:lat_max_index, lon_min_index:lon_max_index] 
    mask_CROCO[mask_CROCO == 0 ] = float("NAN")
    var_CROCO = np.multiply(var_CROCO,mask_CROCO)
    h_CROCO = np.multiply(h_CROCO,mask_CROCO)
    f = np.multiply(f,mask_CROCO)
    pn = np.multiply(pn,mask_CROCO)
    pm = np.multiply(pm,mask_CROCO)
    root.close()
    return var_CROCO, lon_CROCO, lat_CROCO, h_CROCO, pn, pm, f, time_CROCO

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def contour_loc(contour,lon,lat):
    Xloc = np.empty((0))
    Yloc = np.empty((0))
    for k in range(0,len(contour)):
        lat_index = find_nearest(lat,contour[k,1])[1]
        lon_index = find_nearest(lon,contour[k,0])[1]
        Yloc = np.append(Yloc,lat_index)
        Xloc = np.append(Xloc,lon_index)
    return(Xloc,Yloc)


#%% Read in the TIME data

file_base =  '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT/'
Y_min = 2004
Y_max = 2018

# Get lon, lat and topography from the first file

lon_CROCO, lat_CROCO, h_CROCO, pn, pm, f =  extract_CROCO('croco_avg_Y'+str(Y_min)+'M'+str(1)+'.nc',file_base=file_base,var='temp', 
                                 level = 59, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[1:7]

GEOS = np.empty((0,len(lat_CROCO),len(lon_CROCO)))
TIME = np.empty((0))
for i in range(Y_min, Y_max+1):
    for j in range(1,12+1):
        file_name = 'croco_avg_Y'+str(i)+'M'+str(j)+'.nc'
        if file_name != 'croco_avg_Y2018M12.nc':
            print("Processing "  + file_name)
            nc_file = file_base + file_name
            root = Dataset(nc_file)
            time_CROCO = root.variables['time'][:]
            ug = extract_CROCO(file_name,file_base=file_base,var='u', 
                                level = 59, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[0]
            vg = extract_CROCO(file_name,file_base=file_base,var='v', 
                                level = 59, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[0]
            gg = np.sqrt(ug**2+vg**2)
            GEOS = np.concatenate((GEOS, gg), axis=0)
            TIME = np.append(TIME,time_CROCO)
        else:
            print("No Data")
            
# Any repeated times, remove them

uniqueValues, indicesList = np.unique(TIME, return_index=True)
TIME = TIME[indicesList]    
GEOS = GEOS[indicesList,:,:]   
# Convert the time array to a date

TIME = TIME.data   # Masked array appernetly
d0 = datetime.datetime(1990,1,1)
DATE = np.empty((0))
for k in range(0,len(TIME)):
       dt = datetime.timedelta(seconds = TIME[k])
       date  = d0 + dt
       DATE = np.append(DATE,date)

# Break down date array into months so we can get the seasons

Y = np.empty((0))
M = np.empty((0))
for dd in range(0,len(TIME)):
    year = DATE[dd].year
    month = DATE[dd].month
    Y = np.append(Y,year)
    M = np.append(M,month)
    
Summer = np.array([1,2,12])
Winter = np.array([6,7,8])

# Summer

ind_summer = np.empty((0)) 
for ii in range(0,len(Summer)):    
    ind_date = np.where(M == Summer[ii])[0]
    ind_summer = np.append(ind_summer,ind_date)

ind_summer = np.sort(ind_summer)

# Winter
ind_winter = np.empty((0)) 
for ii in range(0,len(Winter)):    
    ind_date = np.where(M == Winter[ii])[0]
    ind_winter = np.append(ind_winter,ind_date)  
    
ind_winter = np.sort(ind_winter)  

# Average our data

GEOS = np.transpose(GEOS,(1,2,0))  # Fix the array orientation such that time is third dimension

#######################
# SAVE the GEOS data
#######################

f = open('GEOS.pckl', 'wb')
pickle.dump(GEOS, f)
f.close()    

#%% 

####################
# LOAD in the EDGES
####################

f = open('EDGE.pckl', 'rb')
EDGE = pickle.load(f)
f.close()

########################
# LOAD vertical velocity
########################

f = open('W.pckl', 'rb')
W = pickle.load(f)
f.close()

########################
# LOAD BVF
########################

f = open('BVF.pckl', 'rb')
BVF = pickle.load(f)
f.close()

#%%

############################
# OBJECTIVE
############################

# Design two subplots of snapshots of geostrophic velocity as well as vertical 
# velocity which will be overliad with the edges to see their correspondence to the jet locations and upwelling
# regions. Choose a year and then do random days with the apprpriate edge

# Pick a year

YY = 2016
ind_YY = np.where(Y == YY)[0]

DATE_YY = DATE[ind_YY]

DATE_start = np.empty((0)) 
for ii in range(0,len(DATE_YY)):    
    date = DATE_YY[ii].replace(day=1)
    DATE_start = np.append(DATE_start,date) 
    
# Only get the unique dates

uniqueValues, indicesList = np.unique(DATE_start, return_index=True)
DATE_start = DATE_start[indicesList]   

DATE_ind = np.empty((0))
for ii in range(0,len(DATE_start)):
    ind_date = np.where(DATE_YY == DATE_start[ii])[0]
    DATE_ind = np.append(DATE_ind,ind_date) 

ind_YY = ind_YY[DATE_ind.astype(int)]

# Current velocity

# define subplot grid
fig, axs = plt.subplots(nrows=4, ncols=3, figsize=(15, 15))
fig.suptitle("Current speed vs edges", fontsize=18, y=0.95)
tickers = np.array(ind_YY)
# loop through tickers and axes
for ticker, ax in zip(tickers, axs.ravel()):
    map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
                llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=ax, resolution='l')
    map.fillcontinents(color='grey')
    map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
    map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
    x,y = map(lon_CROCO,lat_CROCO)
    X, Y = np.meshgrid(x, y)
    im = ax.pcolor(X,Y,GEOS[:,:,ticker], cmap=cmo.speed, vmin = 0, vmax = 1)
    ax.contour(X,Y,EDGE[:,:,ticker], linewidths=0.2, linestyles='solid', colors='black')
    ax.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
    map.fillcontinents(color='grey')

    # chart formatting
    ax.set_title(DATE[ticker].strftime("%Y-%m-%d"))

plt.colorbar(im, shrink=0.5,ax=axs.ravel().tolist()).set_label(label="Speed [m/s]",size=13)
plt.show()

fig.savefig('Currentvsedge_APP1.png', format='png', dpi=600)
plt.close()

# Vertical velcoity

fig, axs = plt.subplots(nrows=4, ncols=3, figsize=(15, 15))
fig.suptitle("Vertical velocity vs edges", fontsize=18, y=0.95)
tickers = np.array(ind_YY)
# loop through tickers and axes
for ticker, ax in zip(tickers, axs.ravel()):
    map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
                llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=ax, resolution='l')
    map.fillcontinents(color='grey')
    map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
    map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
    x,y = map(lon_CROCO,lat_CROCO)
    X, Y = np.meshgrid(x, y)
    im = ax.pcolor(X,Y,W[:,:,ticker], cmap=cmo.balance, vmin = -5e-5, vmax = 5e-5)
    ax.contour(X,Y,EDGE[:,:,ticker], linewidths=0.2, linestyles='solid', colors='black')
    ax.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
    map.fillcontinents(color='grey')

    # chart formatting
    ax.set_title(DATE[ticker].strftime("%Y-%m-%d"))

plt.colorbar(im, shrink=0.5,ax=axs.ravel().tolist()).set_label(label="Velocity [m/s]",size=13)
plt.show()

fig.savefig('WvsedgesAPP2.png', format='png', dpi=600)
plt.close()

# BVF

fig, axs = plt.subplots(nrows=4, ncols=3, figsize=(15, 15))
fig.suptitle("BvF @ 30 m vs edges", fontsize=18, y=0.95)
tickers = np.array(ind_YY)
# loop through tickers and axes
for ticker, ax in zip(tickers, axs.ravel()):
    map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
                llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=ax, resolution='l')
    map.fillcontinents(color='grey')
    map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
    map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
    x,y = map(lon_CROCO,lat_CROCO)
    X, Y = np.meshgrid(x, y)
    im = ax.pcolor(X,Y,BVF[:,:,ticker], cmap=cmo.dense)#, vmin = -5e-5, vmax = 5e-5)
    ax.contour(X,Y,EDGE[:,:,ticker], linewidths=0.2, linestyles='solid', colors='black')
    ax.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
    map.fillcontinents(color='grey')

    # chart formatting
    ax.set_title(DATE[ticker].strftime("%Y-%m-%d"))

plt.colorbar(im, shrink=0.5,ax=axs.ravel().tolist()).set_label(label="Velocity [m/s]",size=13)
plt.show()

fig.savefig('WvsedgesAPP2.png', format='png', dpi=600)
plt.close()