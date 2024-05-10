#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 11:28:24 2022

@author: jono
"""
# The idea is to compute Ekman transport along the coast.

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
import math
import pandas as pd
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
    f = np.multiply(f,mask_CROCO)
    var_CROCO = np.multiply(var_CROCO,mask_CROCO)
    h_CROCO = np.multiply(h_CROCO,mask_CROCO)
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

WIND = np.empty((0,len(lat_CROCO),len(lon_CROCO)))
TIME = np.empty((0))
angle = math.radians(22)
for i in range(Y_min, Y_max+1):
    for j in range(1,12+1):
        file_name = 'croco_avg_Y'+str(i)+'M'+str(j)+'.nc'
        if file_name != 'croco_avg_Y2018M12.nc':
            print("Processing "  + file_name)
            nc_file = file_base + file_name
            root = Dataset(nc_file)
            time_CROCO = root.variables['time'][:]            
            v_wind = extract_CROCO(file_name,file_base=file_base,var='svstr', 
                                      level = False, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[0]
            u_wind = extract_CROCO(file_name,file_base=file_base,var='sustr', 
                                      level = False, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[0]
            # Rotate the vectors in space
            u_r = u_wind * math.cos(angle) - v_wind * math.sin(angle)
            v_r = u_wind * math.sin(angle) + v_wind * math.cos(angle)
            TIME = np.append(TIME,time_CROCO)
            WIND = np.concatenate((WIND, v_r), axis=0)
        else:
            print("No Data")
            
# Any repeated times, remove them

uniqueValues, indicesList = np.unique(TIME, return_index=True)
TIME = TIME[indicesList]    
WIND = WIND[indicesList,:,:]

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

WIND = np.transpose(WIND,(1,2,0))  # Fix the array orientation such that time is third dimension

#%%

##################
# CHECK
##################

# We need to roate the wind vectors to be parrallel with the coast. To do this, we first need to plot
# a quiver plot showing the average summer and winter state and then rotate the vectors

fig, axes = plt.subplots(1,2,figsize=(12,9))

# Normal

axes[0].set_title("Normal")
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0], resolution='l')
map.fillcontinents(color='grey')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
axes[0].quiver(X[::20, ::20], Y[::20, ::20], np.nanmean(u_wind,0)[::20, ::20], np.nanmean(v_wind,0)[::20, ::20],
               pivot='mid', scale=1.5)
axes[0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')

# Rotated

axes[1].set_title("Rotated")
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[1], resolution='l')
map.fillcontinents(color='grey')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[1].quiver(X[::20, ::20], Y[::20, ::20], np.nanmean(u_r,0)[::20, ::20], np.nanmean(v_r,0)[::20, ::20],
               pivot='mid', scale=1.5)
axes[1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')

plt.show()
#fig.savefig('WIND_ROT_APP3.png', format='png', dpi=600)
plt.close()

#%% 

###################################
# CREATE ALL THE MASKS
###################################

#######################################
#      Create the various masks
#######################################

shelf_depth = 500
h_mask = h_CROCO.copy()
h_mask[h_mask > shelf_depth] = float("NAN")
mask_canny = ~np.isnan(h_mask)

#                                   STEP 1:
# Design mask to encapsulate only the SBUS. To calculate the displacement of the fronts
# will use the location of the coast. Find the lon and lat pair for the contour defining
# the coast and then use it as a reference when calculating the displacment of the Shelf-break
# and upwelling front

# Find the location of the coast

ind_coast_X = np.empty((0))
ind_coast_Y = np.empty((0))
for i in range(0,np.size(h_CROCO,0)):
    loc = np.argwhere(np.isnan(h_CROCO[i,:]))
    if len(loc) > 1:
        ind_coast_X = np.append(ind_coast_X,min(loc))
        ind_coast_Y = np.append(ind_coast_Y,i)
    else:
            ind_coast_X = np.append(ind_coast_X,float("NAN"))
            ind_coast_Y = np.append(ind_coast_Y,i)
        
#lon_coast = lon_CROCO[ind_coast_X.astype(int)]
#lat_coast = lat_CROCO[ind_coast_Y.astype(int)]

# Plot the coast contour to check
ind_nan = np.argwhere(np.isnan(ind_coast_X))

fig = plt.figure()
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
for i in range(53, 340):
    plt.plot(x[ind_coast_X[i].astype(int)],
         y[ind_coast_Y[i].astype(int)], 'ok', markersize=5)
plt.show()

# Create the necessary mask

CC_lat = -32.82    # Latitude of Cape Columbine
SBUS_mask = h_mask.copy()

ind_CC = find_nearest(lat_CROCO,CC_lat)[1]

SBUS_mask[0:ind_CC,:] = float("NAN")

# Create specific mask for shelf-break: 500 - 200 m and upwelling front: <200 m

# Shelf-break front

sb_mask = SBUS_mask.copy()
sb_mask[sb_mask < 200]  = float("NAN")
sb_mask[sb_mask > 500] = float("NAN")
sb_mask = ~np.isnan(sb_mask)

# Upwelling front

up_mask = SBUS_mask.copy()
up_mask[up_mask > 200]  = float("NAN")
up_mask = ~np.isnan(up_mask)

# Correct for slight mismatch in topography

sb_mask[221:233,26:35] = True
up_mask[221:233,26:35] = False

#%% 

#######################################
# GET THE ALONGSHORE WIND
#######################################

WIND_COAST = np.empty((len(lat_CROCO),0))
for j in range(0,len(TIME)):
    wind_coast = []
    for i in range(0,len(ind_coast_X)):
        if np.isnan(ind_coast_X[i]):
            out = float("NAN")
            wind_coast.append(out)  # Vertical velocity
            #print("NAN")
        else:
            ind_x = list(range(ind_coast_X[i].astype(int)-4,ind_coast_X[i].astype(int)-1))
            out = np.nanmean(WIND[i,ind_x,j])
            wind_coast.append(out)
    wind_coast = np.array(wind_coast)        
    WIND_COAST = np.concatenate((WIND_COAST, wind_coast[...,np.newaxis]), axis=1)
    #WIND_COAST.append(wind_coast)                 
            
###########################################
# Subdivide the SBUS into three regions
###########################################

# The SBUS will be divided into three: St Helena Bay [-32.82 -31.5]
# Mid-shelf [-31.5 -30.3] and Namaqua [-30.3 -28]

###########################################
# Calculate the offshore EKMAN transport
############################################

# Me_x = tau_y/f
# Need to extract the wind component parallel to the coast

def EKMAN_TRANSPORT(wind,f):
    rho_0 = 1025
    transport = np.multiply(wind,1/f[...,np.newaxis])
    transport = transport/rho_0
    return transport

# St Helena Mask

CO_lat_min = -32.82   # Latitude of Cut min
CO_lat_max = -31.5 
ind_CO_min = find_nearest(lat_CROCO,CO_lat_min)[1]
ind_CO_max = find_nearest(lat_CROCO,CO_lat_max)[1]
ff = f[ind_CO_min:ind_CO_max,0]
WIND_SH = WIND_COAST[ind_CO_min:ind_CO_max,:] 
ME_SH = EKMAN_TRANSPORT(WIND_SH,ff)
ME_SH = np.nanmean(ME_SH,0)

# Mid-shelf

CO_lat_min = -31.5   # Latitude of Cut min
CO_lat_max = -30.3   # Latitude of Cut max 
ind_CO_min = find_nearest(lat_CROCO,CO_lat_min)[1]
ind_CO_max = find_nearest(lat_CROCO,CO_lat_max)[1]
ff = f[ind_CO_min:ind_CO_max,0]
WIND_MS = WIND_COAST[ind_CO_min:ind_CO_max,:]
ME_MS = EKMAN_TRANSPORT(WIND_MS,ff)
ME_MS = np.nanmean(ME_MS,0)
# Namaqua

CO_lat_min = -30.3   # Latitude of Cut min
CO_lat_max = -28   # Latitude of Cut max 
ind_CO_min = find_nearest(lat_CROCO,CO_lat_min)[1]
ind_CO_max = find_nearest(lat_CROCO,CO_lat_max)[1]
ff = f[ind_CO_min:ind_CO_max,0]
WIND_NM = WIND_COAST[ind_CO_min:ind_CO_max,:]
ME_NM = EKMAN_TRANSPORT(WIND_NM,ff) 
ME_NM = np.nanmean(ME_NM,0)

#%%

#####################
# LOAD Data frame
#####################

f = open('ST_HELENA.pckl', 'rb')
ST_HELENA = pickle.load(f)
f.close()

df_shelf = ST_HELENA['Shelf_break']
df_up = ST_HELENA['Upwelling']

# Add data to the dataframes

df_shelf['EKMAN'] = ME_SH
df_up['EKMAN'] = ME_SH

#%%
#######################################################
# Plot combining Ekman transport and vertical velocity
#######################################################

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection 

# Compute the monthly means for each variable

data_columns = ['DIS', 'LEN', 'NUM','GRAD','W', 'EKMAN']

mon_shelf = df_shelf[data_columns].resample('M').mean()
day3_shelf = df_shelf[data_columns].rolling(3).mean()
mon_up = df_up[data_columns].resample('M').mean()
day3_up = df_up[data_columns].rolling(3).mean()

##########################
# MAIN FIGURE
##########################

fig, axes = plt.subplots(2,1,sharex=True, figsize=(18,7))

# Shelf-break
dates = day3_shelf.index
y = np.array(day3_shelf.EKMAN)
c = np.array(day3_shelf.W)
             
axes[0].set_title("(a) Shelf-break region",loc='left', fontsize=16)
#convert dates to numbers first
inxval = mdates.date2num(dates)
years = mdates.YearLocator()
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="RdBu_r", linewidth=1)

# Add the monthly values
axes[0].scatter(x=mon_shelf.index, y=mon_shelf['EKMAN'],s=90,c=mon_shelf['W'],cmap='RdBu_r',
                vmin=-2e-5, vmax=2e-5)
# Add line to show mean
axes[0].axhline(y=np.nanmean(df_shelf.EKMAN), xmin=0, xmax=1, linewidth=2, color='r')
# set color to date values
lc.set_array(c)
line = axes[0].add_collection(lc)
line.set_clim(vmin=-2e-5, vmax=2e-5)
axes[0].tick_params(axis="y", labelsize=14) 
axes[0].grid()
axes[0].set_ylabel("Transport [kg/m^2/s]",fontsize=16)
#axes[0].set_ylim([80, 180])

# Upwelling
dates = day3_up.index
y = np.array(day3_up.EKMAN)
c = np.array(day3_up.W)

axes[1].set_title("(b) Upwelling region",loc='left', fontsize=16)
#convert dates to numbers first
inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="RdBu_r", linewidth=1)
# Add the monthly values
axes[1].scatter(x=mon_up.index, y=mon_up['EKMAN'],s=90,c=mon_up['W'],cmap='RdBu_r',
                vmin =-2e-5, vmax=2e-5)
# Add line to show mean
axes[1].axhline(y=np.nanmean(df_up.EKMAN), xmin=0, xmax=1, linewidth=2, color='r')
# set color to date values
lc.set_array(c)
line = axes[1].add_collection(lc)
axes[1].xaxis.set_major_locator(years)
line.set_clim(vmin=-2e-5, vmax=2e-5)
axes[1].grid()
axes[1].set_xlabel("DATE",fontsize=16)
axes[1].set_ylabel("Transport [kg/m^2/s]",fontsize=16)
#axes[1].set_ylim([30, 100])
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a
plt.colorbar(line, ax=axes.ravel().tolist()).set_label(label="Vertical velocities [m/s]",size=13)
axes[1].xaxis_date()
plt.show()

fig.savefig('NAM_EKMAN&W', format='png', dpi=600)

#%%

###############################################
# CLIMATOLOGY FIGURES
###############################################

# The 3-day time-series are usueful and show the huge range of interannual
# variability, but we are interested in the seasonality of the fronts

# Restrict distance domain: shelf, upwelling
# MS: 80-130; 0-50
# NM: 110-160; 40-90
# SH: 80-120; 30-70

# Construct a subplot with the EKMAN transport vs Displacement and EKMAN Transport vs vertical
# for each region

df_clim_shelf = df_shelf.groupby([df_shelf.index.month]).mean() # CLIM
df_clim_up = df_up.groupby([df_up.index.month]).mean() # CLIM

fig, axes = plt.subplots(2,2,sharex=True, figsize=(18,7))

labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul',
          'Aug','Sep','Oct','Nov','Dec']

# Dsiplacement plot

# Shelf-break
dates = df_clim_shelf.index
y = np.array(df_clim_shelf.DIS)
c = np.array(df_clim_shelf.EKMAN)
             
axes[0,0].set_title("(a) Shelf-break region",loc='left', fontsize=16)
#convert dates to numbers first
inxval = list(range(1, 13))
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="viridis", linewidth=1)

# Add the monthly values
axes[0,0].scatter(x=df_clim_shelf.index, y=df_clim_shelf['DIS'],s=90,c=df_clim_shelf['EKMAN'],cmap='viridis',
                vmin=-500, vmax=100)
# Add line to show mean
axes[0,0].axhline(y=np.nanmean(df_clim_shelf.DIS), xmin=0, xmax=1, linewidth=2, color='r')
# set color to date values
lc.set_array(c)
line = axes[0,0].add_collection(lc)
line.set_clim(vmin=-500, vmax=100)
axes[0,0].tick_params(axis="y", labelsize=14) 
axes[0,0].grid()
axes[0,0].set_ylabel("Distance [km]",fontsize=16)
axes[0,0].set_ylim([80, 120])

# Upwelling
dates = df_clim_up.index
y = np.array(df_clim_up.DIS)
c = np.array(df_clim_up.EKMAN)

axes[0,1].set_title("(b) Upwelling region",loc='left', fontsize=16)
#convert dates to numbers first
#inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="viridis", linewidth=1)
# Add the monthly values
axes[0,1].scatter(x=df_clim_up.index, y=df_clim_up['DIS'],s=90,c=df_clim_up['EKMAN'],cmap='viridis',
                vmin =-500, vmax=100)
# Add line to show mean
axes[0,1].axhline(y=np.nanmean(df_clim_up.DIS), xmin=0, xmax=1, linewidth=2, color='r')
# set color to date values
lc.set_array(c)
line = axes[0,1].add_collection(lc)
line.set_clim(vmin=-500, vmax=100)
axes[0,1].tick_params(axis="y", labelsize=14)
axes[0,1].grid()
axes[0,1].set_ylabel("Distance [km]",fontsize=16)
axes[0,1].set_ylim([30, 70])
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a
plt.colorbar(line, ax=axes[0,:], label="Trasnport [km/m^2/s]")

# Vertical velcoity

# Shelf-break
dates = df_clim_shelf.index
y = np.array(df_clim_shelf.EKMAN)
c = np.array(df_clim_shelf.W)

axes[1,0].set_title("(c) Shelf-break region",loc='left', fontsize=16)
#convert dates to numbers first
#inxval = mdates.date2num(dates)
#years = mdates.YearLocator()
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="RdBu_r", linewidth=1)
# Add the monthly values
axes[1,0].scatter(x=df_clim_shelf.index, y=df_clim_shelf['EKMAN'],s=90,c=df_clim_shelf['W'],cmap='RdBu_r',
                vmin=-5e-6, vmax=5e-6)
# Add a mean line
axes[1,0].axhline(y=np.nanmean(df_clim_shelf.EKMAN), xmin=0, xmax=1, linewidth=2, color='b')
# set color to date values
lc.set_array(c)
line = axes[1,0].add_collection(lc)
line.set_clim(vmin=-5e-6, vmax=5e-6) 
axes[1,0].grid()
axes[1,0].tick_params(axis="x", labelsize=14)
axes[1,0].tick_params(axis="y", labelsize=14)
axes[1,0].xaxis.set_ticks(inxval)
axes[1,0].xaxis.set_ticklabels(labels,rotation=45,ha = 'right')
axes[1,0].set_ylabel("Transport",fontsize=16)

# Upwelling
dates = df_clim_up.index
y = np.array(df_clim_up.EKMAN)
c = np.array(df_clim_up.W)

axes[1,1].set_title("(d) Upwelling region",loc='left', fontsize=16)
#convert dates to numbers first
#inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="RdBu_r", linewidth=1)
# Add the monthly values
axes[1,1].scatter(x=df_clim_up.index, y=df_clim_up['EKMAN'],s=90,c=df_clim_up['W'],cmap='RdBu_r',
                vmin=-5e-6, vmax=5e-6)
# Add a mean line
axes[1,1].axhline(y=np.nanmean(df_clim_up.EKMAN), xmin=0, xmax=1, linewidth=2, color='b')
# set color to date values
lc.set_array(c)
line = axes[1,1].add_collection(lc)
line.set_clim(vmin=-5e-6, vmax=5e-6)
axes[1,1].grid()
axes[1,1].set_ylabel("Transport",fontsize=16)
axes[1,1].xaxis.set_ticks(inxval)
axes[1,1].xaxis.set_ticklabels(labels,rotation=45,ha = 'right')
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a

plt.colorbar(line, ax=axes[1, :], label='Vertical velocities [m/s]' )
plt.show()

fig.savefig('ST_HELENA_EKMAN', format='png', dpi=600)
plt.close()

#%%

######################################
# EKMAN TRANSPORT 
######################################

# One figure showing the EKMAN TRANSPORT for the 
# Upwelling region of the different regions on one plot

# Create a data-frame

from matplotlib.transforms import Affine2D

data = {'Date': DATE,
        'ST_HELENA': ME_SH,
        'MID_SHELF': ME_MS,
        'NAMAQUA': ME_NM}
EKMAN = pd.DataFrame(data) 
EKMAN = EKMAN.set_index('Date')

EKMAN_clim = EKMAN.groupby([EKMAN.index.month]).mean() # CLIM
EKMAN_clim_STD = EKMAN.groupby([EKMAN.index.month]).std() # CLIM

NAME = ["ST_HELENA","MID_SHELF","NAMAQUA"]
MARKER = ["o","v","s"]
COLOR = ["coral","deepskyblue","limegreen"]
labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul',
          'Aug','Sep','Oct','Nov','Dec']

#fig = plt.figure(figsize=(9,8))
fig, ax = plt.subplots(figsize=(9,8))
ax.set_facecolor('lightgrey')

trans1 = Affine2D().translate(-0.2, 0.0) + ax.transData
trans2 = Affine2D().translate(0.2, 0.0) + ax.transData
trans3 = Affine2D().translate(0.0, 0.0) + ax.transData

TRANS = [trans1,trans2,trans3]

for i in range(0,len(NAME)):
    inxval = list(range(1, 13))
#    plt.fill_between(inxval, 
#                     EKMAN_clim[NAME[i]] - EKMAN_clim_STD[NAME[i]], 
#                     EKMAN_clim[NAME[i]] + EKMAN_clim_STD[NAME[i]],
#                           alpha=0.25,
#                           edgecolor='black',
#                           facecolor=COLOR[i]) 
#    plt.plot(inxval,EKMAN_clim[NAME[i]],
#             marker=MARKER[i], 
#             markersize=14,
#             color= COLOR[i],
#             label=NAME[i])
    plt.errorbar(inxval, EKMAN_clim[NAME[i]], marker=MARKER[i], 
                 markersize=14,
                 color= COLOR[i],
                 label=NAME[i],
                 yerr = EKMAN_clim_STD[NAME[i]],
                 fmt="",
                 ecolor = COLOR[i],
                 elinewidth=3,
                 transform = TRANS[i])
    
plt.xticks(inxval, labels, rotation=45)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel("Transport [$m^{2}$ $s^{-1}$]",fontsize=18)
plt.legend(borderpad=1)
plt.grid()
plt.show()


fig.savefig('SBUS_EKMAN', format='png', dpi=600)
plt.close()



