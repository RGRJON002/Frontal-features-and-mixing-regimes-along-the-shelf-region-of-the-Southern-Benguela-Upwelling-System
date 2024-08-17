#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 20:16:37 2022

@author: jono
"""
# This script aims to create a nice figure combining the time-series of the alongshore wind
# stratification and alongshore mixing regimes for a year


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

#%% Read in the TIME data and v component

file_base =  '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT/'
Y_min = 2004
Y_max = 2018

# Get lon, lat and topography from the first file

lon_CROCO, lat_CROCO, h_CROCO, pn, pm, f =  extract_CROCO('croco_avg_Y'+str(Y_min)+'M'+str(1)+'.nc',file_base=file_base,var='temp', 
                                 level = 59, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[1:7]

TIME = np.empty((0))
for i in range(Y_min, Y_max+1):
    for j in range(1,12+1):
        file_name = 'croco_avg_Y'+str(i)+'M'+str(j)+'.nc'
        if file_name != 'croco_avg_Y2018M12.nc':
            print("Processing "  + file_name)
            nc_file = file_base + file_name
            root = Dataset(nc_file)
            time_CROCO = root.variables['time'][:]            
            TIME = np.append(TIME,time_CROCO)
        else:
            print("No Data")
            
# Any repeated times, remove them

uniqueValues, indicesList = np.unique(TIME, return_index=True)
TIME = TIME[indicesList]    

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

# New mask for ohhsore region

g_mask = h_CROCO.copy()
g_mask[g_mask < shelf_depth] = float("NAN")
mask_offshelf = ~np.isnan(g_mask)
del g_mask 

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

###########################################
# Subdivide the SBUS into three regions
###########################################

# The SBUS will be divided into three: St Helena Bay [-32.82 -31.5]
# Mid-shelf [-31.5 -30.3] and Namaqua [-30.3 -28]

# St Helena Mask

cut_mask = SBUS_mask.copy()
CO_lat = -31.5   # Latitude of Cut min 
ind_CO = find_nearest(lat_CROCO,CO_lat)[1]
cut_mask[ind_CO:len(cut_mask),:] = float("NAN")
cut_mask = ~np.isnan(cut_mask)

SH_mask = cut_mask.copy()
del cut_mask

# Mid-shelf

cut_mask = SBUS_mask.copy()
CO_lat_min = -31.5   # Latitude of Cut min
CO_lat_max = -30.3   # Latitude of Cut max 
ind_CO_min = find_nearest(lat_CROCO,CO_lat_min)[1]
ind_CO_max = find_nearest(lat_CROCO,CO_lat_max)[1]
cut_mask[0:ind_CO_min,:] = float("NAN")
cut_mask[ind_CO_max:len(cut_mask),:] = float("NAN")
cut_mask = ~np.isnan(cut_mask)

MS_mask = cut_mask.copy()
del cut_mask
# Namaqua

cut_mask = SBUS_mask.copy()
CO_lat = -30.3   # Latitude of Cut min 
ind_CO = find_nearest(lat_CROCO,CO_lat)[1]
cut_mask[0:ind_CO,:] = float("NAN")
cut_mask = ~np.isnan(cut_mask)

NM_mask = cut_mask.copy()
del cut_mask

#%%

#######################
# LOAD the WIND File
#######################

f = open('/media/data/DIAGNOSTICS/FRONTS/WIND.pckl', 'rb')
WIND_COAST = pickle.load(f)
f.close()


#%%

# Calculate values                
            
###########################################
# Subdivide the SBUS into three regions
###########################################

WIND_up = np.nanmean(WIND_COAST,axis=0)
WIND_shelf = np.nanmean(WIND_COAST,axis=0)

# The SBUS will be divided into three: St Helena Bay [-32.82 -31.5]
# Mid-shelf [-31.5 -30.3] and Namaqua [-30.3 -28]

# St Helena Mask

CO_lat_min = -32.82   # Latitude of Cut min
CO_lat_max = -31.5 
ind_CO_min = find_nearest(lat_CROCO,CO_lat_min)[1]
ind_CO_max = find_nearest(lat_CROCO,CO_lat_max)[1]
WIND_SH = WIND_COAST[ind_CO_min:ind_CO_max,:] 
WIND_SH = np.nanmean(WIND_SH,0)


# Mid-shelf

CO_lat_min = -31.5   # Latitude of Cut min
CO_lat_max = -30.3   # Latitude of Cut max 
ind_CO_min = find_nearest(lat_CROCO,CO_lat_min)[1]
ind_CO_max = find_nearest(lat_CROCO,CO_lat_max)[1]
WIND_MS = WIND_COAST[ind_CO_min:ind_CO_max,:]
WIND_MS = np.nanmean(WIND_MS,0)



# Namaqua

CO_lat_min = -30.3   # Latitude of Cut min
CO_lat_max = -28   # Latitude of Cut max 
ind_CO_min = find_nearest(lat_CROCO,CO_lat_min)[1]
ind_CO_max = find_nearest(lat_CROCO,CO_lat_max)[1]
WIND_NM = WIND_COAST[ind_CO_min:ind_CO_max,:]
WIND_NM = np.nanmean(WIND_NM,0)


#%%

#####################
# LOAD Data frame
#####################

f = open('/media/data/DIAGNOSTICS/FRONTS/ST_HELENA.pckl', 'rb')
ST_HELENA = pickle.load(f)
f.close()

f = open('/media/data/DIAGNOSTICS/FRONTS/MID_SHELF.pckl', 'rb')
MID_SHELF = pickle.load(f)
f.close()

f = open('/media/data/DIAGNOSTICS/FRONTS/NAMAQUA.pckl', 'rb')
NAMAQUA = pickle.load(f)
f.close()

# Add data to the dataframes

ST_HELENA['Shelf_break']['WIND'] = WIND_SH 
ST_HELENA['Upwelling']['WIND'] = WIND_SH


MID_SHELF['Shelf_break']['WIND'] = WIND_MS 
MID_SHELF['Upwelling']['WIND'] = WIND_MS


NAMAQUA['Shelf_break']['WIND'] = WIND_NM 
NAMAQUA['Upwelling']['WIND'] = WIND_NM


#%%

#################################
# STRATIFICATION
#################################
 
# We want to know how stratification affect frontal dynamics. Calculate the avergage frontal BVF for NAMAQUA,
# St Helena and Mid-shelf 

######################################################################
# STEP 1: Calculate BVF_shelf and BVF_up 
######################################################################

# Edge map 

f = open('/media/data/DIAGNOSTICS/FRONTS/EDGE.pckl', 'rb')
EDGE = pickle.load(f)
f.close()

# BVF map

f = open('/media/data/DIAGNOSTICS/FRONTS/BVF.pckl', 'rb')
BVF = pickle.load(f)
f.close()

# Create a once off location of the offshore extent of the SST
BVF_shelf = []
BVF_up = []
for i in range(0,len(TIME)):
    loc_data = np.where(EDGE[:,:,i]*up_mask == 1)
    BVF_up.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF
    del loc_data
    loc_data = np.where(EDGE[:,:,i]*sb_mask == 1)
    BVF_shelf.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF

## Clean up and calculate index

del BVF
del EDGE

#%% Quick plot of the correlations for the shelf-break and upwelling regions for
# 1) the entire time domain
# 2) Summer
# 3 Winter

f = open('/media/data/DIAGNOSTICS/FRONTS/SBUS_up.pckl', 'rb')
SBUS_up = pickle.load(f)
f.close()

f = open('/media/data/DIAGNOSTICS/FRONTS/SBUS_shelf.pckl', 'rb')
SBUS_shelf = pickle.load(f)
f.close()

####################################
# Add wind and stratification
####################################

SBUS_up['WIND'] = WIND_up
SBUS_up['BVF'] = BVF_up

SBUS_shelf['WIND'] = WIND_shelf
SBUS_shelf['BVF'] = BVF_shelf

#%% READ in the FTLE

import scipy.io
mat = scipy.io.loadmat('FTLE_TS.mat')

FTLE_TS_up = np.array(mat['FTLE_TS_up'])
FTLE_TS_shelf = np.array(mat['FTLE_TS_shf'])
FTLE_TS_all = np.array(mat['FTLE_TS_all'])

# Hack to create dates

from datetime import date, timedelta

start_dt = date(2014, 1, 1)
end_dt = date(2014, 12, 31)

# difference between current and previous date
delta = timedelta(days=1)

# store the dates between two dates in a list
FTLE_dates = []

while start_dt <= end_dt:
    # add current date to list by converting  it to iso format
    FTLE_dates.append(start_dt.isoformat())
    # increment start date by timedelta
    start_dt += delta

print('Dates between', start_dt, 'and', end_dt)

FTLE_dates = np.array(FTLE_dates)

# Create a data set

data = {'Date': FTLE_dates,
        'UP': np.squeeze(FTLE_TS_up),
        'SHELF': np.squeeze(FTLE_TS_shelf),
        'ALL': np.squeeze(FTLE_TS_all)}
FTLE_df = pd.DataFrame(data) 

# convert time_date col to datetime64 dtype
FTLE_df['Date'] = pd.to_datetime(FTLE_df['Date'], utc=True)
FTLE_df = FTLE_df.set_index('Date')

FTLE_df_clim = FTLE_df.groupby([FTLE_df.index.month]).mean() # CLIM


#%%

# Two subplots, first showing the wind and second showing the stratification
 
######################################
# WIND Along the coast 
######################################

# One figure showing the WIND along the coast for the 
# Upwelling region of the different regions on one plot

# Create a data-frame

from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection 

data = {'Date': DATE,
        'ST_HELENA': WIND_SH,
        'MID_SHELF': WIND_MS,
        'NAMAQUA': WIND_NM}
WIND_df = pd.DataFrame(data) 
WIND_df = WIND_df.set_index('Date')

WIND_df_clim = WIND_df.groupby([WIND_df.index.month]).mean() # CLIM
WIND_df_clim_STD = WIND_df.groupby([WIND_df.index.month]).std() # CLIM

NAME = ["ST_HELENA","MID_SHELF","NAMAQUA"]
MARKER = ["o","v","s"]
COLOR = ["deepskyblue","limegreen","coral"]
labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul',
          'Aug','Sep','Oct','Nov','Dec']

BVF_up_clim = SBUS_up['BVF'].groupby([SBUS_up.index.month]).mean()
BVF_shelf_clim = SBUS_shelf['BVF'].groupby([SBUS_shelf.index.month]).mean()

BVF_up_clim = np.log10(BVF_up_clim)
BVF_shelf_clim = np.log10(BVF_shelf_clim)

fig, ax = plt.subplots(1,3,figsize=(20,5))

ax[0].set_facecolor('lightgrey')

trans1 = Affine2D().translate(-0.2, 0.0) + ax[0].transData
trans2 = Affine2D().translate(0.0, 0.0) + ax[0].transData
trans3 = Affine2D().translate(0.2, 0.0) + ax[0].transData

TRANS = [trans1,trans2,trans3]

for i in range(0,len(NAME)):
    inxval = list(range(1, 13))
    ax[0].errorbar(inxval, WIND_df_clim[NAME[i]], marker=MARKER[i], 
                 markersize=14,
                 color= COLOR[i],
                 label=NAME[i],
                 yerr = WIND_df_clim_STD[NAME[i]],
                 fmt="",
                 ecolor = COLOR[i],
                 elinewidth=3,
                 transform = TRANS[i])

ax[0].set_xticks(inxval, labels, rotation=45)
ax[0].tick_params(axis='both', which='major', labelsize=16)
ax[0].set_ylabel("Wind-stress [N/$m^{2}$]",fontsize=18)
ax[0].set_title("(a)",fontsize=15,loc='left')
ax[0].legend(borderpad=1,loc='lower right')
ax[0].grid()

ax[1].set_facecolor('lightgrey')
ax[1].plot(inxval,BVF_up_clim,
           color='blue',
           marker='o',
           markersize='14',
           linestyle='--',
           label='Upwelling region')
ax[1].plot(inxval,BVF_shelf_clim,
           color='red',
           marker='^',
           markersize='14',
           linestyle='--',
           label='Shelf-break region')

ax[1].set_xticks(inxval, labels, rotation=45)
ax[1].tick_params(axis='both', which='major', labelsize=16)
ax[1].set_ylabel("log10 BVF [$s^{-2}$]",fontsize=18)
ax[1].set_title("(b)",fontsize=15,loc='left')
ax[1].legend(borderpad=1,loc='lower left')
ax[1].grid()

ax[2].set_facecolor('lightgrey')
ax[2].plot(inxval,FTLE_df_clim.UP,
           color='blue',
           marker='o',
           markersize='14',
           linestyle='--',
           label='Upwelling region')
ax[2].plot(inxval,FTLE_df_clim.SHELF,
           color='red',
           marker='^',
           markersize='14',
           linestyle='--',
           label='Shelf-break region')
ax[2].plot(inxval,FTLE_df_clim.ALL,
           color='black',
           marker='s',
           markersize='14',
           linestyle='--',
           label='Whole shelf')

ax[2].set_xticks(inxval, labels, rotation=45)
ax[2].tick_params(axis='both', which='major', labelsize=16)
ax[2].set_ylabel("FTLE [$days^{-1}$]",fontsize=18)
ax[2].set_title("(c)",fontsize=15,loc='left')
ax[2].legend(borderpad=1,loc='upper left')
ax[2].grid()

plt.tight_layout()

plt.show()

fig.savefig('WIND_STRESS_BVF_FTLE_clim', format='png', dpi=600)
plt.close()