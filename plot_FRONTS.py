#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 12:45:58 2022

@author: jono
"""

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
    mask_CROCO[mask_CROCO == 0 ] = float("NAN")
    var_CROCO = np.multiply(var_CROCO,mask_CROCO)
    h_CROCO = np.multiply(h_CROCO,mask_CROCO)
    root.close()
    return var_CROCO, lon_CROCO, lat_CROCO, h_CROCO, pn, pm, time_CROCO

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

lon_CROCO, lat_CROCO, h_CROCO, pn, pm =  extract_CROCO('croco_avg_Y'+str(Y_min)+'M'+str(1)+'.nc',file_base=file_base,var='temp', 
                                 level = 59, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[1:6]

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

######################
# ALTERNTIVE
######################
# STILL WORK in PROGRESS
# Define a latitude band and find fronts in it
#latitude = -32
#ind_lat = find_nearest(lat_CROCO,latitude)[1]
#SH_mask = SBUS_mask.copy()[ind_lat-10:ind_lat+10,:]

#%%

#######################
# LOAD in the data
#######################

f = open('GRADIENT.pckl', 'rb')
GRADIENT = pickle.load(f)
f.close()

f = open('EDGE.pckl', 'rb')
EDGE = pickle.load(f)
f.close()

f = open('W.pckl', 'rb')
W = pickle.load(f)
f.close()

# Quick plot

# Check contours in relation to masks

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
map.pcolormesh(x,y,W[:,:,95]*up_mask*NM_mask,vmin=-3e-5,vmax=3e-5,cmap=cmo.balance, shading='auto')
pcm = map.contour(X, Y, EDGE[:,:,95]*up_mask*NM_mask, linewidths=1, linestyles='dashed', colors='green')
map.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
plt.colorbar()

#%%

def angleFromCoordinate(lat1, long1, lat2, long2):
    import math
    dLon = (long2 - long1)

    y = math.sin(dLon) * math.cos(lat2)
    x = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dLon)

    brng = math.atan2(y, x)

    brng = math.degrees(brng)
    brng = (brng + 360) % 360
    brng = 360 - brng # count degrees clockwise - remove to make counter-clockwise

    return brng

ang_coast = 22
    
######################
# MAIN 
######################

# Need to ammend the length of them and calculate angle
# SHELF-BREAK FRONT

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
        llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)

DIS_shelf = []
LEN_shelf = []
NUM_shelf = []
GRAD_shelf = []
W_shelf = []
LONG_shelf = []
ACROSS_shelf = []
for i in range(0,len(TIME)): 
    loc_data = np.where(EDGE[:,:,i]*sb_mask*MS_mask == 1)
    W_shelf.append(np.nanmean(W[loc_data[0],loc_data[1],i]))  # Vertical velocity
    GRAD_shelf.append(np.nanmean(GRADIENT.data[loc_data[0],loc_data[1],i]))  # Pixel gradient
    pcm = map.contour(X, Y, EDGE[:,:,i]*sb_mask*MS_mask, linewidths=0.1, linestyles='solid', colors='black')
    pcm = pcm.allsegs[0]
    NUM_shelf.append(len(pcm))  # Number of fronts detected 
    dis_mean = []
    len_mean = []
    long = []
    across = []
    for k in range(0,len(pcm)):
        cc = pcm[k]
        Xloc, Yloc = contour_loc(cc,x,y)
        # New Segment
        
        cc, indices = np.unique(cc, return_index=True, axis=0)
        Xloc = Xloc[indices]
        Yloc = Yloc[indices]
        cc = cc[np.argsort(cc[:, 1])]
        if len(cc) == 0:
            print("EMPTY")
        ANGLE = []
        res = []
        for xx in range(0,len(cc)-1):
            ang = angleFromCoordinate(cc[xx,1], cc[xx,0], cc[xx+1,1], cc[xx+1,0])
            ANGLE.append(ang) 
         
        for idx in range(0, len(ANGLE)):
            if ANGLE[idx] < ang_coast + 45 or ANGLE[idx] > (360 - 45)+ang_coast:
                res.append(idx)
            elif ANGLE[idx] > (90+45)+ang_coast and ANGLE[idx] < (180+45)+ang_coast:
                res.append(idx)
         
        if len(res) >= 0.5*len(ANGLE):
            long.append(1)
        else:
            across.append(1)       
         
        # New segment         
        len_mean.append(len(cc))
        DIS = np.empty((0))
        for j in range(0,len(Yloc)):
            point = (lat_CROCO[Yloc[j].astype(int)].astype(float).tolist(),
                     lon_CROCO[Xloc[j].astype(int)].astype(float).tolist())
            coast = (lat_CROCO[Yloc[j].astype(int)].astype(float).tolist(),
                    lon_CROCO[ind_coast_X[np.argwhere(ind_coast_Y == Yloc[j])]
                                 .astype(int)].data.item())
            dis = distance.distance(point,coast).km
            DIS = np.append(DIS,dis)
        dis_mean.append(np.mean(DIS))
    LONG_shelf.append(np.sum(long))
    ACROSS_shelf.append(np.sum(across))
    DIS_shelf.append(np.mean(dis_mean)) # Displacement of fronts
    LEN_shelf.append(np.mean(len_mean)) # Length of fronts
    
# UPWELLING FRONT

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
        llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)

DIS_up = []
LEN_up = []
NUM_up = []
GRAD_up = []
W_up = []
LONG_up = []
ACROSS_up = []
for i in range(0,len(TIME)): 
    loc_data = np.where(EDGE[:,:,i]*up_mask*MS_mask == 1)
    W_up.append(np.nanmean(W[loc_data[0],loc_data[1],i]))  # Vertical velocity
    GRAD_up.append(np.nanmean(GRADIENT.data[loc_data[0],loc_data[1],i]))  # Pixel gradient
    pcm = map.contour(X, Y, EDGE[:,:,i]*up_mask*MS_mask, linewidths=0.1, linestyles='solid', colors='black')
    pcm = pcm.allsegs[0]
    NUM_up.append(len(pcm))  # Number of fronts detected 
    dis_mean = []
    len_mean = [] 
    long = []
    across = []
    for k in range(0,len(pcm)):
        cc = pcm[k]
        Xloc, Yloc = contour_loc(cc,x,y)
        # New Segment
        
        cc, indices = np.unique(cc, return_index=True, axis=0)
        Xloc = Xloc[indices]
        Yloc = Yloc[indices]
        cc = cc[np.argsort(cc[:, 1])]
        if len(cc) == 0:
            print("EMPTY")
        ANGLE = []
        res = []
        for xx in range(0,len(cc)-1):
            ang = angleFromCoordinate(cc[xx,1], cc[xx,0], cc[xx+1,1], cc[xx+1,0])
            ANGLE.append(ang) 
         
        for idx in range(0, len(ANGLE)):
            if ANGLE[idx] < ang_coast + 45 or ANGLE[idx] > (360 - 45)+ang_coast:
                res.append(idx)
            elif ANGLE[idx] > (90+45)+ang_coast and ANGLE[idx] < (180+45)+ang_coast:
                res.append(idx)
         
        if len(res) >= 0.5*len(ANGLE):
            long.append(1)
        else:
            across.append(1)       
         
        # New segment         
        len_mean.append(len(cc))
        DIS = np.empty((0))
        for j in range(0,len(Yloc)):
            point = (lat_CROCO[Yloc[j].astype(int)].astype(float).tolist(),
                     lon_CROCO[Xloc[j].astype(int)].astype(float).tolist())
            coast = (lat_CROCO[Yloc[j].astype(int)].astype(float).tolist(),
                    lon_CROCO[ind_coast_X[np.argwhere(ind_coast_Y == Yloc[j])]
                                 .astype(int)].data.item())
            dis = distance.distance(point,coast).km
            DIS = np.append(DIS,dis)
        dis_mean.append(np.mean(DIS))
    LONG_up.append(np.sum(long))
    ACROSS_up.append(np.sum(across))    
    DIS_up.append(np.mean(dis_mean)) # Displacement of fronts
    LEN_up.append(np.mean(len_mean)) # Length of fronts

#%%
##################    
# CREATE dataframe
##################
    
################################
# SHELF-BREAK data frame
################################    
    
import pandas as pd    

data = {'Date': DATE,
        'DIS': DIS_shelf,
        'LEN': LEN_shelf,
        'NUM': NUM_shelf,
        'NUM_ACROSS': ACROSS_shelf,
        'NUM_LONG': LONG_shelf,
        'GRAD': GRAD_shelf,
        'W': W_shelf}
df_shelf = pd.DataFrame(data) 
df_shelf = df_shelf.set_index('Date')

###################################
# UPWELLING data frame
###################################

data = {'Date': DATE,
        'DIS': DIS_up,
        'LEN': LEN_up,
        'NUM': NUM_up,
        'NUM_ACROSS': ACROSS_up,
        'NUM_LONG': LONG_up,
        'GRAD': GRAD_up,
        'W': W_up}
df_up = pd.DataFrame(data) 
df_up = df_up.set_index('Date')   

#%%

###################
# SAVE dataframes
###################

MID_SHELF =	{
  "Upwelling": df_up,
  "Shelf_break": df_shelf
}

f = open('MID_SHELF.pckl', 'wb')
pickle.dump(MID_SHELF, f)
f.close()

#%%

#####################
# LOAD Data frame
#####################

f = open('NAMAQUA.pckl', 'rb')
NAMAQUA = pickle.load(f)
f.close()

df_shelf = NAMAQUA['Shelf_break']
df_up = NAMAQUA['Upwelling']

#%%
############################################
# Plot combining displacmenet and gradient
############################################

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection 

# Compute the monthly means for each variable

data_columns = ['DIS', 'LEN', 'NUM','GRAD','W']

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
y = np.array(day3_shelf.DIS)
c = np.array(day3_shelf.GRAD)
             
axes[0].set_title("(a) Shelf-break region",loc='left', fontsize=16)
#convert dates to numbers first
inxval = mdates.date2num(dates)
years = mdates.YearLocator()
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="viridis", linewidth=1)

# Add the monthly values
axes[0].scatter(x=mon_shelf.index, y=mon_shelf['DIS'],s=90,c=mon_shelf['GRAD'],cmap='viridis',
                vmin=0, vmax=0.2)
# Add line to show mean
axes[0].axhline(y=np.nanmean(df_shelf.DIS), xmin=0, xmax=1, linewidth=2, color='r')
# set color to date values
lc.set_array(c)
line = axes[0].add_collection(lc)
line.set_clim(vmin=0, vmax=0.2)
axes[0].tick_params(axis="y", labelsize=14) 
axes[0].grid()
axes[0].set_ylabel("Distance [km]",fontsize=16)
axes[0].set_ylim([80, 180])

# Upwelling
dates = day3_up.index
y = np.array(day3_up.DIS)
c = np.array(day3_up.GRAD)

axes[1].set_title("(b) Upwelling region",loc='left', fontsize=16)
#convert dates to numbers first
inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="viridis", linewidth=1)
# Add the monthly values
axes[1].scatter(x=mon_up.index, y=mon_up['DIS'],s=90,c=mon_up['GRAD'],cmap='viridis',
                vmin =0, vmax=0.2)
# Add line to show mean
axes[1].axhline(y=np.nanmean(df_up.DIS), xmin=0, xmax=1, linewidth=2, color='r')
# set color to date values
lc.set_array(c)
line = axes[1].add_collection(lc)
axes[1].xaxis.set_major_locator(years)
line.set_clim(vmin=0, vmax=0.2)
axes[1].grid()
axes[1].set_xlabel("DATE",fontsize=16)
axes[1].set_ylabel("Distance [km]",fontsize=16)
axes[1].set_ylim([30, 100])
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a
plt.colorbar(line, ax=axes.ravel().tolist()).set_label(label="Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",size=13)
axes[1].xaxis_date()
plt.show()

fig.savefig('NAM_Displacement', format='png', dpi=600)

#%%
##########################################################
# GRADINET vs vertical velocity
###########################################################

##########################
# MAIN FIGURE
##########################

fig, axes = plt.subplots(2,1,sharex=True, figsize=(18,7))

# Shelf-break
dates = day3_shelf.index
y = np.array(day3_shelf.GRAD)
c = np.array(day3_shelf.W)

axes[0].set_title("(a) Shelf-break region",loc='left', fontsize=16)
#convert dates to numbers first
inxval = mdates.date2num(dates)
years = mdates.YearLocator()
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="RdBu_r", linewidth=1)
# Add the monthly values
axes[0].scatter(x=mon_shelf.index, y=mon_shelf['GRAD'],s=90,c=mon_shelf['W'],cmap='RdBu_r',
                vmin=-2e-5, vmax=2e-5)
# Add a mean line
axes[0].axhline(y=np.nanmean(df_shelf.GRAD), xmin=0, xmax=1, linewidth=2, color='b')
# set color to date values
lc.set_array(c)
line = axes[0].add_collection(lc)
line.set_clim(vmin=-2e-5, vmax=2e-5)
axes[0].tick_params(axis="y", labelsize=14) 
axes[0].grid()
axes[0].set_ylabel("Gradient",fontsize=16)

# Upwelling
dates = day3_up.index
y = np.array(day3_up.GRAD)
c = np.array(day3_up.W)

axes[1].set_title("(b) Upwelling region",loc='left', fontsize=16)
#convert dates to numbers first
inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="RdBu_r", linewidth=1)
# Add the monthly values
axes[1].scatter(x=mon_up.index, y=mon_up['GRAD'],s=90,c=mon_up['W'],cmap='RdBu_r',
                vmin=-2e-5, vmax=2e-5)
# Add a mean line
axes[1].axhline(y=np.nanmean(df_up.GRAD), xmin=0, xmax=1, linewidth=2, color='b')
# set color to date values
lc.set_array(c)
line = axes[1].add_collection(lc)
line.set_clim(vmin=-2e-5, vmax=2e-5)
axes[1].grid()
axes[1].xaxis.set_major_locator(years)
axes[1].set_xlabel("DATE",fontsize=16)
axes[1].set_ylabel("Gradient",fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a

plt.colorbar(line, ax=axes.ravel().tolist(), label='Vertical velocities [m/s]' )
axes[1].xaxis_date()
plt.show()

fig.savefig('NAM_Wgrad', format='png', dpi=600)
plt.close()

#%%

###############################################
# CLIMATOLOGY FIGURES
###############################################

# The 3-day time-series are usueful and show the huge range of interannual
# variability, but we are interested in the seasonality of the fronts

# Construct a subplot with the displacement vs gradiet and upwelling vs gradient
# for each region

df_clim_shelf = df_shelf.groupby([df_shelf.index.month]).mean() # CLIM MEAN
df_clim_shelf_STD = df_shelf.groupby([df_shelf.index.month]).std()  # CLIM standard deviation
 
df_clim_up = df_up.groupby([df_up.index.month]).mean() # CLIM
df_clim_up_STD = df_up.groupby([df_up.index.month]).std() # CLIM standard deviation

fig, axes = plt.subplots(2,2,sharex=True, figsize=(18,7))

labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul',
          'Aug','Sep','Oct','Nov','Dec']

# Dsiplacement plot

# Shelf-break
dates = df_clim_shelf.index
y = np.array(df_clim_shelf.DIS)
c = np.array(df_clim_shelf.GRAD)
             
axes[0,0].set_title("(a) Shelf-break region",loc='left', fontsize=16)
#convert dates to numbers first
inxval = list(range(1, 13))
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="viridis", linewidth=1)

# Shade in the displacement variability

y1 = np.array(df_clim_shelf.DIS) - np.array(df_clim_shelf_STD.DIS)
y2 = np.array(df_clim_shelf.DIS) + np.array(df_clim_shelf_STD.DIS)

axes[0,0].fill_between(inxval, y1, y2,
                       alpha=0.3,
                       edgecolor='black',
                       facecolor='grey')       # Transparency of the fill

# Add the monthly values
axes[0,0].scatter(x=df_clim_shelf.index, y=df_clim_shelf['DIS'],s=90,c=df_clim_shelf['GRAD'],cmap='viridis',
                vmin=0, vmax=0.2)
# Add line to show mean
axes[0,0].axhline(y=np.nanmean(df_clim_shelf.DIS), xmin=0, xmax=1, linewidth=2, color='r')
# set color to date values
lc.set_array(c)
line = axes[0,0].add_collection(lc)
line.set_clim(vmin=0, vmax=0.2)
axes[0,0].tick_params(axis="y", labelsize=14) 
axes[0,0].grid()
axes[0,0].set_ylabel("Distance [km]",fontsize=16)
axes[0,0].set_ylim([80, 180])

# Upwelling
dates = df_clim_up.index
y = np.array(df_clim_up.DIS)
c = np.array(df_clim_up.GRAD)

axes[0,1].set_title("(b) Upwelling region",loc='left', fontsize=16)
#convert dates to numbers first
#inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="viridis", linewidth=1)
# Shade in the displacement variability
y1 = np.array(df_clim_up.DIS) - np.array(df_clim_up_STD.DIS)
y2 = np.array(df_clim_up.DIS) + np.array(df_clim_up_STD.DIS)

axes[0,1].fill_between(inxval, y1, y2,
                       alpha=0.3,
                       edgecolor='black',
                       facecolor='grey')       # Transparency of the fill
# Add the monthly values
axes[0,1].scatter(x=df_clim_up.index, y=df_clim_up['DIS'],s=90,c=df_clim_up['GRAD'],cmap='viridis',
                vmin =0, vmax=0.2)
# Add line to show mean
axes[0,1].axhline(y=np.nanmean(df_clim_up.DIS), xmin=0, xmax=1, linewidth=2, color='r')
# set color to date values
lc.set_array(c)
line = axes[0,1].add_collection(lc)
line.set_clim(vmin=0, vmax=0.2)
axes[0,1].tick_params(axis="y", labelsize=14)
axes[0,1].grid()
axes[0,1].set_ylabel("Distance [km]",fontsize=16)
axes[0,1].set_ylim([30, 100])
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a
plt.colorbar(line, ax=axes[0,:], label="Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]")

# Gradient plot

# Shelf-break
dates = df_clim_shelf.index
y = np.array(df_clim_shelf.GRAD)
c = np.array(df_clim_shelf.W)

axes[1,0].set_title("(c) Shelf-break region",loc='left', fontsize=16)
#convert dates to numbers first
#inxval = mdates.date2num(dates)
#years = mdates.YearLocator()
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="RdBu_r", linewidth=1)
# Shade in the displacement variability

y1 = np.array(df_clim_shelf.GRAD) - np.array(df_clim_shelf_STD.GRAD)
y2 = np.array(df_clim_shelf.GRAD) + np.array(df_clim_shelf_STD.GRAD)

axes[1,0].fill_between(inxval, y1, y2,
                       alpha=0.3,
                       edgecolor='black',
                       facecolor='grey')       # Transparency of the fill
# Add the monthly values
axes[1,0].scatter(x=df_clim_shelf.index, y=df_clim_shelf['GRAD'],s=90,c=df_clim_shelf['W'],cmap='RdBu_r',
                vmin=-5e-6, vmax=5e-6)
# Add a mean line
axes[1,0].axhline(y=np.nanmean(df_clim_shelf.GRAD), xmin=0, xmax=1, linewidth=2, color='b')
# set color to date values
lc.set_array(c)
line = axes[1,0].add_collection(lc)
line.set_clim(vmin=-5e-6, vmax=5e-6) 
axes[1,0].grid()
axes[1,0].tick_params(axis="x", labelsize=14)
axes[1,0].tick_params(axis="y", labelsize=14)
axes[1,0].xaxis.set_ticks(inxval)
axes[1,0].xaxis.set_ticklabels(labels,rotation=45,ha = 'right')
axes[1,0].set_ylabel("Gradient",fontsize=16)

# Upwelling
dates = df_clim_up.index
y = np.array(df_clim_up.GRAD)
c = np.array(df_clim_up.W)

axes[1,1].set_title("(d) Upwelling region",loc='left', fontsize=16)
#convert dates to numbers first
#inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="RdBu_r", linewidth=1)
# Shade in the displacement variability

y1 = np.array(df_clim_up.GRAD) - np.array(df_clim_up_STD.GRAD)
y2 = np.array(df_clim_up.GRAD) + np.array(df_clim_up_STD.GRAD)

axes[1,1].fill_between(inxval, y1, y2,
                       alpha=0.3,
                       edgecolor='black',
                       facecolor='grey')       # Transparency of the fill

# Add the monthly values
axes[1,1].scatter(x=df_clim_up.index, y=df_clim_up['GRAD'],s=90,c=df_clim_up['W'],cmap='RdBu_r',
                vmin=-5e-6, vmax=5e-6)
# Add a mean line
axes[1,1].axhline(y=np.nanmean(df_clim_up.GRAD), xmin=0, xmax=1, linewidth=2, color='b')
# set color to date values
lc.set_array(c)
line = axes[1,1].add_collection(lc)
line.set_clim(vmin=-5e-6, vmax=5e-6)
axes[1,1].grid()
axes[1,1].set_ylabel("Gradient",fontsize=16)
axes[1,1].xaxis.set_ticks(inxval)
axes[1,1].xaxis.set_ticklabels(labels,rotation=45,ha = 'right')
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a

plt.colorbar(line, ax=axes[1, :], label='Vertical velocities [m/s]' )
plt.show()

fig.savefig('NAMAQUA_TS', format='png', dpi=600)
plt.close()

#%%

#############################
# Front probability plot
#############################

###################################
# LOAD BVF
###################################

f = open('BVF.pckl', 'rb')
BVF = pickle.load(f)
f.close()

#####################################

prob_sum = np.nansum(EDGE[:,:,ind_summer.astype(int)],2)/len(TIME[ind_summer.astype(int)])
prob_win = np.nansum(EDGE[:,:,ind_winter.astype(int)],2)/len(TIME[ind_winter.astype(int)])

grad_sum = np.nanmean(GRADIENT[:,:,ind_summer.astype(int)],2)
grad_win = np.nanmean(GRADIENT[:,:,ind_winter.astype(int)],2)

w_sum = np.nanmean(W[:,:,ind_summer.astype(int)],2)
w_win = np.nanmean(W[:,:,ind_winter.astype(int)],2)

bvf_sum = np.nanmean(BVF[:,:,ind_summer.astype(int)],2)
bvf_win = np.nanmean(BVF[:,:,ind_winter.astype(int)],2)

bvf_sum = np.log10(bvf_sum)
bvf_win = np.log10(bvf_win)

prob_sum[prob_sum == 0] = float("NAN")
prob_win[prob_win == 0] = float("NAN")

fig, axes = plt.subplots(4,2,sharex=False,constrained_layout = True)

# Probability

# Summer
axes[0,0].set_title("(a) Summer [DJF]", loc='left', size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0,0], resolution='l')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[0,0].contourf(X,Y,prob_sum.data,vmin=0,vmax=0.3,cmap=cmo.dense)
axes[0,0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

# Winter
axes[0,1].set_title("(b) Winter [JJA]",loc='left', size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0,1], resolution='l')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[0,1].contourf(X,Y,prob_win.data,vmin=0,vmax=0.3,cmap=cmo.dense)
axes[0,1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

plt.colorbar(im, ax=axes[0, :], label='Probability' )

# Gradient

# Summer
axes[1,0].set_title("(c)",loc='left', size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[1,0], resolution='l')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[1,0].pcolor(X,Y,grad_sum.data,vmin=0,vmax=0.2,cmap=cmo.thermal)
axes[1,0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

# Winter
axes[1,1].set_title("(d)",loc='left',size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[1,1], resolution='l')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[1,1].pcolor(X,Y,grad_win.data,vmin=0,vmax=0.2,cmap=cmo.thermal)
axes[1,1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

plt.colorbar(im, ax=axes[1, :], label="Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]")

# Vertical Velocity

# Summer
axes[2,0].set_title("(e)",loc='left', size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[2,0], resolution='l')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[2,0].pcolor(X,Y,w_sum.data,vmin=-2e-5,vmax=2e-5,cmap='RdBu_r', shading='nearest')
axes[2,0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

# Winter
axes[2,1].set_title("(f)",loc='left',size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[2,1], resolution='l')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[2,1].pcolor(X,Y,w_win.data,vmin=-2e-5,vmax=2e-5,cmap='RdBu_r', shading='nearest')
axes[2,1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

plt.colorbar(im, ax=axes[2, :], label='Vertical velocities [m/s]')

plt.tick_params(axis='both', which='major', labelsize=12)

# BVF

axes[3,0].set_title("(g)",loc='left', size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[3,0], resolution='l')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[3,0].pcolor(X,Y,bvf_sum.data,vmin=-6,vmax=-2,cmap=cmo.amp, shading='nearest')
axes[3,0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

# Winter
axes[3,1].set_title("(h)",loc='left',size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[3,1], resolution='l')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[3,1].pcolor(X,Y,bvf_win.data,vmin=-6,vmax=-2,cmap=cmo.amp, shading='nearest')
axes[3,1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

plt.colorbar(im, ax=axes[3, :], label='log10 BVF [$s^{-2}$]')  

plt.tick_params(axis='both', which='major', labelsize=12)

#cbar = fig.colorbar(im, shrink = 0.5, ax=axes.ravel().tolist(),
#                    spacing='uniform',
##                    orientation='vertical'
#                   )
#cbar.set_label("Probability",size=12)
fig.set_size_inches(7,12)
plt.show()
fig.savefig('SBUS_prob.png', format='png', dpi=600)
plt.close()