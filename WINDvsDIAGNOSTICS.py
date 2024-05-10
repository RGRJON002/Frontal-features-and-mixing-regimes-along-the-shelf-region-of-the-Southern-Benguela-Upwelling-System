#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 13:57:09 2022

@author: jono
"""

# Compare the alongshore wind intensity at different spatial levels for the frontal
# diagnostics. We also want to look at frontal diagnostics under certain wind conditions
# To assess this, we need to construct some global histograms of wind and from that
# identify wind and then count events.

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
    
#######################
# SAVE the wind data
#######################

#f = open('WIND.pckl', 'wb')
#pickle.dump(WIND_COAST, f)
#f.close()    

#######################
# LOAD the WIND File
#######################

f = open('WIND.pckl', 'rb')
WIND_COAST = pickle.load(f)
f.close()

#%%

########################################
# GLOBAL HISTOGRAM of WIND
########################################

fig, axes = plt.subplots(2,1,figsize=(10,10))

# Summer

axes[0].set_title('Frequency plot Summer',size=13)
axes[0].set_ylabel('Normalized Cumulative Frequency',size=12)
#axes[0].set_xlabel("Wind-Stress [N/m^2]",size=12)
gr = axes[0].hist(np.nanmean(WIND_COAST[:,ind_summer.astype(int)],axis=0), bins=100, density=1, histtype='step', cumulative=True)

val_0_8 = find_nearest(gr[0],0.8)
val_0_9 = find_nearest(gr[0],0.9) 

axes[0].hlines(y=val_0_8[0], xmin=0, xmax=gr[1][val_0_8[1].astype(int)], linewidth=2, color='r',
               label = ('Lower Threshold='+str(round(gr[1][val_0_8[1].astype(int)],3))))
axes[0].vlines(x=gr[1][val_0_8[1].astype(int)],ymin=0, ymax = gr[0][val_0_8[1].astype(int)],
               linewidth=2, color='r')
axes[0].hlines(y=val_0_9[0], xmin=0, xmax=gr[1][val_0_9[1].astype(int)], linewidth=2, color='g',
               label = ('Upper Threshold='+str(round(gr[1][val_0_9[1].astype(int)],3))))
axes[0].vlines(x=gr[1][val_0_9[1].astype(int)],ymin=0, ymax = gr[0][val_0_9[1].astype(int)],linewidth=2, color='g')

leg = axes[0].legend(loc = 'upper left')
axes[0].tick_params(axis='both', which='major', labelsize=12)

# Winter

axes[1].set_title('Frequency plot Winter',size=13)
axes[1].set_ylabel('Normalized Cumulative Frequency',size=12)
axes[1].set_xlabel("Wind-Stress [N/m^2]",size=12)
gr = axes[1].hist(np.nanmean(WIND_COAST[:,ind_winter.astype(int)],axis=0), bins=100, density=1, histtype='step', cumulative=True)

val_0_8 = find_nearest(gr[0],0.8)
val_0_9 = find_nearest(gr[0],0.9) 

axes[1].hlines(y=val_0_8[0], xmin=0, xmax=gr[1][val_0_8[1].astype(int)], linewidth=2, color='r',
               label = ('Lower Threshold='+str(round(gr[1][val_0_8[1].astype(int)],3))))
axes[1].vlines(x=gr[1][val_0_8[1].astype(int)],ymin=0, ymax = gr[0][val_0_8[1].astype(int)],
               linewidth=2, color='r')
axes[1].hlines(y=val_0_9[0], xmin=0, xmax=gr[1][val_0_9[1].astype(int)], linewidth=2, color='g',
               label = ('Upper Threshold='+str(round(gr[1][val_0_9[1].astype(int)],3))))
axes[1].vlines(x=gr[1][val_0_9[1].astype(int)],ymin=0, ymax = gr[0][val_0_9[1].astype(int)],linewidth=2, color='g')

leg = axes[1].legend(loc = 'upper left')
axes[1].tick_params(axis='both', which='major', labelsize=12)


plt.tick_params(axis='both', which='major', labelsize=12)
plt.show()

#per90 = np.percentile(np.nanmean(WIND_COAST[:,ind_summer.astype(int)],axis=0), 90)
#fig.savefig('WIND_HIST_sumwin.png', format='png', dpi=600)

plt.close()

#%%

# Calculate values                
            
###########################################
# Subdivide the SBUS into three regions
###########################################

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

f = open('ST_HELENA.pckl', 'rb')
ST_HELENA = pickle.load(f)
f.close()

f = open('MID_SHELF.pckl', 'rb')
MID_SHELF = pickle.load(f)
f.close()

f = open('NAMAQUA.pckl', 'rb')
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

#fig = plt.figure(figsize=(9,8))
fig, ax = plt.subplots(figsize=(9,8))
ax.set_facecolor('lightgrey')

trans1 = Affine2D().translate(-0.2, 0.0) + ax.transData
trans2 = Affine2D().translate(0.0, 0.0) + ax.transData
trans3 = Affine2D().translate(0.2, 0.0) + ax.transData

TRANS = [trans1,trans2,trans3]

for i in range(0,len(NAME)):
    inxval = list(range(1, 13))
    plt.errorbar(inxval, WIND_df_clim[NAME[i]], marker=MARKER[i], 
                 markersize=14,
                 color= COLOR[i],
                 label=NAME[i],
                 yerr = WIND_df_clim_STD[NAME[i]],
                 fmt="",
                 ecolor = COLOR[i],
                 elinewidth=3,
                 transform = TRANS[i])
    
plt.xticks(inxval, labels, rotation=45)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel("Wind-stress [N/$m^{2}$]",fontsize=18)
plt.legend(borderpad=1)
plt.grid()
plt.show()

#fig.savefig('WIND_STRESS_clim', format='png', dpi=600)
plt.close()
#%% Quick plot of the correlations for the shelf-break and upwelling regions for
# 1) the entire time domain
# 2) Summer
# 3 Winter

fields = ['WIND','W','DIS','LEN','GRAD','NUM']

import seaborn as sn

from scipy.stats import pearsonr
import pandas as pd

def calculate_pvalues(df):
    df = df.dropna()._get_numeric_data()
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            pvalues[r][c] = round(pearsonr(df[r], df[c])[1], 5)
    return pvalues

###########################################
# WHOLE YEAR
###########################################

fig, ax = plt.subplots(3, 2,figsize=(15,12))

fig.suptitle('Shelf-break and Upwelling', fontsize=16)

ax[0,0].set_title('(a) St Helena',size=13)
corrMatrix = ST_HELENA['Shelf_break'][fields].groupby(ST_HELENA['Shelf_break'].index.month).transform(lambda x: x-x.mean()).corr()
PMatrix = calculate_pvalues(ST_HELENA['Shelf_break'][fields].groupby(ST_HELENA['Shelf_break'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[0,1].set_title('(b) St Helena',size=13)
corrMatrix = ST_HELENA['Upwelling'][fields].groupby(ST_HELENA['Upwelling'].index.month).transform(lambda x: x-x.mean()).corr()
PMatrix = calculate_pvalues(ST_HELENA['Upwelling'][fields].groupby(ST_HELENA['Upwelling'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,1],vmin=-0.6, vmax = 0.6, cbar = False)

ax[1,0].set_title('(c) Mid-shelf',size=13)
corrMatrix = MID_SHELF['Shelf_break'][fields].groupby(MID_SHELF['Shelf_break'].index.month).transform(lambda x: x-x.mean()).corr()
PMatrix = calculate_pvalues(MID_SHELF['Shelf_break'][fields].groupby(MID_SHELF['Shelf_break'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[1,1].set_title('(d) Mid-shelf',size=13)
corrMatrix = MID_SHELF['Upwelling'][fields].groupby(MID_SHELF['Upwelling'].index.month).transform(lambda x: x-x.mean()).corr()
PMatrix = calculate_pvalues(MID_SHELF['Upwelling'][fields].groupby(MID_SHELF['Upwelling'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,1],vmin=-0.6, vmax = 0.6, cbar = False)

ax[2,0].set_title('(e) Namaqua',size=13)
corrMatrix = NAMAQUA['Shelf_break'][fields].groupby(NAMAQUA['Shelf_break'].index.month).transform(lambda x: x-x.mean()).corr()
PMatrix = calculate_pvalues(NAMAQUA['Shelf_break'][fields].groupby(NAMAQUA['Shelf_break'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[2,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[2,1].set_title('(f) Namaqua',size=13)
corrMatrix = NAMAQUA['Upwelling'][fields].groupby(NAMAQUA['Upwelling'].index.month).transform(lambda x: x-x.mean()).corr()
PMatrix = calculate_pvalues(NAMAQUA['Upwelling'][fields].groupby(NAMAQUA['Upwelling'].index.month).transform(lambda x: x-x.mean()).corr())
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[2,1],vmin=-0.6, vmax = 0.6, cbar = False)

plt.show()

#fig.savefig('CORR_ALLYEAR', format='png', dpi=600)
plt.close()

##################################
# Summer
##################################

fig, ax = plt.subplots(3, 2,figsize=(15,12))

fig.suptitle('Shelf-break and Upwelling: Summer [DJF]', fontsize=16)

ax[0,0].set_title('(a) St Helena',size=13)
corrMatrix = ST_HELENA['Shelf_break'][fields].groupby(ST_HELENA['Shelf_break'].index.month).transform(lambda x: x-x.mean()).iloc[ind_summer.astype(int)].corr()
PMatrix = calculate_pvalues(ST_HELENA['Shelf_break'][fields].groupby(ST_HELENA['Shelf_break'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[0,1].set_title('(b) St Helena',size=13)
corrMatrix = ST_HELENA['Upwelling'][fields].groupby(ST_HELENA['Upwelling'].index.month).transform(lambda x: x-x.mean()).iloc[ind_summer.astype(int)].corr()
PMatrix = calculate_pvalues(ST_HELENA['Upwelling'][fields].groupby(ST_HELENA['Upwelling'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,1],vmin=-0.6, vmax = 0.6, cbar = False)

ax[1,0].set_title('(c) Mid-shelf',size=13)
corrMatrix = MID_SHELF['Shelf_break'][fields].groupby(MID_SHELF['Shelf_break'].index.month).transform(lambda x: x-x.mean()).iloc[ind_summer.astype(int)].corr()
PMatrix = calculate_pvalues(MID_SHELF['Shelf_break'][fields].groupby(MID_SHELF['Shelf_break'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[1,1].set_title('(d) Mid-shelf',size=13)
corrMatrix = MID_SHELF['Upwelling'][fields].groupby(MID_SHELF['Upwelling'].index.month).transform(lambda x: x-x.mean()).iloc[ind_summer.astype(int)].corr()
PMatrix = calculate_pvalues(MID_SHELF['Upwelling'][fields].groupby(MID_SHELF['Upwelling'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,1],vmin=-0.6, vmax = 0.6, cbar = False)

ax[2,0].set_title('(e) Namaqua',size=13)
corrMatrix = NAMAQUA['Shelf_break'][fields].groupby(NAMAQUA['Shelf_break'].index.month).transform(lambda x: x-x.mean()).iloc[ind_summer.astype(int)].corr()
PMatrix = calculate_pvalues(NAMAQUA['Shelf_break'][fields].groupby(NAMAQUA['Shelf_break'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[2,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[2,1].set_title('(f) Namaqua',size=13)
corrMatrix = NAMAQUA['Upwelling'][fields].groupby(NAMAQUA['Upwelling'].index.month).transform(lambda x: x-x.mean()).corr()
PMatrix = calculate_pvalues(NAMAQUA['Upwelling'][fields].groupby(NAMAQUA['Upwelling'].index.month).transform(lambda x: x-x.mean()).iloc[ind_summer.astype(int)].corr())
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[2,1],vmin=-0.6, vmax = 0.6, cbar = False)

plt.show()

#fig.savefig('CORR_Summer', format='png', dpi=600)
plt.close()

########################################
# Winter
#########################################

fig, ax = plt.subplots(3, 2,figsize=(15,12))

fig.suptitle('Shelf-break and Upwelling: Winter [JJA]', fontsize=16)

ax[0,0].set_title('(a) St Helena',size=13)
corrMatrix = ST_HELENA['Shelf_break'][fields].groupby(ST_HELENA['Shelf_break'].index.month).transform(lambda x: x-x.mean()).iloc[ind_winter.astype(int)].corr()
PMatrix = calculate_pvalues(ST_HELENA['Shelf_break'][fields].groupby(ST_HELENA['Shelf_break'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[0,1].set_title('(b) St Helena',size=13)
corrMatrix = ST_HELENA['Upwelling'][fields].groupby(ST_HELENA['Upwelling'].index.month).transform(lambda x: x-x.mean()).iloc[ind_winter.astype(int)].corr()
PMatrix = calculate_pvalues(ST_HELENA['Upwelling'][fields].groupby(ST_HELENA['Upwelling'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,1],vmin=-0.6, vmax = 0.6, cbar = False)

ax[1,0].set_title('(c) Mid-shelf',size=13)
corrMatrix = MID_SHELF['Shelf_break'][fields].groupby(MID_SHELF['Shelf_break'].index.month).transform(lambda x: x-x.mean()).iloc[ind_winter.astype(int)].corr()
PMatrix = calculate_pvalues(MID_SHELF['Shelf_break'][fields].groupby(MID_SHELF['Shelf_break'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[1,1].set_title('(d) Mid-shelf',size=13)
corrMatrix = MID_SHELF['Upwelling'][fields].groupby(MID_SHELF['Upwelling'].index.month).transform(lambda x: x-x.mean()).iloc[ind_winter.astype(int)].corr()
PMatrix = calculate_pvalues(MID_SHELF['Upwelling'][fields].groupby(MID_SHELF['Upwelling'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,1],vmin=-0.6, vmax = 0.6, cbar = False)

ax[2,0].set_title('(e) Namaqua',size=13)
corrMatrix = NAMAQUA['Shelf_break'][fields].groupby(NAMAQUA['Shelf_break'].index.month).transform(lambda x: x-x.mean()).iloc[ind_winter.astype(int)].corr()
PMatrix = calculate_pvalues(NAMAQUA['Shelf_break'][fields].groupby(NAMAQUA['Shelf_break'].index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[2,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[2,1].set_title('(f) Namaqua',size=13)
corrMatrix = NAMAQUA['Upwelling'][fields].groupby(NAMAQUA['Upwelling'].index.month).transform(lambda x: x-x.mean()).corr()
PMatrix = calculate_pvalues(NAMAQUA['Upwelling'][fields].groupby(NAMAQUA['Upwelling'].index.month).transform(lambda x: x-x.mean()).iloc[ind_winter.astype(int)].corr())
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[2,1],vmin=-0.6, vmax = 0.6, cbar = False)

plt.show()

#fig.savefig('CORR_Winter', format='png', dpi=600)
plt.close()
##############################################################
#%%  
#################################################
# TEST ectreme events
#################################################

Threshold = 0.057 # 80 th percentile
fields = ['WIND','W','DIS','LEN','GRAD','NUM']

##################
# Choose region
##################

region = MID_SHELF['Upwelling'][fields].iloc[ind_summer.astype(int),:]
#########################################################

wind_region = region['WIND'].to_numpy() # Covert wind to np array 
wind_region = wind_region>=Threshold    # Boolean array

myarray = wind_region.astype(float)    # Binary array for the wind events

###############################################################################
def search_sequence_numpy(arr,seq):
    """ Find sequence in an array using NumPy only.

    Parameters
    ----------    
    arr    : input 1D array
    seq    : input 1D array

    Output
    ------    
    Output : 1D Array of indices in the input array that satisfy the 
    matching of input sequence in the input array.
    In case of no match, an empty list is returned.
    """

    # Store sizes of input array and sequence
    Na, Nseq = arr.size, seq.size

    # Range of sequence
    r_seq = np.arange(Nseq)

    # Create a 2D array of sliding indices across the entire length of input array.
    # Match up with the input sequence & get the matching starting indices.
    M = (arr[np.arange(Na-Nseq+1)[:,None] + r_seq] == seq).all(1)

    # Get the range of those indices as final output
    if M.any() >0:
        return np.where(np.convolve(M,np.ones((Nseq),dtype=int))>0)[0]
    else:
        return []         # No match found
###############################################################################

# Define a type of wind event and identify it
seq = np.array([0,0,1,1,1,0,0])  # Two days either side quiescent and then strong between

ind_event = search_sequence_numpy(myarray,seq)  # Search for the events
myevents = myarray[ind_event.astype(int)]    # Index array to get events

# Some basic info

num_total = (myarray == 1).sum()   # Number of extreme events
num_events = sum(myevents)/sum(seq) # Number of events matching criteria

###############################################################################

# Descritize the events into individual. Taking into the account that the convolution
# only flags the signal, sometimes, a value may be repetaed. Hence, we need to account for the consecutive
# signal and add the appropriate indice as a double to account for the two events

def my_sort(ind_event,num_events,seq):
    for i in range(1,num_events.astype(int)+1):
        tmp = ind_event[(i*len(seq))-len(seq):i*len(seq)]
        diff = tmp - tmp[0]
        if len(diff) == len(seq):
            ind_tmp = np.argwhere(diff>len(seq)-1)
            if ind_tmp.any() > 0:
                count_ind = ind_tmp.size
                Toadd = []
                for j in range(1,count_ind+1):
                    Toadd.append(tmp[0]-j)
                Toadd = sorted(Toadd)
                count = 0
                for k in range(1,len(Toadd)+1): 
                    ind_event = np.insert(ind_event, (i*len(seq))-len(seq)-k+count, Toadd[k-1])
                    count = count + 1
        else:
           Toadd = tmp[0] - np.arange(1,(len(seq) - len(tmp))+1,1)
           Toadd = sorted(Toadd)
           count = 0
           for k in range(1,len(Toadd)+1): 
               ind_event = np.insert(ind_event, (i*len(seq))-len(seq)-k+count, Toadd[k-1])
               count = count + 1
        # Check
#        tmp = ind_event[(i*len(seq))-len(seq):i*len(seq)]
#        check = myarray[tmp.astype(int)]
#        if sum(check) > sum(seq):
#            Toadd = []
#            for j in range(1,2+1):  # Padding with zeros
#               Toadd.append(tmp[0]-j)
#            Toadd = sorted(Toadd)
#            count = 0
#            for k in range(1,len(Toadd)+1): 
#                ind_event = np.insert(ind_event, (i*len(seq))-len(seq)-k+count, Toadd[k-1])
#                count = count + 1
    return ind_event 

ind_event = my_sort(ind_event,num_events,seq)
    

# Old version
#for i in range(1,num_events.astype(int)+1):
#    tmp = ind_event[(i*len(seq))-len(seq):i*len(seq)]
#    diff = tmp[1:]-tmp[0:-1]
#    ind_tmp = np.argwhere(diff>1)
#    if ind_tmp.any() > 0:
#        count_ind = ind_tmp.size
#        Toadd = []
#        for j in range(1,count_ind+1):
#            Toadd.append(tmp[0]-j)
#        for k in range(1,len(Toadd)+1):    
#            ind_event = np.insert(ind_event, (i*len(seq))-len(seq)-k, Toadd[k-1])
            
###############################################################################

# Quick look at the correaltions
gr1 = region.iloc[ind_event.astype(int),:]
fig, ax = plt.subplots(figsize=(15,12))
corrMatrix = gr1.corr()
PMatrix = calculate_pvalues(gr1)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax,vmin=-0.6, vmax = 0.6, cbar = False)

################################
# CREATE SOME PLOTS
################################

def STACK_array(array,num_events,seq):
    STACK = np.empty((num_events.astype(int),len(seq)))
    for i in range(1,num_events.astype(int)+1):
        grep = array[(i*len(seq))-len(seq):i*len(seq)]
        STACK[i-1,:] = grep
    return STACK

# Plot the wind-profiles

wind_region = region['WIND'].to_numpy()
wind_plot = wind_region[STACK_array(ind_event,num_events,seq).astype(int)]

fig, ax = plt.subplots(figsize=(15,12))
x = list(range(1, len(seq)+1))
for i in range(0,num_events.astype(int)):
    ax.plot(x,wind_plot[i,:],'-', color='black')
plt.plot(x,np.mean(wind_plot,axis=0),'--',color='red')
plt.xlabel('Days',fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel("Wind-stress [N/$m^{2}$]",fontsize=18)
plt.grid()
plt.show()
plt.close

# Plot another variable

# Displacement
myvar = region['GRAD'].to_numpy()
myvar_plot = myvar[STACK_array(ind_event,num_events,seq).astype(int)]

fig, ax = plt.subplots(figsize=(15,12))
x = list(range(1, len(seq)+1))
for i in range(0,num_events.astype(int)):
    tmp = myvar_plot[i,:]
    tmp_plot = tmp - tmp[0]
    myvar_plot[i,:] = tmp_plot
var_mean = np.nanmean(myvar_plot,axis=0)   
var_std = np.nanstd(myvar_plot,axis=0)
plt.plot(x,var_mean,'-',color='blue')
plt.fill_between(x, var_mean-var_std, var_mean+var_std,
    alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF',
    linewidth=4, linestyle='dashdot', antialiased=True)
plt.xlabel('Days',fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel("Variable",fontsize=18)
plt.grid()
plt.show()
plt.close

#%%
##############################################
# Little function to compute data for plots
##############################################

# When using this function: input must be a regional data frame, indexed correctly
# Output will be in order:
#    (1) The indexed data frame to the extreme events
#    (2) The wind_events that can be plotted
#    (3) The variables mean to be plotted (note that wind will be the first column)
#    (4) The variables std to be plotted (note that wind will be the first column)

def get_my_data(region,seq,Threshold):
    wind_region = region['WIND'].to_numpy() # Covert wind to np array 
    wind_region = wind_region>=Threshold    # Boolean array
    myarray = wind_region.astype(float)    # Binary array for the wind events
    ind_event = search_sequence_numpy(myarray,seq)  # Search for the events
    myevents = myarray[ind_event.astype(int)]    # Index array to get events
    num_total = (myarray == 1).sum()   # Number of extreme events
    num_events = sum(myevents)/sum(seq) # Number of events matching criteria
    print(num_total)
    print(num_events)
    gr1 = region.iloc[ind_event.astype(int),:]  # Get the data for the correlation maps
    ind_event = my_sort(ind_event,num_events,seq)
    wind_region = region['WIND'].to_numpy()
    wind_plot = wind_region[STACK_array(ind_event,num_events,seq).astype(int)]
    VAR_MEAN = np.zeros((len(seq),len(fields)))
    VAR_STD = np.zeros((len(seq),len(fields)))
    for k in range(0,len(fields)):
        myvar = region[fields[k]].to_numpy() 
        myvar_plot = myvar[STACK_array(ind_event,num_events,seq).astype(int)]
        for j in range(0,num_events.astype(int)):
            tmp = myvar_plot[j,:]
            tmp_plot = tmp - tmp[0]
            myvar_plot[j,:] = tmp_plot
        var_mean = np.nanmean(myvar_plot,axis=0)   
        var_std = np.nanstd(myvar_plot,axis=0)
        VAR_MEAN[:,k] = var_mean
        VAR_STD[:,k] = var_std
    return gr1, wind_plot, VAR_MEAN, VAR_STD  

w, x, y, z = get_my_data(region,seq,Threshold)      
#%%
############################
# MAIN PLOTS
############################

# We need to start with the correaltion plots to see if there is a significant 
# relationship between the length of the fronts and wind

#########################################################
# Correlation plot
##########################################################

seq = np.array([0,0,1,1,1,0,0])

#############
# Summer
#############

Threshold = 0.057

fig, ax = plt.subplots(3, 2,figsize=(15,12))

fig.suptitle('Shelf-break and Upwelling: Summer [DJF]', fontsize=16)

ax[0,0].set_title('(a) St Helena',size=13)
region = ST_HELENA['Shelf_break'][fields].iloc[ind_summer.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[0,1].set_title('(b) St Helena',size=13)
region = ST_HELENA['Upwelling'][fields].iloc[ind_summer.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,1],vmin=-0.6, vmax = 0.6, cbar = False)

ax[1,0].set_title('(c) Mid-shelf',size=13)
region = MID_SHELF['Shelf_break'][fields].iloc[ind_summer.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[1,1].set_title('(d) Mid-shelf',size=13)
region = MID_SHELF['Upwelling'][fields].iloc[ind_summer.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,1],vmin=-0.6, vmax = 0.6, cbar = False)

ax[2,0].set_title('(e) Namaqua',size=13)
region = NAMAQUA['Shelf_break'][fields].iloc[ind_summer.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[2,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[2,1].set_title('(f) Namaqua',size=13)
region = NAMAQUA['Upwelling'][fields].iloc[ind_summer.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[2,1],vmin=-0.6, vmax = 0.6, cbar = False)

plt.show()

fig.savefig('CORR_Summer_event_3day', format='png', dpi=600)
plt.close()

####################
# Winter
####################

Threshold = 0.022

fig, ax = plt.subplots(3, 2,figsize=(15,12))

fig.suptitle('Shelf-break and Upwelling: Winter [JJA]', fontsize=16)

ax[0,0].set_title('(a) St Helena',size=13)
region = ST_HELENA['Shelf_break'][fields].iloc[ind_winter.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[0,1].set_title('(b) St Helena',size=13)
region = ST_HELENA['Upwelling'][fields].iloc[ind_winter.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,1],vmin=-0.6, vmax = 0.6, cbar = False)

ax[1,0].set_title('(c) Mid-shelf',size=13)
region = MID_SHELF['Shelf_break'][fields].iloc[ind_winter.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[1,1].set_title('(d) Mid-shelf',size=13)
region = MID_SHELF['Upwelling'][fields].iloc[ind_winter.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,1],vmin=-0.6, vmax = 0.6, cbar = False)

ax[2,0].set_title('(e) Namaqua',size=13)
region = NAMAQUA['Shelf_break'][fields].iloc[ind_winter.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[2,0],vmin=-0.6, vmax = 0.6, cbar = False)

ax[2,1].set_title('(f) Namaqua',size=13)
region = NAMAQUA['Upwelling'][fields].iloc[ind_winter.astype(int)]
out = get_my_data(region,seq,Threshold)[0]
corrMatrix = out.corr()
PMatrix = calculate_pvalues(out)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[2,1],vmin=-0.6, vmax = 0.6, cbar = False)

plt.show()

fig.savefig('CORR_Winter_event_3day', format='png', dpi=600)
plt.close()


#%%

###############################
# Create line graphs
###############################

from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection 

NAME = ["ST_HELENA","MID_SHELF","NAMAQUA"]
MARKER = ["o","v","s"]
COLOR = ["deepskyblue","limegreen","coral"]

fields = ['WIND','DIS','GRAD','LEN','NUM','W']

region1 = ST_HELENA['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_summer.astype(int),:]

seq = np.array([0,0,1,1,1,0,0])  # Two days either side quiescent and then strong between

allregions_summer = [region1,region2,region3]

fig, ax = plt.subplots(figsize=(9,8))
ax.set_facecolor('lightgrey')

trans1 = Affine2D().translate(-0.2, 0.0) + ax.transData
trans2 = Affine2D().translate(0.0, 0.0) + ax.transData
trans3 = Affine2D().translate(0.2, 0.0) + ax.transData

TRANS = [trans1,trans2,trans3]

for i in range(0,len(allregions_summer)):
    inxval = list(range(1, len(seq)+1))
    out_summer = get_my_data(allregions_summer[i],seq,0.057)[2:4]
    myvar = out_summer[0][:,4]
    mystd = out_summer[1][:,4]
    plt.errorbar(inxval, myvar, marker=MARKER[i], 
                 markersize=14,
                 color= COLOR[i],
                 label=NAME[i],
                 yerr = mystd,
                 fmt="",
                 ecolor = COLOR[i],
                 elinewidth=3,
                 transform = TRANS[i])
    
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel("Distnace offshore",fontsize=18)
plt.legend(borderpad=1)
plt.grid()
plt.show()

#fig.savefig('WIND_STRESS_clim', format='png', dpi=600)
plt.close()

#%%

##################
# MAIN Figure
##################

# Function to quickly plot the variables 

def myplot(data_set,Threshold,title,var_index):
    inxval = list(range(1, len(seq)+1))
    ax.set_title(title,loc='left',size=13)
    ax.set_facecolor('lightgrey')
    trans1 = Affine2D().translate(-0.2, 0.0) + ax.transData
    trans2 = Affine2D().translate(0.0, 0.0) + ax.transData
    trans3 = Affine2D().translate(0.2, 0.0) + ax.transData
    TRANS = [trans1,trans2,trans3]
    ax.hlines(y=0, xmin=inxval[0], xmax=inxval[-1], linewidth=1, linestyle='--', color='b')
    for i in range(0,len(allregions_summer)):
        out_summer = get_my_data(data_set[i],seq,Threshold)[2:4]
        myvar = out_summer[0][:,var_index]
        mystd = out_summer[1][:,var_index]
        ax.errorbar(inxval, myvar, marker=MARKER[i], 
                     markersize=14,
                     color= COLOR[i],
                     label=NAME[i],
                     yerr = mystd,
                     fmt="",
                     ecolor = COLOR[i],
                     elinewidth=3,
                     transform = TRANS[i])
    return

# To save time, will do it manually as a lot of different info needs to be added
# For each region: will show the wind profile, distance, gradient and number of fronts
# The panel will be divided into two: left will be summer and right
# will be the winter

# Upwelling region

from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection 

###########################################################
NAME = ["ST_HELENA","MID_SHELF","NAMAQUA"]
MARKER = ["o","v","s"]
COLOR = ["deepskyblue","limegreen","coral"]
fields = ['WIND','DIS','GRAD']

#############################################################
region1 = ST_HELENA['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_summer.astype(int),:]

allregions_summer = [region1,region2,region3]

region1 = ST_HELENA['Upwelling'][fields].iloc[ind_winter.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_winter.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_winter.astype(int),:]

allregions_winter = [region1,region2,region3]
################################################################

Threshold_sum = 0.057
Threshold_win = 0.022
seq = np.array([0,0,1,1,1,0,0])  # Two days either side quiescent and then strong between

###################################################################
fig, axes = plt.subplots(3, 2,sharex=True,figsize=(15,12))

ax = axes
ax = ax[0,0]
ax.set_title('(a)       Summer [DJF]',size=13)
ax.set_facecolor('lightgrey')

trans1 = Affine2D().translate(-0.2, 0.0) + ax.transData
trans2 = Affine2D().translate(0.0, 0.0) + ax.transData
trans3 = Affine2D().translate(0.2, 0.0) + ax.transData

TRANS = [trans1,trans2,trans3]

for i in range(0,len(allregions_summer)):
    inxval = list(range(1, len(seq)+1))
    out_summer = get_my_data(allregions_summer[i],seq,Threshold_sum)[1]
    myvar = np.nanmean(out_summer,axis=0)
    mystd = np.nanstd(out_summer,axis=0)
    ax.errorbar(inxval, myvar, marker=MARKER[i], 
                 markersize=14,
                 color= COLOR[i],
                 label=NAME[i],
                 yerr = mystd,
                 fmt="",
                 ecolor = COLOR[i],
                 elinewidth=3,
                 transform = TRANS[i])

ax.set_ylabel("Wind-stress [N/$m^{2}$]",fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.legend(borderpad=1)
ax.grid()
del ax

ax = axes
ax = ax[0,1]
ax.set_title('(b)       Winter [JJA]',size=13)
ax.set_facecolor('lightgrey')

trans1 = Affine2D().translate(-0.2, 0.0) + ax.transData
trans2 = Affine2D().translate(0.0, 0.0) + ax.transData
trans3 = Affine2D().translate(0.2, 0.0) + ax.transData

TRANS = [trans1,trans2,trans3]

for i in range(0,len(allregions_summer)):
    inxval = list(range(1, len(seq)+1))
    out_summer = get_my_data(allregions_winter[i],seq,Threshold_win)[1]
    myvar = np.nanmean(out_summer,axis=0)
    mystd = np.nanstd(out_summer,axis=0)
    ax.errorbar(inxval, myvar, marker=MARKER[i], 
                 markersize=14,
                 color= COLOR[i],
                 label=NAME[i],
                 yerr = mystd,
                 fmt="",
                 ecolor = COLOR[i],
                 elinewidth=3,
                 transform = TRANS[i])

#ax.set_ylabel("Wind-stress [N/$m^{2}$]",fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid()
del ax

ax = axes
ax = ax[1,0]
myplot(allregions_summer,Threshold_sum,'(c)',1)
ax.set_ylabel('Distance [km]',fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid()
del ax

ax = axes
ax = ax[1,1]
myplot(allregions_winter,Threshold_win,'(d)',1)
#ax.set_ylabel('Distance [km]',fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid()
del ax

ax = axes
ax = ax[2,0]
myplot(allregions_summer,Threshold_sum,'(e)',2)
ax.set_xlabel('Days',fontsize=15)
ax.set_ylabel("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid()
del ax

ax = axes
ax = ax[2,1]
myplot(allregions_winter,Threshold_win,'(e)',2)
ax.set_xlabel('Days',fontsize=15)
#ax.set_ylabel('Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]',fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid()
del ax

fig.savefig('Upwelling_3day.png', format='png', dpi=600)

# Shelf-break

from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection 

###########################################################
NAME = ["ST_HELENA","MID_SHELF","NAMAQUA"]
MARKER = ["o","v","s"]
COLOR = ["deepskyblue","limegreen","coral"]
fields = ['WIND','DIS','GRAD']

#############################################################
region1 = ST_HELENA['Shelf_break'][fields].iloc[ind_summer.astype(int),:]
region2 = MID_SHELF['Shelf_break'][fields].iloc[ind_summer.astype(int),:]
region3 = NAMAQUA['Shelf_break'][fields].iloc[ind_summer.astype(int),:]

allregions_summer = [region1,region2,region3]

region1 = ST_HELENA['Shelf_break'][fields].iloc[ind_winter.astype(int),:]
region2 = MID_SHELF['Shelf_break'][fields].iloc[ind_winter.astype(int),:]
region3 = NAMAQUA['Shelf_break'][fields].iloc[ind_winter.astype(int),:]

allregions_winter = [region1,region2,region3]
################################################################

Threshold_sum = 0.057
Threshold_win = 0.022
seq = np.array([0,0,1,1,1,0,0])  # Two days either side quiescent and then strong between

###################################################################
fig, axes = plt.subplots(3, 2,sharex=True,figsize=(15,12))

ax = axes
ax = ax[0,0]
ax.set_title('(a)       Summer [DJF]',size=13)
ax.set_facecolor('lightgrey')

trans1 = Affine2D().translate(-0.2, 0.0) + ax.transData
trans2 = Affine2D().translate(0.0, 0.0) + ax.transData
trans3 = Affine2D().translate(0.2, 0.0) + ax.transData

TRANS = [trans1,trans2,trans3]

for i in range(0,len(allregions_summer)):
    inxval = list(range(1, len(seq)+1))
    out_summer = get_my_data(allregions_summer[i],seq,Threshold_sum)[1]
    myvar = np.nanmean(out_summer,axis=0)
    mystd = np.nanstd(out_summer,axis=0)
    ax.errorbar(inxval, myvar, marker=MARKER[i], 
                 markersize=14,
                 color= COLOR[i],
                 label=NAME[i],
                 yerr = mystd,
                 fmt="",
                 ecolor = COLOR[i],
                 elinewidth=3,
                 transform = TRANS[i])

ax.set_ylabel("Wind-stress [N/$m^{2}$]",fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.legend(borderpad=1)
ax.grid()
del ax

ax = axes
ax = ax[0,1]
ax.set_title('(b)       Winter [JJA]',size=13)
ax.set_facecolor('lightgrey')

trans1 = Affine2D().translate(-0.2, 0.0) + ax.transData
trans2 = Affine2D().translate(0.0, 0.0) + ax.transData
trans3 = Affine2D().translate(0.2, 0.0) + ax.transData

TRANS = [trans1,trans2,trans3]

for i in range(0,len(allregions_summer)):
    inxval = list(range(1, len(seq)+1))
    out_summer = get_my_data(allregions_winter[i],seq,Threshold_win)[1]
    myvar = np.nanmean(out_summer,axis=0)
    mystd = np.nanstd(out_summer,axis=0)
    ax.errorbar(inxval, myvar, marker=MARKER[i], 
                 markersize=14,
                 color= COLOR[i],
                 label=NAME[i],
                 yerr = mystd,
                 fmt="",
                 ecolor = COLOR[i],
                 elinewidth=3,
                 transform = TRANS[i])

#ax.set_ylabel("Wind-stress [N/$m^{2}$]",fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid()
del ax

ax = axes
ax = ax[1,0]
myplot(allregions_summer,Threshold_sum,'(c)',1)
ax.set_ylabel('Distance [km]',fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid()
del ax

ax = axes
ax = ax[1,1]
myplot(allregions_winter,Threshold_win,'(d)',1)
#ax.set_ylabel('Distance [km]',fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid()
del ax

ax = axes
ax = ax[2,0]
myplot(allregions_summer,Threshold_sum,'(e)',2)
ax.set_xlabel('Days',fontsize=15)
ax.set_ylabel("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid()
del ax

ax = axes
ax = ax[2,1]
myplot(allregions_winter,Threshold_win,'(f)',2)
ax.set_xlabel('Days',fontsize=15)
#ax.set_ylabel('Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]',fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid()
del ax

fig.savefig('Shelf_3day.png', format='png', dpi=600)

#%%

