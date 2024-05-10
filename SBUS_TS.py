#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 13:03:13 2023

@author: jono
"""

# Plot the average displacement and and vertical velcoity/gradeint comosite just for the 
# upwelligna dn shelf break region

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

#STEP 1:
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

#####################
# LOAD Data frame
#####################

f = open('SBUS_up.pckl', 'rb')
SBUS_up = pickle.load(f)
f.close()

f = open('SBUS_shelf.pckl', 'rb')
SBUS_shelf = pickle.load(f)
f.close()

df_shelf = SBUS_shelf
df_up = SBUS_up
#%%

############################################
# Plot combining displacmenet and gradient
############################################

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

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

fig, axes = plt.subplots(2,2,sharex=True, figsize=(18,8))

labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul',
          'Aug','Sep','Oct','Nov','Dec']

# Dsiplacement plot

# Shelf-break
dates = df_clim_shelf.index
y = np.array(df_clim_shelf.DIS)
c = np.array(df_clim_shelf.GRAD)
             
axes[0,0].set_title("(a) Shelf-break region [SBF]",loc='left', fontsize=16)
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
axes[0,0].set_ylim([80, 150])
axes[0,0].yaxis.set_major_locator(MultipleLocator(10.0))

# Upwelling
dates = df_clim_up.index
y = np.array(df_clim_up.DIS)
c = np.array(df_clim_up.GRAD)

axes[0,1].set_title("(b) Upwelling region [UF]",loc='left', fontsize=16)
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
axes[0,1].set_ylim([30, 90])
axes[0,1].yaxis.set_major_locator(MultipleLocator(10.0))
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a
plt.colorbar(line, ax=axes[0,:], label="Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]")

# Gradient plot

# Shelf-break
dates = df_clim_shelf.index
y = np.array(df_clim_shelf.GRAD)
c = np.array(df_clim_shelf.W)

axes[1,0].set_title("(c) Shelf-break region [SBF]",loc='left', fontsize=16)
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
axes[1,0].set_ylim([0.04, 0.18])
axes[1,0].yaxis.set_major_locator(MultipleLocator(0.02))
axes[1,0].set_ylabel("Gradient " +"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",fontsize=16)

# Upwelling
dates = df_clim_up.index
y = np.array(df_clim_up.GRAD)
c = np.array(df_clim_up.W)

axes[1,1].set_title("(d) Upwelling region [UF]",loc='left', fontsize=16)
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
axes[1,1].set_ylabel("Gradient " +"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",fontsize=16)
axes[1,1].xaxis.set_ticks(inxval)
axes[1,1].xaxis.set_ticklabels(labels,rotation=45,ha = 'right')
axes[1,1].set_ylim([0.04, 0.18])
axes[1,1].yaxis.set_major_locator(MultipleLocator(0.02))
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a

plt.colorbar(line, ax=axes[1, :], label='Vertical velocities [m/s]' )
plt.show()

fig.savefig('SBUS_TS', format='png', dpi=600)
plt.close()





