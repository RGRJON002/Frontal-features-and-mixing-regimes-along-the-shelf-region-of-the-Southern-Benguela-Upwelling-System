#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 14:01:44 2022

@author: jono
"""

#### This script is intended to visualize the frontal length and numbers and 
#### as well as the orientation of the fronts

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

#######################
# LOAD in the data
#######################

f = open('ST_HELENA.pckl', 'rb')
ST_HELENA = pickle.load(f)
f.close()

f = open('MID_SHELF.pckl', 'rb')
MID_SHELF = pickle.load(f)
f.close()

f = open('NAMAQUA.pckl', 'rb')
NAMAQUA = pickle.load(f)
f.close()

#%% 

#########################################################################
# Visualize frontal number and lemgth in conjunction with orientation
##########################################################################

#################
# CLIM PLOT
#################

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection 

# Compute the monthly means for each variable

#data_columns = ['DIS', 'LEN', 'NUM', 'NUM_ACROSS','NUM_LONG','GRAD','W']
data_columns = ['DIS', 'LEN', 'NUM','GRAD','W']

fig, ax = plt.subplots(2,2,sharex=False, figsize=(18,12))

labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul',
          'Aug','Sep','Oct','Nov','Dec']

inxval = list(range(1, 13))
labels_shelf = ['St Helena','Mid Shelf','Namaqua']

x = ST_HELENA['Shelf_break'].groupby([ST_HELENA['Shelf_break'].index.month]).mean().index 

######################################
# FRONT NUMBERS
######################################

#####################
# Shelf-break region
#####################

ax[0,0].set_title("(a) Shelf-break region",loc='left', fontsize=16)

y1 = ST_HELENA['Shelf_break'].groupby([ST_HELENA['Shelf_break'].index.month]).mean().NUM 
z1 = MID_SHELF['Shelf_break'].groupby([MID_SHELF['Shelf_break'].index.month]).mean().NUM 
k1 = NAMAQUA['Shelf_break'].groupby([NAMAQUA['Shelf_break'].index.month]).mean().NUM 

#y1 = ST_HELENA['Shelf_break'].groupby([ST_HELENA['Shelf_break'].index.month]).mean().NUM_ACROSS 
#z1 = MID_SHELF['Shelf_break'].groupby([MID_SHELF['Shelf_break'].index.month]).mean().NUM_ACROSS 
#k1 = NAMAQUA['Shelf_break'].groupby([NAMAQUA['Shelf_break'].index.month]).mean().NUM_ACROSS 

#y2 = ST_HELENA['Shelf_break'].groupby([ST_HELENA['Shelf_break'].index.month]).mean().NUM_LONG 
#z2 = MID_SHELF['Shelf_break'].groupby([MID_SHELF['Shelf_break'].index.month]).mean().NUM_LONG 
#k2 = NAMAQUA['Shelf_break'].groupby([NAMAQUA['Shelf_break'].index.month]).mean().NUM_LONG

y_std = ST_HELENA['Shelf_break'].groupby([ST_HELENA['Shelf_break'].index.month]).std().NUM
z_std = MID_SHELF['Shelf_break'].groupby([MID_SHELF['Shelf_break'].index.month]).std().NUM 
k_std = NAMAQUA['Shelf_break'].groupby([NAMAQUA['Shelf_break'].index.month]).std().NUM 

ax[0,0].bar(x-0.2, y1, width=0.2, color='skyblue', align='center', hatch='')
ax[0,0].bar(x, z1, width=0.2, color='lawngreen', align='center', hatch='')
ax[0,0].bar(x+0.2, k1, width=0.2, color='lightcoral', align='center', hatch='')

#ax[0,0].bar(x-0.2, y2, width=0.2, bottom=y1, color='skyblue', hatch='///')
#ax[0,0].bar(x, z2, width=0.2, bottom=z1, color='lawngreen', hatch='///')
#ax[0,0].bar(x+0.2, k2, width=0.2, bottom=k1, color='lightcoral', hatch='///')

ax[0,0].errorbar(x-0.2, ST_HELENA['Shelf_break'].groupby([ST_HELENA['Shelf_break'].index.month]).mean().NUM,
               yerr=y_std, fmt="o", color="blue")

ax[0,0].errorbar(x, MID_SHELF['Shelf_break'].groupby([MID_SHELF['Shelf_break'].index.month]).mean().NUM,
               yerr=z_std, fmt="o", color="green")

ax[0,0].errorbar(x+0.2, NAMAQUA['Shelf_break'].groupby([NAMAQUA['Shelf_break'].index.month]).mean().NUM,
               yerr=k_std, fmt="o", color="red")

ax[0,0].grid()
ax[0,0].set_ylim([0, 12])
ax[0,0].xaxis.set_ticks(inxval)
ax[0,0].xaxis.set_ticklabels(labels,rotation=45,ha = 'right')
ax[0,0].tick_params(axis="x", labelsize=14)
ax[0,0].set_ylabel("Number of fronts",fontsize=16)
ax[0,0].tick_params(axis="y", labelsize=14)

########################
# Upwelling
########################

ax[0,1].set_title("(b) Upwelling region",loc='left', fontsize=16)

y1 = ST_HELENA['Upwelling'].groupby([ST_HELENA['Upwelling'].index.month]).mean().NUM 
z1 = MID_SHELF['Upwelling'].groupby([MID_SHELF['Upwelling'].index.month]).mean().NUM 
k1 = NAMAQUA['Upwelling'].groupby([NAMAQUA['Upwelling'].index.month]).mean().NUM 

#y1 = ST_HELENA['Upwelling'].groupby([ST_HELENA['Upwelling'].index.month]).mean().NUM_ACROSS 
#z1 = MID_SHELF['Upwelling'].groupby([MID_SHELF['Upwelling'].index.month]).mean().NUM_ACROSS 
#k1 = NAMAQUA['Upwelling'].groupby([NAMAQUA['Upwelling'].index.month]).mean().NUM_ACROSS 

#y2 = ST_HELENA['Upwelling'].groupby([ST_HELENA['Upwelling'].index.month]).mean().NUM_LONG 
#z2 = MID_SHELF['Upwelling'].groupby([MID_SHELF['Upwelling'].index.month]).mean().NUM_LONG 
#k2 = NAMAQUA['Upwelling'].groupby([NAMAQUA['Upwelling'].index.month]).mean().NUM_LONG

y_std = ST_HELENA['Upwelling'].groupby([ST_HELENA['Upwelling'].index.month]).std().NUM
z_std = MID_SHELF['Upwelling'].groupby([MID_SHELF['Upwelling'].index.month]).std().NUM 
k_std = NAMAQUA['Upwelling'].groupby([NAMAQUA['Upwelling'].index.month]).std().NUM 

ax[0,1].bar(x-0.2, y1, width=0.2, color='skyblue', align='center', hatch='')
ax[0,1].bar(x, z1, width=0.2, color='lawngreen', align='center', hatch='')
ax[0,1].bar(x+0.2, k1, width=0.2, color='lightcoral', align='center', hatch='')

#ax[0,1].bar(x-0.2, y2, width=0.2, bottom=y1, color='skyblue', hatch='///')
#ax[0,1].bar(x, z2, width=0.2, bottom=z1, color='lawngreen', hatch='///')
#ax[0,1].bar(x+0.2, k2, width=0.2, bottom=k1, color='lightcoral', hatch='///')

ax[0,1].errorbar(x-0.2, ST_HELENA['Upwelling'].groupby([ST_HELENA['Upwelling'].index.month]).mean().NUM,
               yerr=y_std, fmt="o", color="blue")

ax[0,1].errorbar(x, MID_SHELF['Upwelling'].groupby([MID_SHELF['Upwelling'].index.month]).mean().NUM,
               yerr=z_std, fmt="o", color="green")

ax[0,1].errorbar(x+0.2, NAMAQUA['Upwelling'].groupby([NAMAQUA['Upwelling'].index.month]).mean().NUM,
               yerr=k_std, fmt="o", color="red")

ax[0,1].grid()
ax[0,1].set_ylim([0, 18])
ax[0,1].xaxis.set_ticks(inxval)
ax[0,1].xaxis.set_ticklabels(labels,rotation=45,ha = 'right')
ax[0,1].tick_params(axis="x", labelsize=14)
ax[0,1].tick_params(axis="y", labelsize=14)

###########################################
# FRONT LENGTH
###########################################

#####################
# Shelf-break region
#####################

ax[1,0].set_title("(c) Shelf-break region",loc='left', fontsize=16)

y1 = ST_HELENA['Shelf_break'].groupby([ST_HELENA['Shelf_break'].index.month]).mean().LEN 
z1 = MID_SHELF['Shelf_break'].groupby([MID_SHELF['Shelf_break'].index.month]).mean().LEN 
k1 = NAMAQUA['Shelf_break'].groupby([NAMAQUA['Shelf_break'].index.month]).mean().LEN 

y_std = ST_HELENA['Shelf_break'].groupby([ST_HELENA['Shelf_break'].index.month]).std().LEN
z_std = MID_SHELF['Shelf_break'].groupby([MID_SHELF['Shelf_break'].index.month]).std().LEN 
k_std = NAMAQUA['Shelf_break'].groupby([NAMAQUA['Shelf_break'].index.month]).std().LEN 

ax[1,0].bar(x-0.2, y1, width=0.2, color='skyblue', align='center', hatch='')
ax[1,0].bar(x, z1, width=0.2, color='lawngreen', align='center', hatch='')
ax[1,0].bar(x+0.2, k1, width=0.2, color='lightcoral', align='center', hatch='')

ax[1,0].errorbar(x-0.2, ST_HELENA['Shelf_break'].groupby([ST_HELENA['Shelf_break'].index.month]).mean().LEN,
               yerr=y_std, fmt="o", color="blue")

ax[1,0].errorbar(x, MID_SHELF['Shelf_break'].groupby([MID_SHELF['Shelf_break'].index.month]).mean().LEN,
               yerr=z_std, fmt="o", color="green")

ax[1,0].errorbar(x+0.2, NAMAQUA['Shelf_break'].groupby([NAMAQUA['Shelf_break'].index.month]).mean().LEN,
               yerr=k_std, fmt="o", color="red")

ax[1,0].grid()
ax[1,0].set_ylim([0, 70])
ax[1,0].xaxis.set_ticks(inxval)
ax[1,0].xaxis.set_ticklabels(labels,rotation=45,ha = 'right')
ax[1,0].tick_params(axis="x", labelsize=14)
ax[1,0].set_ylabel("Frontal length [pixels]",fontsize=16)
ax[1,0].tick_params(axis="y", labelsize=14)

########################
# Upwelling
########################
ax[1,1].set_title("(d) Upwelling region",loc='left', fontsize=16)

y1 = ST_HELENA['Upwelling'].groupby([ST_HELENA['Upwelling'].index.month]).mean().LEN 
z1 = MID_SHELF['Upwelling'].groupby([MID_SHELF['Upwelling'].index.month]).mean().LEN 
k1 = NAMAQUA['Upwelling'].groupby([NAMAQUA['Upwelling'].index.month]).mean().LEN 

y_std = ST_HELENA['Upwelling'].groupby([ST_HELENA['Upwelling'].index.month]).std().LEN
z_std = MID_SHELF['Upwelling'].groupby([MID_SHELF['Upwelling'].index.month]).std().LEN 
k_std = NAMAQUA['Upwelling'].groupby([NAMAQUA['Upwelling'].index.month]).std().LEN 

ax[1,1].bar(x-0.2, y1, width=0.2, color='skyblue', align='center', hatch='')
ax[1,1].bar(x, z1, width=0.2, color='lawngreen', align='center', hatch='')
ax[1,1].bar(x+0.2, k1, width=0.2, color='lightcoral', align='center', hatch='')

ax[1,1].errorbar(x-0.2, ST_HELENA['Upwelling'].groupby([ST_HELENA['Upwelling'].index.month]).mean().LEN,
               yerr=y_std, fmt="o", color="blue")

ax[1,1].errorbar(x, MID_SHELF['Upwelling'].groupby([MID_SHELF['Upwelling'].index.month]).mean().LEN,
               yerr=z_std, fmt="o", color="green")

ax[1,1].errorbar(x+0.2, NAMAQUA['Upwelling'].groupby([NAMAQUA['Upwelling'].index.month]).mean().LEN,
               yerr=k_std, fmt="o", color="red")

ax[1,1].grid()
ax[1,1].set_ylim([0, 80])
ax[1,1].xaxis.set_ticks(inxval)
ax[1,1].xaxis.set_ticklabels(labels,rotation=45,ha = 'right')

#################################################################################

ax[0,0].legend(labels_shelf,loc='upper left',fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.show()

fig.savefig('Front_num_length_clim', format='png', dpi=600)
plt.close()
