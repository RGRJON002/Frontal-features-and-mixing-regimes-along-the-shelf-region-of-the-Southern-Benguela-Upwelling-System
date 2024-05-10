#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 20:16:37 2022

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
# LOAD the V File
#######################

f = open('V.pckl', 'rb')
V_COAST = pickle.load(f)
f.close()

#######################
# LOAD the WIND File
#######################

f = open('WIND.pckl', 'rb')
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
V_SH = V_COAST[ind_CO_min:ind_CO_max,:] 
V_SH = np.nanmean(V_SH,0)


# Mid-shelf

CO_lat_min = -31.5   # Latitude of Cut min
CO_lat_max = -30.3   # Latitude of Cut max 
ind_CO_min = find_nearest(lat_CROCO,CO_lat_min)[1]
ind_CO_max = find_nearest(lat_CROCO,CO_lat_max)[1]
WIND_MS = WIND_COAST[ind_CO_min:ind_CO_max,:]
WIND_MS = np.nanmean(WIND_MS,0)
V_MS = V_COAST[ind_CO_min:ind_CO_max,:]
V_MS = np.nanmean(V_MS,0)


# Namaqua

CO_lat_min = -30.3   # Latitude of Cut min
CO_lat_max = -28   # Latitude of Cut max 
ind_CO_min = find_nearest(lat_CROCO,CO_lat_min)[1]
ind_CO_max = find_nearest(lat_CROCO,CO_lat_max)[1]
WIND_NM = WIND_COAST[ind_CO_min:ind_CO_max,:]
WIND_NM = np.nanmean(WIND_NM,0)
V_NM = V_COAST[ind_CO_min:ind_CO_max,:]
V_NM = np.nanmean(V_NM,0)

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
ST_HELENA['Shelf_break']['V'] = V_SH 
ST_HELENA['Upwelling']['V'] = V_SH

MID_SHELF['Shelf_break']['WIND'] = WIND_MS 
MID_SHELF['Upwelling']['WIND'] = WIND_MS
MID_SHELF['Shelf_break']['V'] = V_MS 
MID_SHELF['Upwelling']['V'] = V_MS

NAMAQUA['Shelf_break']['WIND'] = WIND_NM 
NAMAQUA['Upwelling']['WIND'] = WIND_NM
NAMAQUA['Shelf_break']['V'] = V_NM 
NAMAQUA['Upwelling']['V'] = V_NM

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

f = open('EDGE.pckl', 'rb')
EDGE = pickle.load(f)
f.close()

# BVF map

f = open('BVF.pckl', 'rb')
BVF = pickle.load(f)
f.close()

# Create a once off location of the offshore extent of the SST
BVF_up_SH = []
BVF_shelf_SH = []
BVF_up_MS = []
BVF_shelf_MS = []
BVF_up_NM = []
BVF_shelf_NM = []
BVF_shelf = []
BVF_up = []
for i in range(0,len(TIME)):
    loc_data = np.where(EDGE[:,:,i]*up_mask*SH_mask == 1)
    BVF_up_SH.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF
    del loc_data
    loc_data = np.where(EDGE[:,:,i]*sb_mask*SH_mask == 1)
    BVF_shelf_SH.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF
    del loc_data
    loc_data = np.where(EDGE[:,:,i]*up_mask*MS_mask == 1)
    BVF_up_MS.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF
    del loc_data
    loc_data = np.where(EDGE[:,:,i]*sb_mask*MS_mask == 1)
    BVF_shelf_MS.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF
    del loc_data
    loc_data = np.where(EDGE[:,:,i]*up_mask*NM_mask == 1)
    BVF_up_NM.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF
    del loc_data
    loc_data = np.where(EDGE[:,:,i]*sb_mask*NM_mask == 1)
    BVF_shelf_NM.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF
    del loc_data
    loc_data = np.where(EDGE[:,:,i]*up_mask == 1)
    BVF_up.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF
    del loc_data
    loc_data = np.where(EDGE[:,:,i]*sb_mask == 1)
    BVF_shelf.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF

## Clean up and calculate index

del BVF
del EDGE

## Add BVF

ST_HELENA['Shelf_break']['BVF'] = BVF_shelf_SH 
ST_HELENA['Upwelling']['BVF'] = BVF_up_SH

MID_SHELF['Shelf_break']['BVF'] = BVF_shelf_MS 
MID_SHELF['Upwelling']['BVF'] = BVF_up_MS

NAMAQUA['Shelf_break']['BVF'] = BVF_shelf_NM 
NAMAQUA['Upwelling']['BVF'] = BVF_up_NM

#%%

# To understand how fronts are affected, we apprecaite that we need to reduce the number of dimensions
# https://www.jcchouinard.com/pca-with-python/
# https://towardsdatascience.com/what-are-pca-loadings-and-biplots-9a7897f2e559
# https://support.minitab.com/en-us/minitab/21/help-and-how-to/statistical-modeling/multivariate/how-to/principal-components/interpret-the-results/all-statistics-and-graphs/
# Following the tutorial, going to do a PCA analysis for summer and winter to see if and how the drivers of the system
# https://nirpyresearch.com/detecting-outliers-using-mahalanobis-distance-pca-python/
# Change. We have a number of issues to tackle
# 1) Different seasons
# 2) Missing values
# 3) Distinguish groups (plot only every nth value to simplify the plot) 

# We must standardize the data as PCA scales the variance of the system relative to the variables
# Ideally, we want to also construct a loading plot and save that as well with a scree to see the contribution

# Liu, T., Krisch, S., Xie, R. C., Hopwood, M. J., Dengler, M., & Achterberg, E. P. (2022). 
# Sediment release in the Benguela Upwelling System dominates trace metal input to the shelf
# and eastern South Atlantic Ocean. Global Biogeochemical Cycles, 36, e2022GB007466. https://doi. org/10.1029/2022GB007466

########################################################## 
# What Diagnostics do we want to work with

fields = ['WIND','W','BVF','DIS','GRAD']

#fields = ['WIND','BVF','DIS','GRAD'] 

#feature_names = ["Wind-stress [N/$m^{2}$]""BVF [$s^{-2}$",
#                 "Distance [km]", "Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",
#                 ]

#feature_names = ["Wind-stress [N/$m^{2}$]","Vertical Velocity [m/s]","Along-shore current [m/s]","BVF [$s^{-2}$",
#                 "Distance [km]","Length [Pixels]", "Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",
#                 "Number of fronts"]

feature_names = fields

target_names = ["St Helena UF","St Helena SBF","Mid-Shelf UF", "Mid-shelf SBF", "Namaqua UF",
                "Namaqua SBF"]

########################
# Step 1 Seasons
########################

# Summer

ST_HELENA_summer_up = ST_HELENA['Upwelling'][fields].iloc[ind_summer.astype(int)]
ST_HELENA_summer_shelf = ST_HELENA['Shelf_break'][fields].iloc[ind_summer.astype(int)]

MID_SHELF_summer_up = MID_SHELF['Upwelling'][fields].iloc[ind_summer.astype(int)]
MID_SHELF_summer_shelf = MID_SHELF['Shelf_break'][fields].iloc[ind_summer.astype(int)]

NAMAQUA_summer_up = NAMAQUA['Upwelling'][fields].iloc[ind_summer.astype(int)]
NAMAQUA_summer_shelf = NAMAQUA['Shelf_break'][fields].iloc[ind_summer.astype(int)]

# Winter

ST_HELENA_winter_up = ST_HELENA['Upwelling'][fields].iloc[ind_winter.astype(int)]
ST_HELENA_winter_shelf = ST_HELENA['Shelf_break'][fields].iloc[ind_winter.astype(int)]

MID_SHELF_winter_up = MID_SHELF['Upwelling'][fields].iloc[ind_winter.astype(int)]
MID_SHELF_winter_shelf = MID_SHELF['Shelf_break'][fields].iloc[ind_winter.astype(int)]

NAMAQUA_winter_up = NAMAQUA['Upwelling'][fields].iloc[ind_winter.astype(int)]
NAMAQUA_winter_shelf = NAMAQUA['Shelf_break'][fields].iloc[ind_winter.astype(int)]


#########################
# Step 2 Missing values
#########################

# To ensure that the summer and winter arrays are all the same size, we need to 
# identify all the missing rows from every region for both fronts and remove all of it
# so as not to give any eroor

def my_ind(df):
    arr = df.to_numpy()
    arr_na = np.isnan(arr)
    arr_na = arr_na.astype(float) 
    arr_na = np.sum(arr_na,axis=1)
    ind = np.argwhere(arr_na >= 1)
    return ind

# Summer

SH_summer_up_ind = my_ind(ST_HELENA_summer_up)
SH_summer_shelf_ind = my_ind(ST_HELENA_summer_shelf)

MS_summer_up_ind = my_ind(MID_SHELF_summer_up)
MS_summer_shelf_ind = my_ind(MID_SHELF_summer_shelf)

NM_summer_up_ind = my_ind(NAMAQUA_summer_up)
NM_summer_shelf_ind = my_ind(NAMAQUA_summer_shelf)

ind = np.concatenate((SH_summer_up_ind, SH_summer_shelf_ind, MS_summer_up_ind,
                      MS_summer_shelf_ind,NM_summer_up_ind,NM_summer_shelf_ind))

uniqueValues, indicesList = np.unique(ind, return_index=True)

ST_HELENA_summer_up = ST_HELENA_summer_up.drop(ST_HELENA_summer_up.index[uniqueValues.astype(int)])
ST_HELENA_summer_shelf = ST_HELENA_summer_shelf.drop(ST_HELENA_summer_shelf.index[uniqueValues.astype(int)])

MID_SHELF_summer_up = MID_SHELF_summer_up.drop(MID_SHELF_summer_up.index[uniqueValues.astype(int)])
MID_SHELF_summer_shelf = MID_SHELF_summer_shelf.drop(MID_SHELF_summer_shelf.index[uniqueValues.astype(int)])

NAMAQUA_summer_up = NAMAQUA_summer_up.drop(NAMAQUA_summer_up.index[uniqueValues.astype(int)])
NAMAQUA_summer_shelf = NAMAQUA_summer_shelf.drop(NAMAQUA_summer_shelf.index[uniqueValues.astype(int)])

# Winter

SH_winter_up_ind = my_ind(ST_HELENA_winter_up)
SH_winter_shelf_ind = my_ind(ST_HELENA_winter_shelf)

MS_winter_up_ind = my_ind(MID_SHELF_winter_up)
MS_winter_shelf_ind = my_ind(MID_SHELF_winter_shelf)

NM_winter_up_ind = my_ind(NAMAQUA_winter_up)
NM_winter_shelf_ind = my_ind(NAMAQUA_winter_shelf)

ind = np.concatenate((SH_winter_up_ind, SH_winter_shelf_ind, MS_winter_up_ind,
                      MS_winter_shelf_ind,NM_winter_up_ind,NM_winter_shelf_ind))

uniqueValues, indicesList = np.unique(ind, return_index=True)

ST_HELENA_winter_up = ST_HELENA_winter_up.drop(ST_HELENA_winter_up.index[uniqueValues.astype(int)])
ST_HELENA_winter_shelf = ST_HELENA_winter_shelf.drop(ST_HELENA_winter_shelf.index[uniqueValues.astype(int)])

MID_SHELF_winter_up = MID_SHELF_winter_up.drop(MID_SHELF_winter_up.index[uniqueValues.astype(int)])
MID_SHELF_winter_shelf = MID_SHELF_winter_shelf.drop(MID_SHELF_winter_shelf.index[uniqueValues.astype(int)])

NAMAQUA_winter_up = NAMAQUA_winter_up.drop(NAMAQUA_winter_up.index[uniqueValues.astype(int)])
NAMAQUA_winter_shelf = NAMAQUA_winter_shelf.drop(NAMAQUA_winter_shelf.index[uniqueValues.astype(int)])

###################################
# Step 3: Need to create targets
###################################

# Summer

mylen = len(ST_HELENA_summer_up)

target_summer = np.concatenate((np.full((mylen,1),0),np.full((mylen,1),1),np.full((mylen,1),2),
                                np.full((mylen,1),3),np.full((mylen,1),4),np.full((mylen,1),5)))

# Winter

mylen = len(ST_HELENA_winter_up)

target_winter = np.concatenate((np.full((mylen,1),0),np.full((mylen,1),1),np.full((mylen,1),2),
                                np.full((mylen,1),3),np.full((mylen,1),4),np.full((mylen,1),5)))


###################################
# Step 4: Create our arrays
###################################

# We now want to convert our dataframes to numpy arrays and stack them

# Summer

ST_HELENA_summer_up = ST_HELENA_summer_up.to_numpy()
ST_HELENA_summer_shelf = ST_HELENA_summer_shelf.to_numpy()

MID_SHELF_summer_up = MID_SHELF_summer_up.to_numpy()
MID_SHELF_summer_shelf = MID_SHELF_summer_shelf.to_numpy()

NAMAQUA_summer_up = NAMAQUA_summer_up.to_numpy()
NAMAQUA_summer_shelf = NAMAQUA_summer_shelf.to_numpy()

SUMMER_ARRAY = np.vstack((ST_HELENA_summer_up,ST_HELENA_summer_shelf,MID_SHELF_summer_up,
                         MID_SHELF_summer_shelf,NAMAQUA_summer_up,NAMAQUA_summer_shelf)) 

# Winter

ST_HELENA_winter_up = ST_HELENA_winter_up.to_numpy()
ST_HELENA_summer_shelf = ST_HELENA_winter_shelf.to_numpy()

MID_SHELF_winter_up = MID_SHELF_winter_up.to_numpy()
MID_SHELF_winter_shelf = MID_SHELF_winter_shelf.to_numpy()

NAMAQUA_winter_up = NAMAQUA_winter_up.to_numpy()
NAMAQUA_winter_shelf = NAMAQUA_winter_shelf.to_numpy()

WINTER_ARRAY = np.vstack((ST_HELENA_winter_up,ST_HELENA_winter_shelf,MID_SHELF_winter_up,
                         MID_SHELF_winter_shelf,NAMAQUA_winter_up,NAMAQUA_winter_shelf)) 

#%%
# Construct a plot for summer and winter on one figure
######################################
# PCA analysis
######################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt  
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D

target_names = {
    0:'St Helena UF',
    1:'St Helena SBF',
    2:'Mid-Shelf UF',
    3:'Mid-shelf SBF',
    4:'Namaqua UF',
    5:'Namaqua SBF'
}

mycolors = ["deepskyblue","deepskyblue","limegreen","limegreen","coral","coral"]
myshapes = ['^','o','^','o','^','o']
#################################
# PICK our field
#################################

X = SUMMER_ARRAY
Y = target_summer

# Percent to plot = 

# data scaling so we have the same STD and mean
x_scaled = StandardScaler().fit_transform(X)
 
pca = PCA(n_components=3)
 
# Fit and transform data
pca_features = pca.fit_transform(x_scaled)

# Explained variance
variance = pca.explained_variance_
variance = variance/len(fields)*100
variance = np.round(variance,2)
 
# Create dataframe
pca_df = pd.DataFrame(data=pca_features, columns=['PC1', 'PC2', 'PC3'])
 
# Apply the targett names
pca_df['target'] = target_summer
pca_df['target'] = pca_df['target'].map(target_names)
 
# Create the scaled PCA dataframe
pca_df_scaled = pca_df.copy()
 
scaler_df = pca_df[['PC1', 'PC2', 'PC3']]
scaler = 1 / (scaler_df.max() - scaler_df.min())
 
for index in scaler.index:
    pca_df_scaled[index] *= scaler[index]
 
# Initialize the 3D graph
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(121, projection='3d')
 
# Define scaled features as arrays
xdata = pca_df_scaled['PC1']
ydata = pca_df_scaled['PC2']
zdata = pca_df_scaled['PC3']

# Plot 3D scatterplot of 
for i in range(0,len(target_names)):
    tmpx = xdata[int((i*(len(xdata)/len(target_names)))):int((i*(len(xdata)/len(target_names))+len(xdata)/len(target_names)))]
    tmpy = ydata[int((i*(len(ydata)/len(target_names)))):int((i*(len(ydata)/len(target_names))+len(ydata)/len(target_names)))]
    tmpz = zdata[int((i*(len(zdata)/len(target_names)))):int((i*(len(zdata)/len(target_names))+len(zdata)/len(target_names)))]
    ax.scatter3D(
        tmpx[::15], 
        tmpy[::15], 
        tmpz[::15], 
        color=mycolors[i],
        marker = myshapes[i],
        label = target_names[i],
        alpha=0.5)

# Define the x, y, z variables
loadings = pca.components_
xs = loadings[0]
ys = loadings[1]
zs = loadings[2]
 
# Plot the loadings
for i, varnames in enumerate(feature_names):
    ax.scatter(xs[i], ys[i], zs[i], s=200)
    ax.text(xs[i] + 0.01,ys[i] + 0.03,zs[i] + 0.05, varnames)

# Plot the arrows
x_arr = np.zeros(len(loadings[0]))
y_arr = z_arr = x_arr
ax.quiver(x_arr, y_arr, z_arr, xs, ys, zs)
 
# Plot title of graph
ax.set_title(f'3D Biplot of Frontal Diagnostics [DJF]')
 
# Plot x, y, z labels
ax.set_xlabel('PC1 ' + '(' + str(variance[0]) + '%' +')', rotation=150)
ax.set_ylabel('P2 ' + '(' + str(variance[1]) + '%' +')')
ax.set_zlabel('P3 ' + '(' + str(variance[2]) + '%' +')', rotation=60)

plt.legend(fontsize=9,loc='upper left')

ax = fig.add_subplot(122, projection='3d')

X = WINTER_ARRAY
Y = target_winter

# Percent to plot = 

# data scaling so we have the same STD and mean
x_scaled = StandardScaler().fit_transform(X)
 
pca = PCA(n_components=3)
 
# Fit and transform data
pca_features = pca.fit_transform(x_scaled)

# Explained variance
variance = pca.explained_variance_
variance = variance/len(fields)*100
variance = np.round(variance,2)
 
# Create dataframe
pca_df = pd.DataFrame(data=pca_features, columns=['PC1', 'PC2', 'PC3'])
 
# Apply the targett names
pca_df['target'] = target_winter
pca_df['target'] = pca_df['target'].map(target_names)
 
# Create the scaled PCA dataframe
pca_df_scaled = pca_df.copy()
 
scaler_df = pca_df[['PC1', 'PC2', 'PC3']]
scaler = 1 / (scaler_df.max() - scaler_df.min())
 
for index in scaler.index:
    pca_df_scaled[index] *= scaler[index]

# Define scaled features as arrays
xdata = pca_df_scaled['PC1']
ydata = pca_df_scaled['PC2']
zdata = pca_df_scaled['PC3']

# Plot 3D scatterplot of 
for i in range(0,len(target_names)):
    tmpx = xdata[int((i*(len(xdata)/len(target_names)))):int((i*(len(xdata)/len(target_names))+len(xdata)/len(target_names)))]
    tmpy = ydata[int((i*(len(ydata)/len(target_names)))):int((i*(len(ydata)/len(target_names))+len(ydata)/len(target_names)))]
    tmpz = zdata[int((i*(len(zdata)/len(target_names)))):int((i*(len(zdata)/len(target_names))+len(zdata)/len(target_names)))]
    ax.scatter3D(
        tmpx[::25], 
        tmpy[::25], 
        tmpz[::25], 
        color=mycolors[i],
        marker = myshapes[i],
        label = target_names[i],
        alpha=0.5)

# Define the x, y, z variables
loadings = pca.components_
xs = loadings[0]
ys = loadings[1]
zs = loadings[2]
 
# Plot the loadings
for i, varnames in enumerate(feature_names):
    ax.scatter(xs[i], ys[i], zs[i], s=200)
    ax.text(xs[i] + 0.01,ys[i] + 0.03,zs[i] + 0.05, varnames)

# Plot the arrows
x_arr = np.zeros(len(loadings[0]))
y_arr = z_arr = x_arr
ax.quiver(x_arr, y_arr, z_arr, xs, ys, zs)
 
# Plot title of graph
ax.set_title(f'3D Biplot of Frontal Diagnostics [JJA]')
 
# Plot x, y, z labels
ax.set_xlabel('PC1 ' + '(' + str(variance[0]) + '%' +')', rotation=150)
ax.set_ylabel('P2 ' + '(' + str(variance[1]) + '%' +')')
ax.set_zlabel('P3 ' + '(' + str(variance[2]) + '%' +')', rotation=60)

plt.show()

#fig.savefig('PCA', format='png', dpi=600)

#%%

# Loop over the time-domain 

from geopy import distance

f = open('EDGE.pckl', 'rb')
EDGE = pickle.load(f)
f.close()

f = open('W.pckl', 'rb')
W = pickle.load(f)
f.close()

f = open('GRADIENT.pckl', 'rb')
GRADIENT = pickle.load(f)
f.close()

#############################
# OUTPUTS
#############################

# With this algorithm, it is important to store as much information as
# possible. Consequnetly, the folloing will be recorded per time-step
# 1) Mean displacement in km from the coastline of the fronts
# 2) Mean length of fronts
# 3) Number of front detected
# 4) The vertical velocities
# 5) The pixel SST gradient
# 6) Orientation ( in prep)

# First get data for a given mask

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
for i in range(0,len(TIME)): 
    loc_data = np.where(EDGE[:,:,i]*sb_mask == 1)
    W_shelf.append(np.nanmean(W[loc_data[0],loc_data[1],i]))  # Vertical velocity
    GRAD_shelf.append(np.nanmean(GRADIENT.data[loc_data[0],loc_data[1],i]))  # Pixel gradient
    pcm = map.contour(X, Y, EDGE[:,:,i]*sb_mask, linewidths=0.1, linestyles='solid', colors='black')
    pcm = pcm.allsegs[0]
    NUM_shelf.append(len(pcm))  # Number of fronts detected 
    dis_mean = []
    len_mean = []
    for k in range(0,len(pcm)):
        cc = pcm[k]
        Xloc, Yloc = contour_loc(cc,x,y)
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
for i in range(0,len(TIME)): 
    loc_data = np.where(EDGE[:,:,i]*up_mask == 1)
    W_up.append(np.nanmean(W[loc_data[0],loc_data[1],i]))  # Vertical velocity
    GRAD_up.append(np.nanmean(GRADIENT.data[loc_data[0],loc_data[1],i]))  # Pixel gradient
    pcm = map.contour(X, Y, EDGE[:,:,i]*up_mask, linewidths=0.1, linestyles='solid', colors='black')
    pcm = pcm.allsegs[0]
    NUM_up.append(len(pcm))  # Number of fronts detected 
    dis_mean = []
    len_mean = []
    for k in range(0,len(pcm)):
        cc = pcm[k]
        Xloc, Yloc = contour_loc(cc,x,y)
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
    DIS_up.append(np.mean(dis_mean)) # Displacement of fronts
    LEN_up.append(np.mean(len_mean)) # Length of fronts
    
################################
# SHELF-BREAK data frame
################################    
    
import pandas as pd    

data = {'Date': DATE,
        'DIS': DIS_shelf,
        'LEN': LEN_shelf,
        'NUM': NUM_shelf,
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
        'GRAD': GRAD_up,
        'W': W_up}
df_up = pd.DataFrame(data) 
df_up = df_up.set_index('Date')   

################
# SAVE
###############

import pickle

f = open('SBUS_up.pckl', 'wb')
pickle.dump(df_up, f)
f.close()

f = open('SBUS_shelf.pckl', 'wb')
pickle.dump(df_shelf, f)
f.close()   

#%% Quick plot of the correlations for the shelf-break and upwelling regions for
# 1) the entire time domain
# 2) Summer
# 3 Winter

f = open('SBUS_up.pckl', 'rb')
SBUS_up = pickle.load(f)
f.close()

f = open('SBUS_shelf.pckl', 'rb')
SBUS_shelf = pickle.load(f)
f.close()

####################################
# Add wind and stratification
####################################

SBUS_up['WIND'] = WIND_up
SBUS_up['BVF'] = BVF_up

SBUS_shelf['WIND'] = WIND_shelf
SBUS_shelf['BVF'] = BVF_shelf

#########################################

fields = ['WIND','W','BVF','DIS','LEN','GRAD','NUM']

import seaborn as sn

from scipy.stats import pearsonr

def calculate_pvalues(df):
    df = df.dropna()._get_numeric_data()
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            pvalues[r][c] = round(pearsonr(df[r], df[c])[1], 5)
    return pvalues

##################################
# Summer
##################################

fig, ax = plt.subplots(2, 2,figsize=(7,7))
fig.tight_layout()
cbar_ax = fig.add_axes([.97, .38, .04, .3])

ax[0,0].set_title('(a) Shelf-break region [DJF]',size=13)
corrMatrix = SBUS_shelf[fields].groupby(SBUS_shelf.index.month).transform(lambda x: x-x.mean()).iloc[ind_summer.astype(int)].corr()
PMatrix = calculate_pvalues(SBUS_shelf[fields].groupby(SBUS_shelf.index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,0],vmin=-0.9, vmax = 0.9, linewidth=.5, cbar = True, cbar_ax=cbar_ax, cmap="RdYlGn")
ax[0,0].set_facecolor('lightgrey')

ax[0,1].set_title('(b) Upwelling region [DJF]',size=13)
corrMatrix = SBUS_up[fields].groupby(SBUS_up.index.month).transform(lambda x: x-x.mean()).iloc[ind_summer.astype(int)].corr()
PMatrix = calculate_pvalues(SBUS_up[fields].groupby(SBUS_up.index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[0,1],vmin=-0.9, vmax = 0.9, linewidth=.5, cbar = True, cbar_ax=cbar_ax, cmap="RdYlGn")
ax[0,1].set_facecolor('lightgrey')

ax[1,0].set_title('(c) Shelf-break region [JJA]',size=13)
corrMatrix = SBUS_shelf[fields].groupby(SBUS_shelf.index.month).transform(lambda x: x-x.mean()).iloc[ind_winter.astype(int)].corr()
PMatrix = calculate_pvalues(SBUS_shelf[fields].groupby(SBUS_shelf.index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,0],vmin=-0.9, vmax = 0.9, linewidth=.5, cbar = True, cbar_ax=cbar_ax, cmap="RdYlGn")
ax[1,0].set_facecolor('lightgrey')

ax[1,1].set_title('(d) Upwelling region [JJA]',size=13)
corrMatrix = SBUS_up[fields].groupby(SBUS_up.index.month).transform(lambda x: x-x.mean()).iloc[ind_winter.astype(int)].corr()
PMatrix = calculate_pvalues(SBUS_up[fields].groupby(SBUS_up.index.month).transform(lambda x: x-x.mean()))
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
mask[np.invert(np.tril(PMatrix<0.05))] = True
sn.heatmap(corrMatrix, mask=mask, annot=True, ax = ax[1,1],vmin=-0.9, vmax = 0.9, linewidth=.5, cbar = True, cbar_ax=cbar_ax, cmap="RdYlGn")
ax[1,1].set_facecolor('lightgrey')

plt.show()

fig.savefig('CORR_SBUS', format='png', bbox_inches='tight', dpi=600)
plt.close()

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

fig, ax = plt.subplots(1,2,figsize=(15,5))
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
plt.grid()

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
plt.tick_params(axis='both', which='major', labelsize=16)
ax[1].set_ylabel("log10 BVF [$s^{-2}$]",fontsize=18)
ax[1].set_title("(b)",fontsize=15,loc='left')
ax[1].legend(borderpad=1,loc='lower left')
plt.grid()

plt.show()

#fig.savefig('WIND_STRESS_and_BVF_clim', format='png', dpi=600)
plt.close()

#%%

# Conduct a t-test for the upwelling and shelf region to see if the average frontal length is different and in which season are fronts
# longer. Furthermore, doing the same test, is there a difference in frontal length.

import scipy.stats as stats

# Length

# Upwelling

X_len = SBUS_up['LEN'].iloc[ind_summer.astype(int)]
Y_len = SBUS_up['LEN'].iloc[ind_winter.astype(int)]

X_len = X_len.dropna()
Y_len = Y_len.dropna()

print(np.nanmean(X_len))
print(np.nanmean(Y_len))

# Perform the two sample t-test with equal variances
stats.ttest_ind(a=X_len, b=Y_len, equal_var=True)

# Shelf

X_len = SBUS_shelf['LEN'].iloc[ind_summer.astype(int)]
Y_len = SBUS_shelf['LEN'].iloc[ind_winter.astype(int)]

X_len = X_len.dropna()
Y_len = Y_len.dropna()

print(np.nanmean(X_len))
print(np.nanmean(Y_len))

# Perform the two sample t-test with equal variances
stats.ttest_ind(a=X_len, b=Y_len, equal_var=True)

###############################################################################

# Number

# Upwelling

X_len = SBUS_up['NUM'].iloc[ind_summer.astype(int)]
Y_len = SBUS_up['NUM'].iloc[ind_winter.astype(int)]

X_len = X_len.dropna()
Y_len = Y_len.dropna()

print(np.nanmean(X_len))
print(np.nanmean(Y_len))

# Perform the two sample t-test with equal variances
stats.ttest_ind(a=X_len, b=Y_len, equal_var=True)

# Shelf

X_len = SBUS_shelf['NUM'].iloc[ind_summer.astype(int)]
Y_len = SBUS_shelf['NUM'].iloc[ind_winter.astype(int)]

X_len = X_len.dropna()
Y_len = Y_len.dropna()

print(np.nanmean(X_len))
print(np.nanmean(Y_len))

# Perform the two sample t-test with equal variances
stats.ttest_ind(a=X_len, b=Y_len, equal_var=True)

###############################################################################

# Gradient

X_len = SBUS_up['GRAD'].iloc[ind_summer.astype(int)]
Y_len = SBUS_up['GRAD'].iloc[ind_winter.astype(int)]

X_len = X_len.dropna()
Y_len = Y_len.dropna()

print(np.nanmean(X_len))
print(np.nanmean(Y_len))

# Perform the two sample t-test with equal variances
stats.ttest_ind(a=X_len, b=Y_len, equal_var=True)

# Shelf

X_len = SBUS_shelf['GRAD'].iloc[ind_summer.astype(int)]
Y_len = SBUS_shelf['GRAD'].iloc[ind_winter.astype(int)]

X_len = X_len.dropna()
Y_len = Y_len.dropna()

print(np.nanmean(X_len))
print(np.nanmean(Y_len))

# Perform the two sample t-test with equal variances
stats.ttest_ind(a=X_len, b=Y_len, equal_var=True)
