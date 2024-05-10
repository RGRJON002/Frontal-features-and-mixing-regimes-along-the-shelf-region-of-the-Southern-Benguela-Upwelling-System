#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 19:23:12 2023

@author: jono
"""
# Create a PCA plot just showing the shelf and Upwelling region

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

# New mask for offshore region

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

#%%

#######################
# LOAD the WIND File
#######################

f = open('WIND.pckl', 'rb')
WIND_COAST = pickle.load(f)
f.close()

#%%

# Calculate values                
            
###########################################
# Just get region north of Cape Columbine
###########################################

WIND_up = np.nanmean(WIND_COAST,axis=0)
WIND_shelf = np.nanmean(WIND_COAST,axis=0)


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

#%%

####################################
# Add wind and stratification
####################################

SBUS_up['WIND'] = WIND_up
SBUS_up['BVF'] = BVF_up

SBUS_shelf['WIND'] = WIND_shelf
SBUS_shelf['BVF'] = BVF_shelf

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

feature_names = fields

target_names = ["SBUS UF","SBUS SBF"]

########################
# Step 1 Seasons
########################

# Summer

SBUS_summer_up = SBUS_up[fields].iloc[ind_summer.astype(int)]
SBUS_summer_shelf = SBUS_shelf[fields].iloc[ind_summer.astype(int)]

# Winter

SBUS_winter_up = SBUS_up[fields].iloc[ind_winter.astype(int)]
SBUS_winter_shelf = SBUS_shelf[fields].iloc[ind_winter.astype(int)]

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

SBUS_summer_up_ind = my_ind(SBUS_summer_up)
SBUS_summer_shelf_ind = my_ind(SBUS_summer_shelf)

ind = np.concatenate((SBUS_summer_up_ind, SBUS_summer_shelf_ind))
uniqueValues, indicesList = np.unique(ind, return_index=True)

SBUS_summer_up = SBUS_summer_up.drop(SBUS_summer_up.index[uniqueValues.astype(int)])
SBUS_summer_shelf = SBUS_summer_shelf.drop(SBUS_summer_shelf.index[uniqueValues.astype(int)])

# Winter

SBUS_winter_up_ind = my_ind(SBUS_winter_up)
SBUS_winter_shelf_ind = my_ind(SBUS_winter_shelf)

ind = np.concatenate((SBUS_winter_up_ind, SBUS_winter_shelf_ind))
uniqueValues, indicesList = np.unique(ind, return_index=True)

SBUS_winter_up = SBUS_winter_up.drop(SBUS_winter_up.index[uniqueValues.astype(int)])
SBUS_winter_shelf = SBUS_winter_shelf.drop(SBUS_winter_shelf.index[uniqueValues.astype(int)])

###################################
# Step 3: Need to create targets
###################################

# Summer

mylen = len(SBUS_summer_up)

target_summer = np.concatenate((np.full((mylen,1),0),np.full((mylen,1),1)))

# Winter

mylen = len(SBUS_winter_up)

target_winter = np.concatenate((np.full((mylen,1),0),np.full((mylen,1),1)))


###################################
# Step 4: Create our arrays
###################################

# We now want to convert our dataframes to numpy arrays and stack them

# Summer

SBUS_summer_up = SBUS_summer_up.to_numpy()
SBUS_summer_shelf = SBUS_summer_shelf.to_numpy()

SUMMER_ARRAY = np.vstack((SBUS_summer_up,SBUS_summer_shelf)) 

# Winter

SBUS_winter_up = SBUS_winter_up.to_numpy()
SBUS_summer_shelf = SBUS_winter_shelf.to_numpy()

WINTER_ARRAY = np.vstack((SBUS_winter_up,SBUS_winter_shelf))

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
    0:'SBUS UF',
    1:'SBUS SBF',
    }

mycolors = ["deepskyblue","coral"]
myshapes = ['^','o']

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

plt.legend(fontsize=12,loc='upper left')

#fig.savefig('PCA_summer', format='png', dpi=600) 

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

fig.savefig('PCA_SBUS', format='png', dpi=600) 










