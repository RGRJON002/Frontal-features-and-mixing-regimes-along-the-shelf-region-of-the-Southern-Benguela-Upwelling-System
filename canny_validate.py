#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 09:50:20 2022

@author: jono
"""

##### Validate Frontal detection code against OSTIA product

# Read in the needed libraries

import os
os.environ['PROJ_LIB'] = '/home/jono/anaconda3/share/proj/'
from mpl_toolkits.basemap import Basemap
import matplotlib
import cmocean
import cmocean.cm as cmo
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage.util import random_noise
from skimage import feature
import datetime
from netCDF4 import Dataset
import pandas as pd 
from scipy.interpolate import griddata   

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
    return var_CROCO, lon_CROCO, lat_CROCO, mask_CROCO, h_CROCO, pn, pm, time_CROCO

#############################################
# Read in the SST for CROCO 
#############################################

file_base =  '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT/'
Y_min = 2004
Y_max = 2018

# Get lon, lat and topography from the first file

lon_CROCO, lat_CROCO, mask_CROCO, h_CROCO, pn, pm =  extract_CROCO('croco_avg_Y'+str(Y_min)+'M'+str(1)+'.nc',file_base=file_base,var='temp', 
                                 level = 59, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[1:7]

SST = np.empty((0,len(lat_CROCO),len(lon_CROCO)))
TIME = np.empty((0))
for i in range(Y_min, Y_max+1):
    for j in range(1,12+1):
        file_name = 'croco_avg_Y'+str(i)+'M'+str(j)+'.nc'
        if file_name != 'croco_avg_Y2018M12.nc':
            print("Processing "  + file_name)
            nc_file = file_base + file_name
            root = Dataset(nc_file)
            time_CROCO = root.variables['time'][:]
            sst_CROCO = extract_CROCO(file_name,file_base=file_base,var='temp', 
                                  level = 59, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[0]
            SST = np.concatenate((SST, sst_CROCO), axis=0)
            TIME = np.append(TIME,time_CROCO)
        else:
            print('No data')
        
# Any repeated times, remove them

uniqueValues, indicesList = np.unique(TIME, return_index=True)
TIME = TIME[indicesList]
SST = SST[indicesList,:,:]

SST = np.transpose(SST,(1,2,0))  # Fix the array orientation such that time is third dimension

# Convert the time array to a date

TIME = TIME.data   # Masked array appernetly
d0 = datetime.datetime(1990,1,1)
DATE = np.empty((0))
for k in range(0,len(TIME)):
       dt = datetime.timedelta(seconds = TIME[k])
       date  = d0 + dt
       DATE = np.append(DATE,date)
       
###########################################
# Read in the data for OSTIA 
###########################################

file_OSTIA = '/media/data/OBS_DATA/REMMS/OSTIA_1993_2019_mm.nc'
root = Dataset(file_OSTIA)
lon_OSTIA = root.variables['lon'][:]  
lat_OSTIA = root.variables['lat'][:]
time_OSTIA = root.variables['time'][:].astype(float)

# Index OSTIA with CROCO location

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

# Slightly increase domain to account for all points

lat_min_index = find_nearest(lat_OSTIA,min(lat_CROCO))[1]
lat_max_index = find_nearest(lat_OSTIA,max(lat_CROCO))[1]
lon_min_index = find_nearest(lon_OSTIA,min(lon_CROCO))[1]
lon_max_index = find_nearest(lon_OSTIA,max(lon_CROCO))[1]

lon_OSTIA = lon_OSTIA[lon_min_index-2:lon_max_index+2]
lat_OSTIA = lat_OSTIA[lat_min_index-2:lat_max_index+2]

SST_OSTIA = root.variables['analysed_sst'][:,lat_min_index-2:lat_max_index+2,lon_min_index-2:lon_max_index+2]
SST_OSTIA = SST_OSTIA - 273.15   # Convert from Kelvin

# Convert the time array to a date

TIME = time_OSTIA.data   # Masked array appernetly
d0 = datetime.datetime(1981,1,1)
DATE_OSTIA = np.empty((0))
for k in range(0,len(TIME)):
       dt = datetime.timedelta(seconds = TIME[k])
       date  = d0 + dt
       DATE_OSTIA = np.append(DATE_OSTIA,date)
       
# Subset file to time-frame of CROCO, get the appropriate years

Y = np.empty((0))
for dd in range(0,len(DATE_OSTIA)):
    year = DATE_OSTIA[dd].year
    Y = np.append(Y,year)

ind_min_date = min(np.where(Y_min == Y )[0])
ind_max_date = max(np.where(Y_max == Y)[0])

DATE_OSTIA = DATE_OSTIA[ind_min_date:ind_max_date]
SST_OSTIA = SST_OSTIA[ind_min_date:ind_max_date,:,:]
SST_OSTIA = np.transpose(SST_OSTIA,(1,2,0))  # Format time to 3rd dimension

#############################################
# Interpolate OSTIA to CROCO
#############################################
x_in_grid, y_in_grid = np.meshgrid(lon_OSTIA,lat_OSTIA)
x_out, y_out = np.meshgrid(lon_CROCO,lat_CROCO)

SST_interp = np.empty((len(lat_CROCO),len(lon_CROCO),0))
for i in range(0,len(DATE_OSTIA)):
    array = np.ma.masked_invalid(SST_OSTIA[:,:,i])
    x_in = x_in_grid[~array.mask]
    y_in = y_in_grid[~array.mask]
    gr = np.multiply(griddata((x_in, y_in), array[~array.mask].reshape(-1),(x_out, y_out), method='linear'),mask_CROCO)
    SST_interp = np.concatenate((SST_interp, gr[...,np.newaxis]), axis=2)

# Break down date array into months so we can get the seasons
Summer = np.array([1,2,12])
Winter = np.array([6,7,8])

#### CROCO

M = np.empty((0))
for dd in range(0,len(DATE)):
    month = DATE[dd].month
    M = np.append(M,month)
    
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

CROCO_summer = np.nanmean(SST[:,:,ind_summer.astype(int)], axis=2)
CROCO_winter = np.nanmean(SST[:,:,ind_winter.astype(int)], axis=2)

### OSITA

M = np.empty((0))
for dd in range(0,len(DATE_OSTIA)):
    month = DATE_OSTIA[dd].month
    M = np.append(M,month)
    
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

OSTIA_summer = np.nanmean(SST_interp[:,:,ind_summer.astype(int)], axis=2)
OSTIA_winter = np.nanmean(SST_interp[:,:,ind_winter.astype(int)], axis=2)

##################################
# Find the fronts
##################################

# Construct a figure of normalized gradients for CROCO and OSTIA in summer and 
# winter

def normdata(data):
    Z_data = (data - np.nanmin(data))/(np.nanmax(data)-np.nanmin(data))
    return Z_data

def Magnitude(data):
    x_grad,y_grad = np.gradient(data)
    x_grad = x_grad * (pm*1000)
    y_grad = y_grad * (pn*1000)
    mag = np.sqrt(x_grad**2 + y_grad**2)
    return mag

CROCO_sum_grad = Magnitude(CROCO_summer)
CROCO_win_grad = Magnitude(CROCO_winter)

OSTIA_sum_grad = Magnitude(OSTIA_summer)
OSTIA_win_grad = Magnitude(OSTIA_winter)

# Visualise the gradiants

# Summer example

fig, axes = plt.subplots(2,2)

# CROCO Summer
axes[0,0].set_title("(a) CROCO [DJF]",size=13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0,0], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[0,0].pcolormesh(X,Y,CROCO_sum_grad,vmin=0,vmax=0.2,cmap=cmo.thermal, shading='auto')
axes[0,0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')

# OSTIA summer
axes[0,1].set_title("(b) OSTIA [DJF]",size=13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0,1], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[0,1].pcolormesh(X,Y,OSTIA_sum_grad,vmin=0,vmax=0.2,cmap=cmo.thermal, shading='auto')
axes[0,1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')

# CROCO Winter
axes[1,0].set_title("(c) CROCO [JJA]",size=13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[1,0], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[1,0].pcolormesh(X,Y,CROCO_win_grad,vmin=0,vmax=0.2,cmap=cmo.thermal, shading='auto')
axes[1,0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')

# OSTIA summer
axes[1,1].set_title("(d) OSTIA [JJA]",size=13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[1,1], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[1,1].pcolormesh(X,Y,OSTIA_win_grad,vmin=0,vmax=0.2,cmap=cmo.thermal, shading='auto')
axes[1,1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')

plt.tight_layout()
plt.tick_params(axis='both', which='major', labelsize=12)
fig.set_size_inches(8,7)

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axes.flat:
    ax.label_outer()

cbar = fig.colorbar(im, ax=axes.ravel().tolist(),
                    spacing='uniform',
                    orientation='vertical')
cbar.set_label("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",size=12)
plt.show()
fig.savefig('SST_gradient_OSTIAvsCROCO.png', format='png', dpi=600)
plt.close()

# Create a mask for the shelf region to apply the filter

shelf_depth = 500
h_mask = h_CROCO.copy()
h_mask[h_mask > shelf_depth] = float("NAN")
mask_canny = ~np.isnan(h_mask)

###################################
# Find the edges
###################################

edges_sum_CROCO = feature.canny(CROCO_summer,sigma=2,low_threshold=0.8,
                       high_threshold=0.9,mask=mask_canny,use_quantiles=True)

edges_win_CROCO = feature.canny(CROCO_winter,sigma=2,low_threshold=0.8,
                       high_threshold=0.9,mask=mask_canny,use_quantiles=True)

edges_sum_OSTIA = feature.canny(OSTIA_summer,sigma=2,low_threshold=0.8,
                       high_threshold=0.9,mask=mask_canny,use_quantiles=True)

edges_win_OSTIA = feature.canny(OSTIA_winter,sigma=2,low_threshold=0.8,
                       high_threshold=0.9,mask=mask_canny,use_quantiles=True)

## Plot the edges

fig, axes = plt.subplots(2,2)

# CROCO Summer
axes[0,0].set_title("(a) CROCO [DJF]",size=13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0,0], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[0,0].pcolormesh(X,Y,CROCO_summer,vmin=12,vmax=24,cmap=cmo.thermal, shading='auto')
axes[0,0].contour(X,Y,edges_sum_CROCO, linewidths=1.0, linestyles='solid', colors='white')
axes[0,0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='lightgreen')

# OSTIA summer
axes[0,1].set_title("(b) OSTIA [DJF]",size=13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0,1], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[0,1].pcolormesh(X,Y,OSTIA_summer,vmin=12,vmax=24,cmap=cmo.thermal, shading='auto')
axes[0,1].contour(X,Y,edges_sum_OSTIA, linewidths=1.0, linestyles='solid', colors='white')
axes[0,1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='lightgreen')

# CROCO Winter
axes[1,0].set_title("(c) CROCO [JJA]",size=13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[1,0], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[1,0].pcolormesh(X,Y,CROCO_winter,vmin=12,vmax=24,cmap=cmo.thermal, shading='auto')
axes[1,0].contour(X,Y,edges_win_CROCO, linewidths=1.0, linestyles='solid', colors='white')
axes[1,0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='lightgreen')

# OSTIA summer
axes[1,1].set_title("(d) OSTIA [JJA]",size=13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[1,1], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[1,1].pcolormesh(X,Y,OSTIA_winter,vmin=12,vmax=24,cmap=cmo.thermal, shading='auto')
axes[1,1].contour(X,Y,edges_win_OSTIA, linewidths=1.0, linestyles='solid', colors='white')
axes[1,1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='lightgreen')

plt.tight_layout()
plt.tick_params(axis='both', which='major', labelsize=12)
fig.set_size_inches(8,7)

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axes.flat:
    ax.label_outer()

cbar = fig.colorbar(im, ax=axes.ravel().tolist(),
                    spacing='uniform',
                    orientation='vertical')
cbar.set_label("SST "+"[" + u"\N{DEGREE SIGN}"+"C"+"]",size=12)
plt.show()
fig.savefig('SST_fronts_OSTIAvsCROCO.png', format='png', dpi=600)
plt.close()


#%% Some statistics

# Compute the mean SST gradient defined for the fronts in OSTIA vs the model
# Standard deviation. Cannot score RMSD or correaltion as both would have to be
# the same length

# Model

loc_data = np.where(edges_sum_CROCO == 1)
MODEL_SST_FRONTS_sum = CROCO_sum_grad[loc_data[0],loc_data[1]] 
print('mean of Summer model is')
print(np.nanmean(MODEL_SST_FRONTS_sum))
print('std of Summer model is')
print(np.nanstd(MODEL_SST_FRONTS_sum))

loc_data = np.where(edges_win_CROCO == 1)
MODEL_SST_FRONTS_win = CROCO_win_grad[loc_data[0],loc_data[1]] 
print('mean of Winter model is')
print(np.nanmean(MODEL_SST_FRONTS_win))
print('std of Winter model is')
print(np.nanstd(MODEL_SST_FRONTS_win))

# OSTIA

loc_data = np.where(edges_sum_OSTIA == 1)
OSTIA_SST_FRONTS_sum = OSTIA_sum_grad[loc_data[0],loc_data[1]] 
print('mean of Summer obs is')
print(np.nanmean(OSTIA_SST_FRONTS_sum))
print('std of Summer obs is')
print(np.nanstd(OSTIA_SST_FRONTS_sum))

loc_data = np.where(edges_win_OSTIA == 1)
OSTIA_SST_FRONTS_win = OSTIA_win_grad[loc_data[0],loc_data[1]] 
print('mean of Winter obs is')
print(np.nanmean(OSTIA_SST_FRONTS_win))
print('std of Winter obs is')
print(np.nanstd(OSTIA_SST_FRONTS_win))






