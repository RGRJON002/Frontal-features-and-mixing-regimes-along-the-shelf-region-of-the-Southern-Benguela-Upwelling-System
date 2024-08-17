#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 10:29:47 2024

@author: jono
"""

# This script computes frontal EDGES of only significant fronts using ABSOLUTE thresholds
# instead of quantiles. This does prevent the study of frontal diagnostics in multiple seasons
# But 
# Order of CANNY output_mask, isobel, jsobel, magnitude

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

#%% Read in the SST and vertical velcocity

file_base =  '/media/jono/SBUS/SBUS_3km/CHPC_OUTPUT/'
Y_min = 2004
Y_max = 2018

# Get lon, lat and topography from the first file

lon_CROCO, lat_CROCO, h_CROCO, pn, pm =  extract_CROCO('croco_avg_Y'+str(Y_min)+'M'+str(1)+'.nc',file_base=file_base,var='temp', 
                                 level = 59, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[1:6]

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
            print("No Data")
            
# Any repeated times, remove them

uniqueValues, indicesList = np.unique(TIME, return_index=True)
TIME = TIME[indicesList]
SST = SST[indicesList,:,:]      
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

SST = np.transpose(SST,(1,2,0))  # Fix the array orientation such that time is third dimension

sst_summer = np.nanmean(SST[:,:,ind_summer.astype(int)], axis=2)
sst_winter = np.nanmean(SST[:,:,ind_winter.astype(int)], axis=2)
sst_all = np.mean(SST,axis=2)

# Create a mask for the shelf region to apply the filter

shelf_depth = 500
h_mask = h_CROCO.copy()
h_mask[h_mask > shelf_depth] = float("NAN")
mask_canny = ~np.isnan(h_mask)

#%%
####################################
# IMPROVE the CANNY algorithm
####################################

# The Canny edge detection algorithm has three parameters that have to be set:
# all of which will have drastically different effects on the ability to detect
# genuine fronts and answer the scientific question investigating the migration
# of fronts along the SBUS. 
# The three parameters are: 
# 1) The vlue of sigma for the gaussian filter
# 2) The threshholds used for detection
# 3) The usage of a mask and quantiles
# All these need to be addressed to create the best and most suitable version
# of Canny for the SBUS

######################################
# VALUE of SIGMA
######################################

# G = exp(-(x^2+y^2)/2sigma^2) is the Canny filter
# The choice of sigma defines the length scale over which we want to define our
# fronts. From Oram et al. (2008), the choice of sigma is dependant upon the
# desired edge scale relative to the resolution of the image. So, for low values
# of sigma, we would be looking for fine scale features whereas a high value of sigma
# will look for broader features. Incidentally, the choice of sigma will impact
# the ability to detect desired features in an image. Therefore, will test some
# values of sigma to see how it influences front detection, both in summer, and winter

# Normalize gradients for a summer and winter example

from scipy.ndimage import gaussian_filter

# Gradients are more intense in summer, so use a summer example

#Summer

x_grad,y_grad = np.gradient(SST[:,:,ind_summer[1].astype(int)])
x_grad = x_grad * (pm*1000)
y_grad = y_grad * (pn*1000)
sst_grad = np.sqrt(x_grad**2 + y_grad**2)
sst_grad_sum = sst_grad

fig, axes = plt.subplots(2,2)

# Original
axes[0,0].set_title("(a) No filter", size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0,0], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[0,0].pcolormesh(X,Y,sst_grad_sum,vmin=0,vmax=0.2,cmap=cmo.thermal, shading='auto')
axes[0,0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')

# Low filter
axes[0,1].set_title("(b) Sigma = 1",size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0,1], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[0,1].pcolormesh(X,Y,gaussian_filter(sst_grad_sum,sigma=1),vmin=0,vmax=0.2,cmap=cmo.thermal, shading='auto')
axes[0,1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

# Moderate filter
axes[1,0].set_title("(c) Sigma = 3", size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[1,0], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[1,0].pcolormesh(X,Y,gaussian_filter(sst_grad_sum,sigma=3),vmin=0,vmax=0.2,cmap=cmo.thermal, shading='auto')
axes[1,0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

# High filter

axes[1,1].set_title("(d) Sigma = 5", size = 13)
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[1,1], resolution='l')
map.fillcontinents(color='grey')
#map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawparallels([-34., -32., -30.],labels=[1,0,0,0]) # draw parallels
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[1,1].pcolormesh(X,Y,gaussian_filter(sst_grad_sum,sigma=5),vmin=0,vmax=0.2,cmap=cmo.thermal, shading='auto')
axes[1,1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
map.fillcontinents(color='grey')

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
fig.savefig('Summer_exampleofsigma.png', format='png', dpi=600)
plt.close()

# The results: as sigma increases, The number of fronts detected, decreases but
# fronts become longer and more connected. Using a sigma value of 3 means the 
# identification of frontal features will have a length scale of 9 km, which is
# a reasonable estimate concerning the fact that fronts and filaments are mesoscale
# features. But, 3 seems to blur a lot out, so maybe sigma 2 is a good comprimise
#%%

############################
# Thresholds of Detection
############################

# Much like sigma, the value of the thresholds influence the detection of the gradiets
# Adopting the approach by Oram et al. (2008), the quantiles are used to compute, the 
# gradient.

# Create the gradient plots for summer and winter
#%%

#Summer

x_grad,y_grad = np.gradient(sst_summer)
x_grad = x_grad * (pm*1000)
y_grad = y_grad * (pn*1000)
sst_grad = np.sqrt(x_grad**2 + y_grad**2)
sst_grad_sum = sst_grad.flatten()

#Winter

x_grad,y_grad = np.gradient(sst_winter)
x_grad = x_grad * (pm*1000)
y_grad = y_grad * (pn*1000)
sst_grad = np.sqrt(x_grad**2 + y_grad**2)
sst_grad_win = sst_grad.flatten()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

# Plot showing the detection values: important to have

fig, axes = plt.subplots(1,2,figsize=(12,5))

axes[0].set_title('(a) Summer [DJF]',size=13)
axes[0].set_ylabel('Normalized Cumulative Frequency',size=12)
axes[0].set_xlabel("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",size=12)
gr = axes[0].hist(sst_grad_sum, bins=100, density=1, histtype='step', cumulative=True)

val_0_8 = find_nearest(gr[0],0.8)
val_0_9 = find_nearest(gr[0],0.9) 

axes[0].hlines(y=val_0_8[0], xmin=0, xmax=gr[1][val_0_8[1].astype(int)], linewidth=2, color='r',
               label = ('Lower Threshold='+str(round(gr[1][val_0_8[1].astype(int)],3))))
axes[0].vlines(x=gr[1][val_0_8[1].astype(int)],ymin=0, ymax = gr[0][val_0_8[1].astype(int)],
               linewidth=2, color='r')
axes[0].hlines(y=val_0_9[0], xmin=0, xmax=gr[1][val_0_9[1].astype(int)], linewidth=2, color='g',
               label = ('Upper Threshold='+str(round(gr[1][val_0_9[1].astype(int)],3))))
axes[0].vlines(x=gr[1][val_0_9[1].astype(int)],ymin=0, ymax = gr[0][val_0_9[1].astype(int)],linewidth=2, color='g')

leg = axes[0].legend()
axes[0].tick_params(axis='both', which='major', labelsize=12)

axes[1].set_title('(b) Winter [JJA]',size=13)
axes[1].set_ylabel('Normalized Cumulative Frequency',size=12)
axes[1].set_xlabel("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",size=12)
gr = axes[1].hist(sst_grad_win, bins=100, density=1, histtype='step', cumulative=True)

val_0_8 = find_nearest(gr[0],0.8)
val_0_9 = find_nearest(gr[0],0.9) 

axes[1].hlines(y=val_0_8[0], xmin=0, xmax=gr[1][val_0_8[1].astype(int)], linewidth=2, color='r',
               label = ('Lower Threshold='+str(round(gr[1][val_0_8[1].astype(int)],3))))
axes[1].vlines(x=gr[1][val_0_8[1].astype(int)],ymin=0, ymax = gr[0][val_0_8[1].astype(int)],
               linewidth=2, color='r')
axes[1].hlines(y=val_0_9[0], xmin=0, xmax=gr[1][val_0_9[1].astype(int)], linewidth=2, color='g',
               label = ('Upper Threshold='+str(round(gr[1][val_0_9[1].astype(int)],3))))
axes[1].vlines(x=gr[1][val_0_9[1].astype(int)],ymin=0, ymax = gr[0][val_0_9[1].astype(int)],linewidth=2, color='g')

leg = axes[1].legend()
plt.tick_params(axis='both', which='major', labelsize=12)
plt.show()

# Calculate values

fig.savefig('Quantiles.png', format='png', dpi=600)
plt.close() 
#%%
###############################################################################
# DO the WHILE YEAR average
###############################################################################

# ALL

x_grad,y_grad = np.gradient(sst_all)
x_grad = x_grad * (pm*1000)
y_grad = y_grad * (pn*1000)
sst_grad = np.sqrt(x_grad**2 + y_grad**2)
sst_grad_all = sst_grad.flatten()

fig, axes = plt.subplots(1,1,figsize=(5,5))

axes.set_title('Annual Mean',size=13)
axes.set_ylabel('Normalized Cumulative Frequency',size=12)
axes.set_xlabel("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",size=12)
gr = axes.hist(sst_grad_all, bins=100, density=1, histtype='step', cumulative=True)

val_0_8 = find_nearest(gr[0],0.8)
val_0_9 = find_nearest(gr[0],0.9) 

axes.hlines(y=val_0_8[0], xmin=0, xmax=gr[1][val_0_8[1].astype(int)], linewidth=2, color='r',
               label = ('Lower Threshold='+str(round(gr[1][val_0_8[1].astype(int)],3))))
axes.vlines(x=gr[1][val_0_8[1].astype(int)],ymin=0, ymax = gr[0][val_0_8[1].astype(int)],
               linewidth=2, color='r')
axes.hlines(y=val_0_9[0], xmin=0, xmax=gr[1][val_0_9[1].astype(int)], linewidth=2, color='g',
               label = ('Upper Threshold='+str(round(gr[1][val_0_9[1].astype(int)],3))))
axes.vlines(x=gr[1][val_0_9[1].astype(int)],ymin=0, ymax = gr[0][val_0_9[1].astype(int)],linewidth=2, color='g')

leg = axes.legend()
axes.tick_params(axis='both', which='major', labelsize=12)

fig.savefig('ABSOLUTE.png', format='png', dpi=600)
plt.close() 

# From this, the peak change in gradient occurs over the 80 and 90th percentile
# Compare the new threshold detection relative to the standard with a sigma of 2
# Used quantiles to create new Canny detection
# Apply quantiles and use threshold of [0.8,0.9] with a sigma value of 2
# The agrrement between the SST fronts and vertical velocities are quite good, but,
# as stated by Miltadou et al. a short coming of the Canny Edge detection is
# that it can underestimate front detection. To improve, the detection of frontal
# features, the image can be enhanced using a histogram equalisation image
# enhancement

#%%
#####################################################
# MAIN PROGRAM
#####################################################

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

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

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

# Ckeck the mask

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
pcm = map.pcolormesh(x,y,SST[:,:,8]*mask_canny,vmin=12,vmax=22,cmap=cmo.thermal, shading='auto')

#%%
# Loop over the SST and get the edges
# As we ae defining absolute threshold values, they are as follows: 0.05 for lower and 
# 0.099 for upper 

def contour_loc(contour,lon,lat):
    Xloc = np.empty((0))
    Yloc = np.empty((0))
    for k in range(0,len(contour)):
        lat_index = find_nearest(lat,contour[k,1])[1]
        lon_index = find_nearest(lon,contour[k,0])[1]
        Yloc = np.append(Yloc,lat_index)
        Xloc = np.append(Xloc,lon_index)
    return(Xloc,Yloc)

EDGE_ABSOLUTE = np.empty((len(lat_CROCO),len(lon_CROCO),0))

for k in range(0, len(TIME)):
    edge= feature.canny(SST[:,:,k],sigma=2,low_threshold=0.05,
                           high_threshold=0.099,mask=mask_canny,use_quantiles=False)[0]
    
    EDGE_ABSOLUTE = np.concatenate((EDGE_ABSOLUTE, edge[...,np.newaxis]), axis=2)
    
# Loop over the SST and get the gradient of SST

#GRADIENT = np.empty((len(lat_CROCO),len(lon_CROCO),0))

#for k in range(0, len(TIME)):
#    x_grad,y_grad = np.gradient(SST[:,:,k])
#    x_grad = x_grad * (pm*1000)
#    y_grad = y_grad * (pn*1000)
#    sst_grad = np.sqrt(x_grad**2 + y_grad**2)
#    GRADIENT = np.concatenate((GRADIENT, sst_grad[...,np.newaxis]), axis=2)
    
#%%

############################
# SAVE AND LOAD VARIABLES
############################

# SAVE

import pickle

f = open('EDGE_ABS.pckl', 'wb')
pickle.dump(EDGE_ABSOLUTE, f)
f.close()

# LOAD

f = open('/media/data/DIAGNOSTICS/FRONTS/GRADIENT.pckl', 'rb')
GRADIENT = pickle.load(f)
f.close()

f = open('EDGE.pckl', 'rb')
EDGE = pickle.load(f)
f.close()

#%%

####################
# CLEAR SOME SPACE
#####################

# SST is very large so need to clear it to process the other variables

del SST, sst_CROCO 

#%%
###################################
# CHECK CONTOURS
###################################

ind = 5

# Check contours in relation to masks

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
im = map.pcolormesh(x,y,GRADIENT[:,:,ind]*mask_canny,vmin=0,vmax=0.3,cmap=cmo.thermal, shading='auto')
pcm = map.contour(X, Y, EDGE_ABSOLUTE[:,:,ind]*mask_canny, linewidths=1, linestyles='dashed', colors='green')
map.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
cbar = fig.colorbar(im,
                    spacing='uniform',
                    orientation='vertical')
cbar.set_label("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",size=12)
plt.show()

