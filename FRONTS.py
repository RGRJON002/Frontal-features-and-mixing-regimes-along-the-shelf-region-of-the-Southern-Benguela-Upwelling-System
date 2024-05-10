#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 14:06:18 2022

@author: jono
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 14:21:15 2022

@author: jono
"""

# %% Code to calculate frontal displacement along the SBUS

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

# Function to normalize data
def normdata(data):
    Z_data = (data - np.nanmin(data))/(np.nanmax(data)-np.nanmin(data))
    return Z_data

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

###################################
# HISTOGRAM ENHANCEMENT
###################################

# Histogram Equalization is one of the fundamental tools in the image processing 
# toolkit. It’s a technique for adjusting the pixel values in an image to enhance
# the contrast by making those intensities more equal across the board. Typically,
# the histogram of an image will have something close to a normal distribution, 
# but equalization aims for a uniform distribution.

# As we have a SST matrix, we will convert it to pixel intensity and then run our
# Canny detection

SST_enh = normdata(SST[:,:,8])     # Normalize SST
SST_enh = np.floor(SST_enh * 255)  # Convert to pixel

from skimage import exposure

img_enh = exposure.equalize_hist(SST_enh.data,nbins=256,mask=mask_canny)
edges_enh = feature.canny(img_enh,sigma=2,low_threshold=0.8,
                       high_threshold=0.9,mask=mask_canny, use_quantiles=True)

# Convert the enhanced image back to SST  

SST_enh = img_enh*(np.nanmax(SST[:,:,8])-np.nanmin(SST[:,:,8])) + np.nanmin(SST[:,:,8])

# Overlay the enhanced edges on the enhanced image and standard

fig, axes = plt.subplots(1,2,figsize=(12,5))
fig.tight_layout()

# Normal SST

axes[0].set_title("SST normal: sigma =2, quantiles")
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0], resolution='l')
map.fillcontinents(color='grey')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[0].pcolormesh(X,Y,SST[:,:,8],vmin=12,vmax=22,cmap=cmo.thermal, shading='auto')
axes[0].contour(X, Y, edges_enh, linewidths=0.1, linestyles='dashed', colors='black')
axes[0].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')

# Enhanced SST
axes[1].set_title("SST enhanced: sigma =2, quantiles")
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[1], resolution='l')
map.fillcontinents(color='grey')
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
im = axes[1].pcolormesh(X,Y,SST_enh,vmin=12,vmax=22,cmap=cmo.thermal, shading='auto')
axes[1].contour(X, Y, edges_enh, linewidths=0.1, linestyles='dashed', colors='black')
axes[1].contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')

plt.show()
fig.savefig('SST_enhanced.png', format='png', dpi=600)
plt.close()

# The fronts become more connected in their structures, will not implement yet
# But may be good idea going forward if required

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
pcm = map.pcolormesh(x,y,SST[:,:,8]*sb_mask,vmin=12,vmax=22,cmap=cmo.thermal, shading='auto')

#%%
# Loop over the SST and get the edges

def contour_loc(contour,lon,lat):
    Xloc = np.empty((0))
    Yloc = np.empty((0))
    for k in range(0,len(contour)):
        lat_index = find_nearest(lat,contour[k,1])[1]
        lon_index = find_nearest(lon,contour[k,0])[1]
        Yloc = np.append(Yloc,lat_index)
        Xloc = np.append(Xloc,lon_index)
    return(Xloc,Yloc)

EDGE = np.empty((len(lat_CROCO),len(lon_CROCO),0))

for k in range(0, len(TIME)):
    edge = feature.canny(SST[:,:,k],sigma=2,low_threshold=0.8,
                           high_threshold=0.9,mask=mask_canny,use_quantiles=True)
    EDGE = np.concatenate((EDGE, edge[...,np.newaxis]), axis=2)
    
# Loop over the SST and get the gradient of SST

GRADIENT = np.empty((len(lat_CROCO),len(lon_CROCO),0))

for k in range(0, len(TIME)):
    x_grad,y_grad = np.gradient(SST[:,:,k])
    x_grad = x_grad * (pm*1000)
    y_grad = y_grad * (pn*1000)
    sst_grad = np.sqrt(x_grad**2 + y_grad**2)
    GRADIENT = np.concatenate((GRADIENT, sst_grad[...,np.newaxis]), axis=2)
    
#%%

############################
# SAVE AND LOAD VARIABLES
############################

# SAVE

import pickle

f = open('GRADIENT.pckl', 'wb')
pickle.dump(GRADIENT, f)
f.close()

f = open('EDGE.pckl', 'wb')
pickle.dump(EDGE, f)
f.close()

# LOAD

f = open('GRADIENT.pckl', 'rb')
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

# Check contours in relation to masks

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
map.pcolormesh(x,y,GRADIENT[:,:,95]*up_mask*SH_mask,vmin=0,vmax=0.3,cmap=cmo.thermal, shading='auto')
pcm = map.contour(X, Y, EDGE[:,:,95]*up_mask*SH_mask, linewidths=1, linestyles='dashed', colors='green')
map.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')

#%%

# Loop over the time-domain 

from geopy import distance

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
#W_shelf = []
for i in range(0,len(TIME)): 
    loc_data = np.where(EDGE[:,:,i]*sb_mask*SH_mask == 1)
#    W_shelf.append(np.nanmean(W[loc_data[0],loc_data[1],i]))  # Vertical velocity
    GRAD_shelf.append(np.nanmean(GRADIENT.data[loc_data[0],loc_data[1],i]))  # Pixel gradient
    pcm = map.contour(X, Y, EDGE[:,:,i]*sb_mask*SH_mask, linewidths=0.1, linestyles='solid', colors='black')
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
#W_up = []
for i in range(0,len(TIME)): 
    loc_data = np.where(EDGE[:,:,i]*up_mask*SH_mask == 1)
#    W_up.append(np.nanmean(W[loc_data[0],loc_data[1],i]))  # Vertical velocity
    GRAD_up.append(np.nanmean(GRADIENT.data[loc_data[0],loc_data[1],i]))  # Pixel gradient
    pcm = map.contour(X, Y, EDGE[:,:,i]*up_mask*SH_mask, linewidths=0.1, linestyles='solid', colors='black')
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
    
#%%    
    
# Create some basic time-series to visualize the data

# Dsiplacement

fig = plt.figure(1, figsize=(12,5))
axes  = fig.add_subplot(111)
axes.set_title("Time-series of frontal displacement relative to coast",fontsize=16)
plt.xlabel("DATE",fontsize=16)
plt.ylabel("Distance [km]",fontsize=16)
plt.plot(DATE, DIS_shelf)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()

fig.savefig('Time_sereis_DIS', format='png', dpi=600)
plt.close()

# Length

fig = plt.figure(1, figsize=(12,5))
axes  = fig.add_subplot(111)
axes.set_title("Time-series of front length",fontsize=16)
plt.xlabel("DATE",fontsize=16)
plt.ylabel("Length [pixels]",fontsize=16)
plt.plot(DATE, LEN_shelf)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()

fig.savefig('Time_sereis_LEN', format='png', dpi=600)
plt.close()

# Number of detected fronts

fig = plt.figure(1, figsize=(12,5))
axes  = fig.add_subplot(111)
axes.set_title("Time-series of detected fronts",fontsize=16)
plt.xlabel("DATE",fontsize=16)
plt.ylabel("Count",fontsize=16)
plt.plot(DATE, NUM_shelf)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()

fig.savefig('Time_sereis_NUM', format='png', dpi=600)
plt.close()

# Gradient

fig = plt.figure(1, figsize=(12,5))
axes  = fig.add_subplot(111)
axes.set_title("Time-series of frontal pixel gradeint",fontsize=16)
plt.xlabel("DATE",fontsize=16)
plt.ylabel("Gradient [pixel]",fontsize=16)
plt.plot(DATE, GRAD_shelf)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()

fig.savefig('Time_sereis_GRAD', format='png', dpi=600)
plt.close()

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
        'GRAD': GRAD_shelf}
df_shelf = pd.DataFrame(data) 
df_shelf = df_shelf.set_index('Date')

###################################
# UPWELLING data frame
###################################

data = {'Date': DATE,
        'DIS': DIS_up,
        'LEN': LEN_up,
        'NUM': NUM_up,
        'GRAD': GRAD_up}
df_up = pd.DataFrame(data) 
df_up = df_up.set_index('Date')   


#%%

############################################
# Plot combining displacmenet and gradient
############################################

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection 

# Compute the monthly means for each variable

data_columns = ['DIS', 'LEN', 'NUM','GRAD']

mon_shelf = df_shelf[data_columns].resample('M').mean()
day3 = df_shelf[data_columns].rolling(3).mean()
mon_up = df_up[data_columns].resample('M').mean()

##########################
# MAIN FIGURE
##########################

fig, axes = plt.subplots(2,1,sharex=True, figsize=(12,5))

# Shelf-break
dates = DATE
y = np.array(day3.DIS)
c = np.array(day3.GRAD)
             
axes[0].set_title("(a)",loc='left', fontsize=16)
#convert dates to numbers first
inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="viridis", linewidth=1)
# Add the monthly values
axes[0].scatter(x=mon_shelf.index, y=mon_shelf['DIS'],s=90,c=mon_shelf['GRAD'],cmap='viridis')
# set color to date values
lc.set_array(c)
line = axes[0].add_collection(lc)
axes[0].tick_params(axis="y", labelsize=14) 
axes[0].grid()
axes[0].set_ylabel("Distance [km]",fontsize=16)
axes[0].set_ylim([60, 150])

# Upwelling
dates = DATE
y = np.array(DIS_up)
c = np.array(GRAD_up)

axes[1].set_title("(b)",loc='left', fontsize=16)
#convert dates to numbers first
inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="viridis", linewidth=1)
# Add the monthly values
axes[1].scatter(x=mon_up.index, y=mon_up['DIS'],s=90,c=mon_up['GRAD'],cmap='viridis')
# set color to date values
lc.set_array(c)
line = axes[1].add_collection(lc)
axes[1].grid()
axes[1].set_xlabel("DATE",fontsize=16)
axes[1].set_ylabel("Distance [km]",fontsize=16)
axes[1].set_ylim([20, 100])
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a
plt.colorbar(line, ax=axes.ravel().tolist()).set_label(label="Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",size=13)
axes[1].xaxis_date()

plt.show()
fig.savefig('SH_Fronts', format='png', dpi=600)


#%%

clim = df_shelf.groupby([df_shelf.index.month]).mean() # CLIM

anom = df_shelf.groupby(df_shelf.index.month).transform(lambda x: x-x.mean()) # ANOM

#%%
##########################################################
# GRADINET vs vertical velocity
###########################################################

##########################
# MAIN FIGURE
##########################

fig, axes = plt.subplots(2,1,sharex=True, figsize=(12,5))

# Shelf-break
dates = DATE
y = np.array(GRAD_shelf)
c = np.array(W_shelf)

axes[0].set_title("(a)",loc='left', fontsize=16)
#convert dates to numbers first
inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="RdBu_r", linewidth=1)
# Add the monthly values
axes[0].scatter(x=mon_shelf.index, y=mon_shelf['GRAD'],s=90,c=mon_shelf['W'],cmap='RdBu_r')
# set color to date values
lc.set_array(c)
line = axes[0].add_collection(lc)
axes[0].tick_params(axis="y", labelsize=14) 
axes[0].grid()
axes[0].set_ylabel("Gradient",fontsize=16)
#axes[0].set_ylim([60, 150])

# Upwelling
dates = DATE
y = np.array(GRAD_up)
c = np.array(W_up)

axes[1].set_title("(b)",loc='left', fontsize=16)
#convert dates to numbers first
inxval = mdates.date2num(dates)
points = np.array([inxval, y]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1],points[1:]], axis=1)
lc = LineCollection(segments, cmap="RdBu_r", linewidth=1)
# Add the monthly values
axes[1].scatter(x=mon_up.index, y=mon_up['GRAD'],s=90,c=mon_up['W'],cmap='RdBu_r')
# set color to date values
lc.set_array(c)
line = axes[1].add_collection(lc)
axes[1].grid()
axes[1].set_xlabel("DATE",fontsize=16)
axes[1].set_ylabel("Gradient",fontsize=16)
#axes[1].set_ylim([20, 100])
plt.tick_params(axis='both', which='major', labelsize=14)

# Colorbar a

plt.colorbar(line, ax=axes.ravel().tolist(), label='Vertical velocities [m/day]' )
axes[1].xaxis_date()
plt.show()

fig.savefig('Time_sereis_GRADvsW', format='png', dpi=600)
plt.close()


# Compute the correlation

from scipy import stats
# Y and Z are numpy arrays or lists of variables
ind_nan = np.argwhere(np.isnan(DIS_shelf))
a = np.delete(WIND_shelf, ind_nan)
b = np.delete(DIS_shelf, ind_nan)
stats.pearsonr(month['Dis'], month['Wind'])

#%%

# Segment to create a small gif of SST and EDGES

# Define the meta data for the movieimport numpy as np
import matplotlib.animation as manimation
import os

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='SST_gif', artist='Matplotlib',
                comment='Animation of SST fronts along the SBUS')
writer = FFMpegWriter(fps=0.5, metadata=metadata)

# Initialize the movie
fig = plt.figure()

with writer.saving(fig, "SBUS_SST.mp4", 100):
    for i in range(0,len(TIME),5):
        map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
                    llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
        map.fillcontinents(color='grey')
        x,y = map(lon_CROCO,lat_CROCO)
        X, Y = np.meshgrid(x, y)
        plt.title("SST" + " " + str(DATE[i]))
        plt.pcolormesh(X,Y,SST[:,:,i], vmin=12, vmax = 22, cmap=matplotlib.cm.nipy_spectral) 
        plt.clim(12, 22)
        cbar = plt.colorbar()     
        cbar.set_label(r'SST')
        plt.contour(X, Y, EDGE[:,:,i], linewidths=0.1, linestyles='dashed', colors='black')
        plt.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
        writer.grab_frame()
        plt.clf()

import imageio
import os, sys

class TargetFormat(object):
    GIF = ".gif"
    MP4 = ".mp4"
    AVI = ".avi"

def convertFile(inputpath, targetFormat):
    """Reference: http://imageio.readthedocs.io/en/latest/examples.html#convert-a-movie"""
    outputpath = os.path.splitext(inputpath)[0] + targetFormat
    print("converting\r\n\t{0}\r\nto\r\n\t{1}".format(inputpath, outputpath))

    reader = imageio.get_reader(inputpath)
    fps = reader.get_meta_data()['fps']

    writer = imageio.get_writer(outputpath, fps=fps)
    for i,im in enumerate(reader):
        sys.stdout.write("\rframe {0}".format(i))
        sys.stdout.flush()
        writer.append_data(im)
    print("\r\nFinalizing...")
    writer.close()
    print("Done.")

convertFile("/home/jono/SBUS_SST.mp4", TargetFormat.GIF)

# Upwelling front

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
        llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=axes[0], resolution='l')
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
displacement = []
for i in range(0,len(TIME)):    
    pcm = map.contour(X, Y, EDGE[:,:,i]*up_mask, linewidths=0.1, linestyles='solid', colors='black')
    pcm = pcm.allsegs[0]
    dis_mean = []
    for k in range(0,len(pcm)):
        cc = pcm[k]
        Xloc, Yloc = contour_loc(cc,x,y)
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
    displacement.append(np.mean(dis_mean))


###################################################################
# Calculate mean displacement of the front for summer and winter
###################################################################

# Apply the Canny algorithm to mean SST for summer

edges_sum_shelf = feature.canny(np.nanmean(SST[:,:,ind_summer.astype(int)], axis = 2),sigma=2,low_threshold=0.8,
                       high_threshold=0.9,mask=mask_canny, use_quantiles=True)

# Quick view
fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
X, Y = np.meshgrid(x, y)
im = map.contourf(X,Y,np.nanmean(SST[:,:,ind_summer.astype(int)], axis = 2),vmin=12,vmax=22,cmap=cmo.thermal, shading='auto')
cbar = plt.colorbar(im)
cbar.set_label(u"\N{DEGREE SIGN}"+"C")
map.contour(X, Y, edges_sum_shelf, linewidths=0.1, linestyles='dashed', colors='black')
map.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
plt.colorbar(im, ax=axes.ravel().tolist(),
                    extend='both',
                    extendfrac='auto',
                    ticks=bounds,
                    spacing='uniform',
                    orientation='vertical')

cbar.set_label(u"\N{DEGREE SIGN}"+"C")
    
# Get the shelf-break isotherm summer

# Consider it as the average isotherm between the 200 and 500 m isobaths

shelf_iso_sum = np.nanmean(SST[:,:,ind_summer.astype(int)], axis = 2)*sb_mask
shelf_iso_sum[shelf_iso_sum == 0] = float("NAN")
shelf_iso_sum = np.nanmean(shelf_iso_sum)

# According to Lamont et al. (2015): position is 17 summer, 15, winter. 
# Accounting for temperature bias of ~0.6-1 deg, good approximation

# Get the shelf-break isotherm summer 

shelf_iso_win = np.nanmean(SST[:,:,ind_winter.astype(int)], axis = 2)*sb_mask
shelf_iso_win[shelf_iso_win == 0] = float("NAN")
shelf_iso_win = np.nanmean(shelf_iso_win)

# Store the location of the contours

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
up_iso = map.contour(X, Y, np.nanmean(SST[:,:,ind_summer.astype(int)], axis = 2),levels=[shelf_iso_sum],
            linewidths=2.5, linestyles='dashed', colors='yellow')
sb_iso= map.contour(X, Y, np.nanmean(SST[:,:,ind_winter.astype(int)], axis = 2),levels=[shelf_iso_win],
            linewidths=2.5, linestyles='dashed', colors='green')

# Loop over the contours and combine them together as a single array and plot

up_iso = up_iso.allsegs
plt.show()
fig.savefig('SST_W_canny.png', format='png', dpi=600)
plt.close()
X_up = np.empty((0))
Y_up = np.empty((0))
for i in range(0,len(up_iso[0])):
    Xloc, Yloc = contour_loc(up_iso[0][i],lon_CROCO,lat_CROCO)
    X_up = np.append(X_up,Xloc)
    Y_up = np.append(Y_up,Yloc)
    
sb_iso = sb_iso.allsegs

X_sb = np.empty((0))
Y_sb = np.empty((0))
for i in range(0,len(sb_iso[0])):
    Xloc, Yloc = contour_loc(sb_iso[0][i],lon_CROCO,lat_CROCO)
    X_sb = np.append(X_sb,Xloc)
    Y_sb = np.append(Y_sb,Yloc)
    


# Get the upwelling isotherm

Cs = plt.contour(x, y, edges_mean, linewidths=0.1, linestyles='solid', colors='black')
mycontours = Cs.allsegs[0]

upwelling_iso = mycontours[14]

def contour_loc(contour,lon,lat):
    Xloc = np.empty((0))
    Yloc = np.empty((0))
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx], idx

    for k in range(0,len(contour)):
        lat_index = find_nearest(lat,contour[k,1])[1]
        lon_index = find_nearest(lon,contour[k,0])[1]
        Yloc = np.append(Yloc,lat_index)
        Xloc = np.append(Xloc,lon_index)
    return(Xloc,Yloc)

Xloc, Yloc = contour_loc(upwelling_iso,lon_CROCO,lat_CROCO)
upwelling_iso = np.nanmean(sst_mean[Yloc.astype(int),Xloc.astype(int)])

# According to Lamont et al. (2015) Columbine described by 15 deg isotherm

# Plot the contours

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
            llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
im = map.pcolormesh(x,y,sst_mean,vmin=12,vmax=22,cmap=cmo.thermal, shading='auto')
cbar = plt.colorbar(im)
cbar.set_label(u"\N{DEGREE SIGN}"+"C")
map.contour(X, Y, edges_mean, linewidths=0.1, linestyles='dashed', colors='black')
map.contour(X, Y, sst_mean,levels=[upwelling_iso],
            linewidths=2.5, linestyles='dashed', colors='yellow')
map.contour(X, Y, sst_mean,levels=[shelf_break_iso],
            linewidths=2.5, linestyles='dashed', colors='green')
map.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='white')
plt.colorbar(im, ax=axes.ravel().tolist(),
                    extend='both',
                    extendfrac='auto',
                    ticks=bounds,
                    spacing='uniform',
                    orientation='vertical')
fig.show()
cbar.set_label(u"\N{DEGREE SIGN}"+"C")
fig.savefig('SST_fronts.png', dpi=300)
plt.close()

# Add the contours to a reference grid called newgrid, here is where the contours will be found
    
newgrid = np.zeros([np.size(edges_orig,0),np.size(edges_orig,1)])

newgrid[Y_sb.astype(int),X_sb.astype(int)] = 2
newgrid[Y_up.astype(int),X_up.astype(int)] = 1