#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 16:02:49 2022

@author: jono
"""

# This script aims to extract the nearshore surface currents and see if we have any near shore reversal during the wind events

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

V = np.empty((0,len(lat_CROCO),len(lon_CROCO)))
TIME = np.empty((0))
for i in range(Y_min, Y_max+1):
    for j in range(1,12+1):
        file_name = 'croco_avg_Y'+str(i)+'M'+str(j)+'.nc'
        if file_name != 'croco_avg_Y2018M12.nc':
            print("Processing "  + file_name)
            nc_file = file_base + file_name
            root = Dataset(nc_file)
            time_CROCO = root.variables['time'][:]            
            v = extract_CROCO(file_name,file_base=file_base,var='v', 
                                      level = 59, lat_min = -36., lat_max = -28., lon_min = 15.,  lon_max = 20.)[0]
            TIME = np.append(TIME,time_CROCO)
            V = np.concatenate((V, v), axis=0)
        else:
            print("No Data")
            
# Any repeated times, remove them

uniqueValues, indicesList = np.unique(TIME, return_index=True)
TIME = TIME[indicesList]    
V = V[indicesList,:,:]

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

V = np.transpose(V,(1,2,0))  # Fix the array orientation such that time is third dimension

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

############################################
# Get the near-coastal meridional currents
############################################

V_COAST = np.empty((len(lat_CROCO),0))
for j in range(0,len(TIME)):
    v_coast = []
    for i in range(0,len(ind_coast_X)):
        if np.isnan(ind_coast_X[i]):
            out = float("NAN")
            v_coast.append(out)  # Vertical velocity
            #print("NAN")
        else:
            ind_x = list(range(ind_coast_X[i].astype(int)-4,ind_coast_X[i].astype(int)-1))
            out = np.nanmean(V[i,ind_x,j])
            v_coast.append(out)
    v_coast = np.array(v_coast)        
    V_COAST = np.concatenate((V_COAST, v_coast[...,np.newaxis]), axis=1)
   
#######################
# SAVE the V data
#######################

#f = open('V.pckl', 'wb')
#pickle.dump(V_COAST, f)
#f.close()    

#######################
# LOAD the V File
#######################

f = open('V.pckl', 'rb')
V_COAST = pickle.load(f)
f.close()

del V
del v

#%%

#################################
# STRATIFICATION
#################################
 
# We want to know how stratification affect frontal dynamics. Herein we make some assumptions. We assume the open ocean
# is always less stratified and constant whereas the stratification under the frontal edges is representative of the stratification
# along the shelf. The shelf is distinguished between the UF and SBF to investigate them differently.
# Some basics of the metric, we calculate the BVF in the offshore domain, 500 km from the 500 m isobath and mean it as the offshore
# The offshore BVF will be BVF_off
# Next, the BVF values need to be extracted from the fontal positions. No need to assign to different regions, only interested in the value
# We will have two: BVF_up for the upwelling region and BVF_shelf for the shelf_break
# We could in fact average these two values of BVF_shelf and BVF_in but in general: we want to subtract from it
# the value of BVF_off to see the difference along the shelf region and whether significant changes in stratiification influence the \
# current dynamics

######################################################################
# STEP 1: Calculate offshore BVF_off as well as BVF_in and BVF_up
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

BVF_up = []
BVF_shelf = []
BVF_off = []
for i in range(0,len(TIME)):
    loc_data = np.where(EDGE[:,:,i]*up_mask*NM_mask == 1)
    BVF_up.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF
     
    loc_data = np.where(EDGE[:,:,i]*sb_mask == 1)
    BVF_shelf.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF 
    
    BVF_off.append(np.nanmean(BVF[:,:,i]*mask_offshelf))
    
## Clean up and calculate index

del BVF
del EDGE

#BVF_up_index = np.array(BVF_up) - np.array(BVF_off)
#BVF_shelf_index = np.array(BVF_shelf) - np.array(BVF_off)

BVF_up_index = np.array(BVF_up)
BVF_shelf_index = np.array(BVF_shelf)

BVF_up_index = np.log10(BVF_up_index)
BVF_shelf_index = np.log10(BVF_shelf_index)
BVF_off_index = np.log10(np.array(BVF_off))

#%%

#######################
# LOAD the WIND File
#######################

f = open('WIND.pckl', 'rb')
WIND_COAST = pickle.load(f)
f.close()

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

###################################
# MY FUNCTIONS
###################################

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
    return ind_event 

def STACK_array(array,num_events,seq):
    STACK = np.empty((num_events.astype(int),len(seq)))
    for i in range(1,num_events.astype(int)+1):
        grep = array[(i*len(seq))-len(seq):i*len(seq)]
        STACK[i-1,:] = grep
    return STACK


# When using this function: input must be a regional data frame, indexed correctly
# Output will be in order:
#    (1) The indexed data frame to the extreme events

def get_my_data(region,seq,Threshold):
    wind_region = region['WIND'].to_numpy() # Covert wind to np array 
    wind_region = wind_region>=Threshold    # Boolean array
    myarray = wind_region.astype(float)    # Binary array for the wind events
    ind_event = search_sequence_numpy(myarray,seq)  # Search for the events
    myevents = myarray[ind_event.astype(int)]    # Index array to get events
    num_events = sum(myevents)/sum(seq) # Number of events matching criteria
    ind_event = my_sort(ind_event,num_events,seq)
    wind_region = region['WIND'].to_numpy()
    wind_plot = wind_region[STACK_array(ind_event,num_events,seq).astype(int)]
    ind_event_region = STACK_array(ind_event,num_events,seq).astype(int)
    return  ind_event_region, wind_plot 

#%%

# For the various regions we are going to find the extreme events for summer
# and winter. From there, we will then index the events against the V data
# To simplify things, we will find all the extreme events for the all the different
# regions and comdine them together as a single array fro summer and winter
# removing replicates. 
# Two Hovmullers can then be designed showing the evolution of the wind-events
# along the shore.

# Define fields 

fields = ['WIND']

##################
# Choose region
##################

region1 = ST_HELENA['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_summer.astype(int),:]

allregions_summer = [region1,region2,region3]

region1 = ST_HELENA['Upwelling'][fields].iloc[ind_winter.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_winter.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_winter.astype(int),:]

allregions_winter = [region1,region2,region3]
################################################################
# Define the thresholds

Threshold_sum = 0.057
Threshold_win = 0.022
seq = np.array([0,0,1,1,1,0,0])  # Two days either side quiescent and then strong between

##################################################################

# Use the functions to loop over the regions and store the necessray data

IND_SUM = []
IND_WIN = []
WIND_SUM = []
WIND_WIN = []
for i in range(0,len(allregions_summer)):
    ind_sum = get_my_data(allregions_summer[i],seq,Threshold_sum)[0]
    ind_win = get_my_data(allregions_winter[i],seq,Threshold_win)[0]
    wind_sum = get_my_data(allregions_summer[i],seq,Threshold_sum)[1]
    wind_win = get_my_data(allregions_winter[i],seq,Threshold_win)[1]
    IND_SUM.append(ind_sum)
    IND_WIN.append(ind_win)
    WIND_SUM.append(wind_sum)
    WIND_WIN.append(wind_win)
    
IND_SUM = np.concatenate(IND_SUM, axis=0)
IND_WIN = np.concatenate(IND_WIN, axis=0)
WIND_SUM = np.concatenate(WIND_SUM, axis=0)
WIND_WIN = np.concatenate(WIND_WIN, axis=0)

IND_SUM,ind_sum = np.unique(IND_SUM, axis=0,return_index=True)  # Only unique events across the regions
IND_WIN,ind_win = np.unique(IND_WIN, axis=0,return_index=True)

WIND_SUM = WIND_SUM[ind_sum.astype(int),:]
WIND_WIN = WIND_WIN[ind_win.astype(int),:]

############################################################################

# From plotting the data as is, the hovmuller does show the wave like motion
# of events with the current reversals
#V_SUM = np.gradient(V_SUM,axis=2)
#V_WIN = np.gradient(V_WIN,axis=2)
#WIND_SUM = np.gradient(WIND_SUM,axis=1)
#WIND_WIN = np.gradient(WIND_WIN,axis=1)

V_SUM = V_COAST[:,ind_summer.astype(int)][:,IND_SUM.astype(int)]
V_WIN = V_COAST[:,ind_winter.astype(int)][:,IND_WIN.astype(int)]

DATE_sum = DATE[ind_summer.astype(int)][IND_SUM.astype(int)]
DATE_win = DATE[ind_winter.astype(int)][IND_WIN.astype(int)]

# Add the stratification
BVF_SUM_up = BVF_up_index[ind_summer.astype(int)][IND_SUM.astype(int)]
BVF_WIN_up = BVF_up_index[ind_winter.astype(int)][IND_WIN.astype(int)]

BVF_SUM_shelf = BVF_shelf_index[ind_summer.astype(int)][IND_SUM.astype(int)]
BVF_WIN_shelf = BVF_shelf_index[ind_winter.astype(int)][IND_WIN.astype(int)]

#############################################################################

# From our arrays: our velocity data is 3-D and our winds are 2-D (being stacked)
# We now do not want to plot every event: but for the mean time we will add a variable
# for which to change to alter this property

plot_n = 1

WIND_SUM = WIND_SUM[::plot_n,:]
WIND_WIN = WIND_WIN[::plot_n,:]

V_SUM = V_SUM[:,::plot_n,:]
V_WIN = V_WIN[:,::plot_n,:]

BVF_SUM_up = BVF_SUM_up[::plot_n]
BVF_WIN_up = BVF_WIN_up[::plot_n]
BVF_SUM_shelf = BVF_SUM_shelf[::plot_n]
BVF_WIN_shelf = BVF_WIN_shelf[::plot_n]

num_sum = WIND_SUM.shape[0]  # Need to remeber the number of events were are plotting
num_win = WIND_WIN.shape[0]  # Need to remeber the number of events were are plotting

# Restructure and flatten the arrays

WIND_SUM = WIND_SUM.flatten()
WIND_WIN = WIND_WIN.flatten()

V_SUM = V_SUM.transpose(0,1,2).reshape(len(lat_CROCO),-1)
V_WIN = V_WIN.transpose(0,1,2).reshape(len(lat_CROCO),-1)

#%%

varplot = V_SUM
windplot = WIND_SUM
BVF_var = BVF_SUM_up.flatten()
num_plot = num_sum
dates = np.arange(0,varplot.shape[1], 1, dtype=int)

# Set up the axes with gridspec

fig = plt.figure(figsize=(12, 12))
grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
main_ax = fig.add_subplot(grid[:-1, 1:])
y_map = fig.add_subplot(grid[:-1, 0], sharey=main_ax)
x_wind = fig.add_subplot(grid[-1, 1:], sharex=main_ax)

# Hovmuller on the main axes
gr = main_ax.pcolor(dates,lat_CROCO,varplot,vmin=-0.3,vmax=0.3,cmap=cmo.balance)
main_ax.contour(dates,lat_CROCO,varplot,levels=[0], linewidths=1.0, linestyles='dashed', colors='black')
gr.set_clim(-0.5, 0.5)
main_ax.tick_params(axis='both', which='major', labelsize=12)

cbar_ax = fig.add_axes([0.92, 0.4, 0.03, 0.45])
cb = fig.colorbar(gr, cax=cbar_ax,pad=0.2)
cb.set_label(label='Velocity [m/s]', size=12)
cb.ax.tick_params(labelsize=11)

# plots on the attached axes
# Wind

x_wind.plot(dates,windplot,'--')
x_wind.set_title('(b)',size=13,loc='left')

for i in range(0,int(num_plot/2)):
    if np.mod(num_plot,2) == 0:
        x_wind.axvspan((i*2*len(seq)),(i*2*len(seq))+len(seq),color='grey', alpha=0.75, lw=0)
    else:
        x_wind.axvspan((i*2*len(seq))+len(seq),(i*2*len(seq))+len(seq)+len(seq),color='grey', alpha=0.75, lw=0)

x_wind.set_xticks(np.arange(3, len(dates)-3, len(seq)))     
x_wind.set_xticklabels(np.arange(1,num_plot+1), fontsize=12)
x_wind.set_xlabel('Wind-events', fontsize = 12)
x_wind.set_ylabel('[N/$m^{2}$]',fontsize=12)
x_wind.set_ylim([-0.05, 0.15])

# twin object for two different y-axis on the sample plot
ax2=x_wind.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(dates, BVF_var,color="red",linestyle = '--')
ax2.set_ylabel("log10 ΔBVF [$s^{-2}$]",fontsize=12)
#ax2.set_ylim([-4, 10])
ax2.tick_params(axis='y', which='major', labelsize=12)

# Map

mymask = h_CROCO>=0
mymask = mymask.astype(float)
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)
y_map.pcolor(X,Y,mymask,cmap='Greys_r')
y_map.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=2.0, linestyles='dashed', colors='Green')
y_map.set_title('(a)',size=13)

y_map.tick_params(axis='both', which='major', labelsize=11)

plt.show()

fig.savefig('Hovmuller_summer_velocity', format='png', dpi=600,bbox_inches='tight')

#%% Create a time-series plot with the wind events for the upwelling and shelf-regions showing the frontal movement

# Define fields 

fields = ['GRAD','DIS']

##################
# Choose region
##################

# Upwelling

region1 = ST_HELENA['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_summer.astype(int),:]

allregions_summer = [region1,region2,region3]

region1 = ST_HELENA['Upwelling'][fields].iloc[ind_winter.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_winter.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_winter.astype(int),:]

allregions_winter = [region1,region2,region3]

GRAD_SUM_up = np.zeros((len(allregions_summer),len(IND_SUM.flatten().astype(int))))
GRAD_WIN_up = np.zeros((len(allregions_winter),len(IND_WIN.flatten().astype(int))))
DIS_SUM_up = np.zeros((len(allregions_summer),len(IND_SUM.flatten().astype(int))))
DIS_WIN_up = np.zeros((len(allregions_winter),len(IND_WIN.flatten().astype(int))))
for i in range(0,len(allregions_summer)):
    # Summer
    grad_sum = allregions_summer[i]['GRAD'].to_numpy()
    grad_sum = grad_sum[IND_SUM.astype(int)]
    dis_sum = allregions_summer[i]['DIS'].to_numpy()
    dis_sum = dis_sum[IND_SUM.astype(int)]
    for j in range(0,len(IND_SUM)):
        tmp = grad_sum[j,:]
        tmp_plot = tmp - np.nanmean(tmp)
        grad_sum[j,:] = tmp_plot
        tmp = dis_sum[j,:]
        tmp_plot = tmp - np.nanmean(tmp)
        dis_sum[j,:] = tmp_plot
    GRAD_SUM_up[i,:] = grad_sum.flatten()
    DIS_SUM_up[i,:] = dis_sum.flatten()
    # Winter
    grad_win = allregions_winter[i]['GRAD'].to_numpy()
    grad_win = grad_win[IND_WIN.astype(int)]
    dis_win = allregions_winter[i]['DIS'].to_numpy()
    dis_win = dis_win[IND_WIN.astype(int)]
    for j in range(0,len(IND_WIN)):
        tmp = grad_win[j,:]
        tmp_plot = tmp - np.nanmean(tmp)
        grad_win[j,:] = tmp_plot
        tmp = dis_win[j,:]
        tmp_plot = tmp - np.nanmean(tmp)
        dis_win[j,:] = tmp_plot
    GRAD_WIN_up[i,:] = grad_win.flatten()
    DIS_WIN_up[i,:] = dis_win.flatten()
    
# Shelf-break

region1 = ST_HELENA['Shelf_break'][fields].iloc[ind_summer.astype(int),:]
region2 = MID_SHELF['Shelf_break'][fields].iloc[ind_summer.astype(int),:]
region3 = NAMAQUA['Shelf_break'][fields].iloc[ind_summer.astype(int),:]

allregions_summer = [region1,region2,region3]

region1 = ST_HELENA['Shelf_break'][fields].iloc[ind_winter.astype(int),:]
region2 = MID_SHELF['Shelf_break'][fields].iloc[ind_winter.astype(int),:]
region3 = NAMAQUA['Shelf_break'][fields].iloc[ind_winter.astype(int),:]

allregions_winter = [region1,region2,region3]

GRAD_SUM_shelf = np.zeros((len(allregions_summer),len(IND_SUM.flatten().astype(int))))
GRAD_WIN_shelf = np.zeros((len(allregions_winter),len(IND_WIN.flatten().astype(int))))
DIS_SUM_shelf = np.zeros((len(allregions_summer),len(IND_SUM.flatten().astype(int))))
DIS_WIN_shelf = np.zeros((len(allregions_winter),len(IND_WIN.flatten().astype(int))))
for i in range(0,len(allregions_summer)):
    # Summer
    grad_sum = allregions_summer[i]['GRAD'].to_numpy()
    grad_sum = grad_sum[IND_SUM.astype(int)]
    dis_sum = allregions_summer[i]['DIS'].to_numpy()
    dis_sum = dis_sum[IND_SUM.astype(int)]
    for j in range(0,len(IND_SUM)):
        tmp = grad_sum[j,:]
        tmp_plot = tmp - np.nanmean(tmp)
        grad_sum[j,:] = tmp_plot
        tmp = dis_sum[j,:]
        tmp_plot = tmp - np.nanmean(tmp)
        dis_sum[j,:] = tmp_plot
    GRAD_SUM_shelf[i,:] = grad_sum.flatten()
    DIS_SUM_shelf[i,:] = dis_sum.flatten()
    # Winter
    grad_win = allregions_winter[i]['GRAD'].to_numpy()
    grad_win = grad_win[IND_WIN.astype(int)]
    dis_win = allregions_winter[i]['DIS'].to_numpy()
    dis_win = dis_win[IND_WIN.astype(int)]
    for j in range(0,len(IND_WIN)):
        tmp = grad_win[j,:]
        tmp_plot = tmp - np.nanmean(tmp)
        grad_win[j,:] = tmp_plot
        tmp = dis_win[j,:]
        tmp_plot = tmp - np.nanmean(tmp)
        dis_win[j,:] = tmp_plot
    GRAD_WIN_shelf[i,:] = grad_win.flatten()
    DIS_WIN_shelf[i,:] = dis_win.flatten()

#######################################################################

# Take the mean of the fronts
# Positions    
DIS_SUM_up = np.nanmean(DIS_SUM_up,axis=0)
DIS_SUM_shelf = np.nanmean(DIS_SUM_shelf,axis=0)
DIS_WIN_up = np.nanmean(DIS_WIN_up,axis=0)
DIS_WIN_shelf = np.nanmean(DIS_WIN_shelf,axis=0)
# Gradient
GRAD_SUM_up = np.nanmean(GRAD_SUM_up,axis=0)
GRAD_SUM_shelf = np.nanmean(GRAD_SUM_shelf,axis=0)
GRAD_WIN_up = np.nanmean(GRAD_WIN_up,axis=0)
GRAD_WIN_shelf = np.nanmean(GRAD_WIN_shelf,axis=0)

######################################################################## 

######################################
# Mean displqcement

mean_DIS_sum_up = []
mean_DIS_sum_shelf = []
for i in range(0,int(int(len(IND_SUM.flatten()))/len(seq))):
    tmp = DIS_SUM_up[(i*len(seq)):(i*len(seq))+len(seq)]
    mean_DIS_sum_up = np.append(mean_DIS_sum_up,np.ptp(tmp))
    tmp = DIS_SUM_shelf[(i*len(seq)):(i*len(seq))+len(seq)]
    mean_DIS_sum_shelf = np.append(mean_DIS_sum_shelf,np.ptp(tmp))
    
mean_DIS_win_up = []
mean_DIS_win_shelf = []
for i in range(0,int(int(len(IND_WIN.flatten()))/len(seq))):
    tmp = DIS_WIN_up[(i*len(seq)):(i*len(seq))+len(seq)]
    mean_DIS_win_up = np.append(mean_DIS_win_up,np.ptp(tmp))
    tmp = DIS_WIN_shelf[(i*len(seq)):(i*len(seq))+len(seq)]
    mean_DIS_win_shelf = np.append(mean_DIS_win_shelf,np.ptp(tmp))
    
mean_DIS_sum_up = np.nanmean(mean_DIS_sum_up)    
mean_DIS_sum_shelf = np.nanmean(mean_DIS_sum_shelf)

mean_DIS_win_up = np.nanmean(mean_DIS_win_up)    
mean_DIS_win_shelf = np.nanmean(mean_DIS_win_shelf)
#########################################################################
# Create a time-series plot

fig,ax = plt.subplots(4,1,figsize=(20,10))
# make a plot of summer
# Distance vs wind
dates = np.arange(0,len(IND_SUM.flatten()), 1, dtype=int)
ax[0].plot(dates,
        WIND_SUM.flatten(),
        color="black", 
        linestyle="-",
        marker = "o")
# set x-axis label
#ax[0].set_xlabel("Event", fontsize = 14)
# set y-axis label
#ax[0].set_ylabel("Wind-stress [N/$m^{2}$]",
#              color="black",
#              fontsize=14)
ax[0].tick_params(axis='y', which='major', labelsize=11)

for i in range(0,int(int(len(IND_SUM))/2)):
    if np.mod(int(len(IND_SUM)),2) == 0:
        ax[0].axvspan((i*2*len(seq)),(i*2*len(seq))+len(seq),color='grey', alpha=0.75, lw=0)
    else:
        ax[0].axvspan((i*2*len(seq))+len(seq),(i*2*len(seq))+len(seq)+len(seq),color='grey', alpha=0.75, lw=0)

# twin object for two different y-axis on the sample plot
ax2=ax[0].twinx()
# make a plot with different y-axis using second axis object
ax2.plot(dates, DIS_SUM_up,color="blue",linestyle="-",marker = "o")
ax2.set_ylabel("Displacement [km]",color="blue",fontsize=14)
ax2.axhspan(0,-20,color='red').set_alpha(0.2)
ax2.axhspan(0,20,color='green').set_alpha(0.2)
ax[0].set_xticks(np.arange(3, len(dates)-3, len(seq)))     
ax[0].set_xticklabels(np.arange(1,len(IND_SUM)+1), fontsize=11)
ax2.set_ylim([-20, 20])
ax2.tick_params(axis='y', which='major', labelsize=12)
ax[0].set_title("(a)",loc='left',fontsize = 14)

ax3 = ax[0].twinx()
ax3.spines["right"].set_position(("axes", 1.05))
ax3.plot(dates, BVF_SUM_up.flatten(),color="red",linestyle='-',marker = "o")
ax3.set_ylabel("log10 ΔBVF [$s^{-2}$]",color='red',fontsize=13)
ax3.yaxis.set_label_position('right')
#gradient vs wind

ax[1].plot(dates,
        WIND_SUM.flatten(),
        color="black", 
        marker="o")
# set x-axis label
#ax[1].set_xlabel("Event", fontsize = 14)
# set y-axis label
#ax[1].set_ylabel("Wind-stress [N/$m^{2}$]",
#              color="black",
#              fontsize=14)
ax[1].tick_params(axis='y', which='major', labelsize=11)
for i in range(0,int(int(len(IND_SUM))/2)):
    if np.mod(int(len(IND_SUM)),2) == 0:
        ax[1].axvspan((i*2*len(seq)),(i*2*len(seq))+len(seq),color='grey', alpha=0.75, lw=0)
    else:
        ax[1].axvspan((i*2*len(seq))+len(seq),(i*2*len(seq))+len(seq)+len(seq),color='grey', alpha=0.75, lw=0)

# twin object for two different y-axis on the sample plot
ax2=ax[1].twinx()
# make a plot with different y-axis using second axis object
ax2.plot(dates, GRAD_SUM_up,color="blue",marker="o")
ax2.set_ylabel("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",color="blue",fontsize=14)
ax2.axhspan(0,-0.1,color='red').set_alpha(0.2)
ax2.axhspan(0,0.1,color='green').set_alpha(0.2)
ax[1].set_xticks(np.arange(3, len(dates)-3, len(seq)))     
ax[1].set_xticklabels(np.arange(1,len(IND_SUM)+1), fontsize=11)
ax2.set_ylim([-0.1, 0.1])
ax2.tick_params(axis='y', which='major', labelsize=12)
ax[1].set_title("(b)",loc='left',fontsize = 14)

ax3 = ax[1].twinx()
ax3.spines["right"].set_position(("axes", 1.07))
ax3.plot(dates, BVF_SUM_up.flatten(),color="red",linestyle='-',marker = "o")
ax3.set_ylabel("log10 ΔBVF [$s^{-2}$]",color='red',fontsize=13)
ax3.yaxis.set_label_position('right')

# make a plot of winter
# Distance vs wind

dates = np.arange(0,len(IND_WIN.flatten()), 1, dtype=int)
ax[2].plot(dates,
        WIND_WIN.flatten(),
        color="black", 
        marker="o")
# set x-axis label
#ax[2].set_xlabel("Event", fontsize = 14)
# set y-axis label
ax[2].set_ylabel("Wind-stress ['[N/$m^{2}$]",
              color="black",
              fontsize=15)
ax[2].yaxis.set_label_coords(-0.04, 1.2)
ax[2].tick_params(axis='y', which='major', labelsize=11)

for i in range(0,int(int(len(IND_WIN))/2)):
    if np.mod(int(len(IND_WIN)),2) == 0:
        ax[2].axvspan((i*2*len(seq)),(i*2*len(seq))+len(seq),color='grey', alpha=0.75, lw=0)
    else:
        ax[2].axvspan((i*2*len(seq))+len(seq),(i*2*len(seq))+len(seq)+len(seq),color='grey', alpha=0.75, lw=0)

# twin object for two different y-axis on the sample plot
ax2=ax[2].twinx()
# make a plot with different y-axis using second axis object
ax2.plot(dates, DIS_WIN_up,color="blue",marker="o")
ax2.set_ylabel("Displacement [km]",color="blue",fontsize=14)
ax2.axhspan(0,-20,color='red').set_alpha(0.2)
ax2.axhspan(0,20,color='green').set_alpha(0.2)

ax[2].set_xticks(np.arange(3, len(dates)-3, len(seq)))     
ax[2].set_xticklabels(np.arange(1,len(IND_WIN)+1), fontsize=11)
ax2.set_ylim([-20, 20])
ax2.tick_params(axis='y', which='major', labelsize=12)
ax[2].set_title("(c)",loc='left',fontsize = 14)

ax3 = ax[2].twinx()
ax3.spines["right"].set_position(("axes", 1.05))
ax3.plot(dates, BVF_WIN_up.flatten(),color="red",linestyle='-',marker = "o")
ax3.set_ylabel("log10 ΔBVF [$s^{-2}$]",color='red',fontsize=13)
ax3.yaxis.set_label_position('right')
# Gradeint vs wind

ax[3].plot(dates,
        WIND_WIN.flatten(),
        color="black", 
        marker="o")
# set x-axis label
ax[3].set_xlabel("Event", fontsize = 15)
# set y-axis label
#ax[3].set_ylabel("Wind-stress ['[N/$m^{2}$]",
#              color="black",
#              fontsize=14)
ax[3].tick_params(axis='y', which='major', labelsize=11)

for i in range(0,int(int(len(IND_WIN))/2)):
    if np.mod(int(len(IND_WIN)),2) == 0:
        ax[3].axvspan((i*2*len(seq)),(i*2*len(seq))+len(seq),color='grey', alpha=0.75, lw=0)
    else:
        ax[3].axvspan((i*2*len(seq))+len(seq),(i*2*len(seq))+len(seq)+len(seq),color='grey', alpha=0.75, lw=0)

# twin object for two different y-axis on the sample plot
ax2=ax[3].twinx()
# make a plot with different y-axis using second axis object
ax2.plot(dates, GRAD_WIN_up,color="blue",marker="o")
ax2.set_ylabel("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",color="blue",fontsize=14)
ax2.axhspan(0,-0.1,color='red').set_alpha(0.2)
ax2.axhspan(0,0.1,color='green').set_alpha(0.2)
ax[3].set_xticks(np.arange(3, len(dates)-3, len(seq)))     
ax[3].set_xticklabels(np.arange(1,len(IND_WIN)+1), fontsize=11)
ax2.set_ylim([-0.1, 0.1])
ax2.tick_params(axis='y', which='major', labelsize=12)
ax[3].set_title("(d)",loc='left',fontsize = 14)

ax3 = ax[3].twinx()
ax3.spines["right"].set_position(("axes", 1.07))
ax3.plot(dates, BVF_WIN_up.flatten(),color="red",linestyle='-',marker = "o")
ax3.set_ylabel("log10 ΔBVF [$s^{-2}$]",color='red',fontsize=13)
ax3.yaxis.set_label_position('right')

plt.show()
fig.savefig('Summer_WINTER_TS_3day_event_up', format='png', dpi=600,bbox_inches='tight')

#%%
 
# Create a time-series plot

fig,ax = plt.subplots(4,1,figsize=(20,10))
# make a plot of summer
# Distance vs wind
dates = np.arange(0,len(IND_SUM.flatten()), 1, dtype=int)
ax[0].plot(dates,
        WIND_SUM.flatten(),
        color="black", 
        marker="o")
# set x-axis label
#ax[0].set_xlabel("Event", fontsize = 14)
# set y-axis label
#ax[0].set_ylabel("Wind-stress [N/$m^{2}$]",
#              color="black",
#              fontsize=14)
ax[0].tick_params(axis='y', which='major', labelsize=11)

for i in range(0,int(int(len(IND_SUM))/2)):
    if np.mod(int(len(IND_SUM)),2) == 0:
        ax[0].axvspan((i*2*len(seq)),(i*2*len(seq))+len(seq),color='grey', alpha=0.75, lw=0)
    else:
        ax[0].axvspan((i*2*len(seq))+len(seq),(i*2*len(seq))+len(seq)+len(seq),color='grey', alpha=0.75, lw=0)

# twin object for two different y-axis on the sample plot
ax2=ax[0].twinx()
# make a plot with different y-axis using second axis object
ax2.plot(dates, DIS_SUM_shelf,color="blue",marker="o")
ax2.set_ylabel("Displacement [km]",color="blue",fontsize=14)
ax2.axhspan(0,-20,color='red').set_alpha(0.2)
ax2.axhspan(0,20,color='green').set_alpha(0.2)
ax[0].set_xticks(np.arange(3, len(dates)-3, len(seq)))     
ax[0].set_xticklabels(np.arange(1,len(IND_SUM)+1), fontsize=11)
ax2.set_ylim([-20, 20])
ax2.tick_params(axis='y', which='major', labelsize=12)
ax[0].set_title("(a)",loc='left',fontsize = 14)

ax3 = ax[0].twinx()
ax3.spines["right"].set_position(("axes", 1.05))
ax3.plot(dates, BVF_SUM_shelf.flatten(),color="red",linestyle='-',marker = "o")
ax3.set_ylabel("log10 ΔBVF [$s^{-2}$]",color='red',fontsize=13)
ax3.yaxis.set_label_position('right')
#gradient vs wind

ax[1].plot(dates,
        WIND_SUM.flatten(),
        color="black", 
        marker="o")
# set x-axis label
#ax[1].set_xlabel("Event", fontsize = 14)
# set y-axis label
#ax[1].set_ylabel("Wind-stress [N/$m^{2}$]",
#              color="black",
#              fontsize=14)
ax[1].tick_params(axis='y', which='major', labelsize=11)
for i in range(0,int(int(len(IND_SUM))/2)):
    if np.mod(int(len(IND_SUM)),2) == 0:
        ax[1].axvspan((i*2*len(seq)),(i*2*len(seq))+len(seq),color='grey', alpha=0.75, lw=0)
    else:
        ax[1].axvspan((i*2*len(seq))+len(seq),(i*2*len(seq))+len(seq)+len(seq),color='grey', alpha=0.75, lw=0)

# twin object for two different y-axis on the sample plot
ax2=ax[1].twinx()
# make a plot with different y-axis using second axis object
ax2.plot(dates, GRAD_SUM_shelf,color="blue",marker="o")
ax2.set_ylabel("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",color="blue",fontsize=14)
ax2.axhspan(0,-0.1,color='red').set_alpha(0.2)
ax2.axhspan(0,0.1,color='green').set_alpha(0.2)
ax[1].set_xticks(np.arange(3, len(dates)-3, len(seq)))     
ax[1].set_xticklabels(np.arange(1,len(IND_SUM)+1), fontsize=11)
ax2.set_ylim([-0.1, 0.1])
ax2.tick_params(axis='y', which='major', labelsize=12)
ax[1].set_title("(b)",loc='left',fontsize = 14)

ax3 = ax[1].twinx()
ax3.spines["right"].set_position(("axes", 1.07))
ax3.plot(dates, BVF_SUM_shelf.flatten(),color="red",linestyle='-',marker = "o")
ax3.set_ylabel("log10 ΔBVF [$s^{-2}$]",color='red',fontsize=13)
ax3.yaxis.set_label_position('right')
# make a plot of winter
# Distance vs wind

dates = np.arange(0,len(IND_WIN.flatten()), 1, dtype=int)
ax[2].plot(dates,
        WIND_WIN.flatten(),
        color="black", 
        marker="o")
# set x-axis label
#ax[2].set_xlabel("Event", fontsize = 14)
# set y-axis label
ax[2].set_ylabel("Wind-stress ['[N/$m^{2}$]",
              color="black",
              fontsize=15)
ax[2].yaxis.set_label_coords(-0.04, 1.2)
ax[2].tick_params(axis='y', which='major', labelsize=11)

for i in range(0,int(int(len(IND_WIN))/2)):
    if np.mod(int(len(IND_WIN)),2) == 0:
        ax[2].axvspan((i*2*len(seq)),(i*2*len(seq))+len(seq),color='grey', alpha=0.75, lw=0)
    else:
        ax[2].axvspan((i*2*len(seq))+len(seq),(i*2*len(seq))+len(seq)+len(seq),color='grey', alpha=0.75, lw=0)

# twin object for two different y-axis on the sample plot
ax2=ax[2].twinx()
# make a plot with different y-axis using second axis object
ax2.plot(dates, DIS_WIN_shelf,color="blue",marker="o")
ax2.set_ylabel("Displacement [km]",color="blue",fontsize=14)
ax2.axhspan(0,-20,color='red').set_alpha(0.2)
ax2.axhspan(0,20,color='green').set_alpha(0.2)

ax[2].set_xticks(np.arange(3, len(dates)-3, len(seq)))     
ax[2].set_xticklabels(np.arange(1,len(IND_WIN)+1), fontsize=11)
ax2.set_ylim([-20, 20])
ax2.tick_params(axis='y', which='major', labelsize=12)
ax[2].set_title("(c)",loc='left',fontsize = 14)

ax3 = ax[2].twinx()
ax3.spines["right"].set_position(("axes", 1.05))
ax3.plot(dates, BVF_WIN_shelf.flatten(),color="red",linestyle='-',marker = "o")
ax3.set_ylabel("log10 ΔBVF [$s^{-2}$]",color='red',fontsize=13)
ax3.yaxis.set_label_position('right')
# Gradeint vs wind

ax[3].plot(dates,
        WIND_WIN.flatten(),
        color="black", 
        marker="o")
# set x-axis label
ax[3].set_xlabel("Event", fontsize = 15)
# set y-axis label
#ax[3].set_ylabel("Wind-stress ['[N/$m^{2}$]",
#              color="black",
#              fontsize=14)
ax[3].tick_params(axis='y', which='major', labelsize=11)

for i in range(0,int(int(len(IND_WIN))/2)):
    if np.mod(int(len(IND_WIN)),2) == 0:
        ax[3].axvspan((i*2*len(seq)),(i*2*len(seq))+len(seq),color='grey', alpha=0.75, lw=0)
    else:
        ax[3].axvspan((i*2*len(seq))+len(seq),(i*2*len(seq))+len(seq)+len(seq),color='grey', alpha=0.75, lw=0)

# twin object for two different y-axis on the sample plot
ax2=ax[3].twinx()
# make a plot with different y-axis using second axis object
ax2.plot(dates, GRAD_WIN_shelf,color="blue",marker="o")
ax2.set_ylabel("Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]",color="blue",fontsize=14)
ax2.axhspan(0,-0.1,color='red').set_alpha(0.2)
ax2.axhspan(0,0.1,color='green').set_alpha(0.2)
ax[3].set_xticks(np.arange(3, len(dates)-3, len(seq)))     
ax[3].set_xticklabels(np.arange(1,len(IND_WIN)+1), fontsize=11)
ax2.set_ylim([-0.1, 0.1])
ax2.tick_params(axis='y', which='major', labelsize=12)
ax[3].set_title("(d)",loc='left',fontsize = 14)

ax3 = ax[3].twinx()
ax3.spines["right"].set_position(("axes", 1.07))
ax3.plot(dates, BVF_WIN_shelf.flatten(),color="red",linestyle='-',marker = "o")
ax3.set_ylabel("log10 ΔBVF [$s^{-2}$]",color='red',fontsize=13)
ax3.yaxis.set_label_position('right')


plt.show()
fig.savefig('Summer_WINTER_TS_3day_event_shelf', format='png', dpi=600,bbox_inches='tight')

#%%
 
# Calculate the average position of the front in summer vs winter

fields = ['GRAD','DIS']

##################
# Choose region
##################

# Upwelling

region1 = ST_HELENA['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_summer.astype(int),:]

allregions_summer = [region1,region2,region3]

region1 = ST_HELENA['Upwelling'][fields].iloc[ind_winter.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_winter.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_winter.astype(int),:]

allregions_winter = [region1,region2,region3]

GRAD_SUM_up = np.zeros((len(allregions_summer),len(IND_SUM.flatten().astype(int))))
GRAD_WIN_up = np.zeros((len(allregions_winter),len(IND_WIN.flatten().astype(int))))
DIS_SUM_up = np.zeros((len(allregions_summer),len(IND_SUM.flatten().astype(int))))
DIS_WIN_up = np.zeros((len(allregions_winter),len(IND_WIN.flatten().astype(int))))
for i in range(0,len(allregions_summer)):
    # Summer
    grad_sum = allregions_summer[i]['GRAD'].to_numpy()
    grad_sum = grad_sum[IND_SUM.astype(int)]
    dis_sum = allregions_summer[i]['DIS'].to_numpy()
    dis_sum = dis_sum[IND_SUM.astype(int)]
    for j in range(0,len(IND_SUM)):
        tmp = grad_sum[j,:]
        tmp_plot = tmp 
        grad_sum[j,:] = tmp_plot
        tmp = dis_sum[j,:]
        tmp_plot = tmp 
        dis_sum[j,:] = tmp_plot
    GRAD_SUM_up[i,:] = grad_sum.flatten()
    DIS_SUM_up[i,:] = dis_sum.flatten()
    # Winter
    grad_win = allregions_winter[i]['GRAD'].to_numpy()
    grad_win = grad_win[IND_WIN.astype(int)]
    dis_win = allregions_winter[i]['DIS'].to_numpy()
    dis_win = dis_win[IND_WIN.astype(int)]
    for j in range(0,len(IND_WIN)):
        tmp = grad_win[j,:]
        tmp_plot = tmp 
        grad_win[j,:] = tmp_plot
        tmp = dis_win[j,:]
        tmp_plot = tmp 
        dis_win[j,:] = tmp_plot
    GRAD_WIN_up[i,:] = grad_win.flatten()
    DIS_WIN_up[i,:] = dis_win.flatten()
    
# Shelf-break

region1 = ST_HELENA['Shelf_break'][fields].iloc[ind_summer.astype(int),:]
region2 = MID_SHELF['Shelf_break'][fields].iloc[ind_summer.astype(int),:]
region3 = NAMAQUA['Shelf_break'][fields].iloc[ind_summer.astype(int),:]

allregions_summer = [region1,region2,region3]

region1 = ST_HELENA['Shelf_break'][fields].iloc[ind_winter.astype(int),:]
region2 = MID_SHELF['Shelf_break'][fields].iloc[ind_winter.astype(int),:]
region3 = NAMAQUA['Shelf_break'][fields].iloc[ind_winter.astype(int),:]

allregions_winter = [region1,region2,region3]

GRAD_SUM_shelf = np.zeros((len(allregions_summer),len(IND_SUM.flatten().astype(int))))
GRAD_WIN_shelf = np.zeros((len(allregions_winter),len(IND_WIN.flatten().astype(int))))
DIS_SUM_shelf = np.zeros((len(allregions_summer),len(IND_SUM.flatten().astype(int))))
DIS_WIN_shelf = np.zeros((len(allregions_winter),len(IND_WIN.flatten().astype(int))))
for i in range(0,len(allregions_summer)):
    # Summer
    grad_sum = allregions_summer[i]['GRAD'].to_numpy()
    grad_sum = grad_sum[IND_SUM.astype(int)]
    dis_sum = allregions_summer[i]['DIS'].to_numpy()
    dis_sum = dis_sum[IND_SUM.astype(int)]
    for j in range(0,len(IND_SUM)):
        tmp = grad_sum[j,:]
        tmp_plot = tmp 
        grad_sum[j,:] = tmp_plot
        tmp = dis_sum[j,:]
        tmp_plot = tmp 
        dis_sum[j,:] = tmp_plot
    GRAD_SUM_shelf[i,:] = grad_sum.flatten()
    DIS_SUM_shelf[i,:] = dis_sum.flatten()
    # Winter
    grad_win = allregions_winter[i]['GRAD'].to_numpy()
    grad_win = grad_win[IND_WIN.astype(int)]
    dis_win = allregions_winter[i]['DIS'].to_numpy()
    dis_win = dis_win[IND_WIN.astype(int)]
    for j in range(0,len(IND_WIN)):
        tmp = grad_win[j,:]
        tmp_plot = tmp 
        grad_win[j,:] = tmp_plot
        tmp = dis_win[j,:]
        tmp_plot = tmp 
        dis_win[j,:] = tmp_plot
    GRAD_WIN_shelf[i,:] = grad_win.flatten()
    DIS_WIN_shelf[i,:] = dis_win.flatten()
    
##
# Mean positions

up_mean_sum = np.nanmean(DIS_SUM_up.flatten())
up_mean_win = np.nanmean(DIS_WIN_up.flatten())

shelf_mean_sum = np.nanmean(DIS_SUM_shelf.flatten())
shelf_mean_win = np.nanmean(DIS_WIN_shelf.flatten())

#%% Look at some of the events

####################
# LOAD in the EDGES
####################

f = open('EDGE.pckl', 'rb')
EDGE = pickle.load(f)
EDGE_SUM = EDGE[:,:,ind_summer.astype(int)]
EDGE_WIN = EDGE[:,:,ind_winter.astype(int)]
del EDGE
f.close()

#####################
# LOAD in GRADIENT
#####################

f = open('GRADIENT.pckl', 'rb')
GRADIENT = pickle.load(f)
GRADIENT_SUM = GRADIENT[:,:,ind_summer.astype(int)]
GRADIENT_WIN = GRADIENT[:,:,ind_winter.astype(int)]
del GRADIENT
f.close()

############################
# LOAD currents
############################

f = open('W.pckl', 'rb')
W = pickle.load(f)
W_SUM = W[:,:,ind_summer.astype(int)]
W_WIN = W[:,:,ind_winter.astype(int)]
del W
f.close()

#%%

fields = ['WIND']

##################
# Choose region
##################

region1 = ST_HELENA['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_summer.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_summer.astype(int),:]

allregions_summer = [region1,region2,region3]

region1 = ST_HELENA['Upwelling'][fields].iloc[ind_winter.astype(int),:]
region2 = MID_SHELF['Upwelling'][fields].iloc[ind_winter.astype(int),:]
region3 = NAMAQUA['Upwelling'][fields].iloc[ind_winter.astype(int),:]

allregions_winter = [region1,region2,region3]

# Re-run this code to reset the matrix structure

IND_SUM = []
IND_WIN = []
WIND_SUM = []
WIND_WIN = []
for i in range(0,len(allregions_summer)):
    ind_sum = get_my_data(allregions_summer[i],seq,Threshold_sum)[0]
    ind_win = get_my_data(allregions_winter[i],seq,Threshold_win)[0]
    wind_sum = get_my_data(allregions_summer[i],seq,Threshold_sum)[1]
    wind_win = get_my_data(allregions_winter[i],seq,Threshold_win)[1]
    IND_SUM.append(ind_sum)
    IND_WIN.append(ind_win)
    WIND_SUM.append(wind_sum)
    WIND_WIN.append(wind_win)
    
IND_SUM = np.concatenate(IND_SUM, axis=0)
IND_WIN = np.concatenate(IND_WIN, axis=0)
WIND_SUM = np.concatenate(WIND_SUM, axis=0)
WIND_WIN = np.concatenate(WIND_WIN, axis=0)

IND_SUM,ind_sum = np.unique(IND_SUM, axis=0,return_index=True)  # Only unique events across the regions
IND_WIN,ind_win = np.unique(IND_WIN, axis=0,return_index=True)

WIND_SUM = WIND_SUM[ind_sum.astype(int),:]
WIND_WIN = WIND_WIN[ind_win.astype(int),:]

############################
# Figure
#############################

# The objective is to create from our wind events a plot showing the evolution
# of an upwelling event for summer and winter
# Will overlay the current velocity as well as the SST gradient along with the 
# gradient edges and the respective wind profile to put the whole story together

############################################
# Summer
############################################

# Strong (20)
# Weak (13)
my_ind = 20

ind_YY = IND_SUM[my_ind,:].flatten()
windplot = WIND_SUM[my_ind,:].flatten()

fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(3,7)
for n in range(0,7):
    # add a new subplot iteratively
    ax1 = fig.add_subplot(gs[0, n])
    map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
                  llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=ax1, resolution='l')
    map.fillcontinents(color='grey')
    map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
    map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
    x,y = map(lon_CROCO,lat_CROCO)
    X, Y = np.meshgrid(x, y)
    im1 = ax1.pcolor(X,Y,GRADIENT_SUM[:,:,ind_YY[n]], cmap=cmo.thermal, vmin = 0, vmax = 0.2)
    ax1.contour(X,Y,EDGE_SUM[:,:,ind_YY[n]], linewidths=0.2, linestyles='solid', colors='white')
    ax1.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='lightgreen')
    map.fillcontinents(color='grey')
    
    ax2 = fig.add_subplot(gs[1, n])
    map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
                  llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=ax2, resolution='l')
    map.fillcontinents(color='grey')
    map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
    map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
    x,y = map(lon_CROCO,lat_CROCO)
    X, Y = np.meshgrid(x, y)
    im2 = ax2.pcolor(X,Y,W_SUM[:,:,ind_YY[n]], cmap=cmo.balance, vmin = -5e-5, vmax = 5e-5)
#    ax2.contour(X,Y,EDGE_SUM[:,:,ind_YY[n]], linewidths=0.2, linestyles='solid', colors='yellow')
    ax2.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='green')
    map.fillcontinents(color='grey')
    
ax3 = fig.add_subplot(gs[2, 2:5])
ax3.plot(np.arange(1,7+1, 1, dtype=int),windplot,'--')
ax3.set_xlabel('Days', fontsize = 12)
ax3.set_ylabel('Wind-stress [N/$m^{2}$]',fontsize=12)
ax3.set_ylim([0, 0.2])

ax4 = ax3.twinx()
ax4.spines["right"].set_position(("axes", 1.07))
ax4.plot(np.arange(1,7+1, 1, dtype=int), BVF_SUM_shelf[my_ind,:],color="red",linestyle='-',marker = "o")
ax4.plot(np.arange(1,7+1, 1, dtype=int),BVF_SUM_up[my_ind,:],color="blue",linestyle='-',marker = "o")
ax4.set_ylabel("log10 ΔBVF [$s^{-2}$]",color='red',fontsize=13)
ax4.yaxis.set_label_position('right')

cbar_ax = fig.add_axes([0.94, 0.69, 0.02, 0.15])
cb = fig.colorbar(im1, cax=cbar_ax,pad=0.2)
cb.set_label(label="Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]", size=12)
cb.ax.tick_params(labelsize=11)

cbar_ax = fig.add_axes([0.94, 0.42, 0.02, 0.15])
cb = fig.colorbar(im2, cax=cbar_ax,pad=0.2)
cb.set_label(label="Vertical velocities [m/s]", size=12)
cb.ax.tick_params(labelsize=11)


plt.show()

fig.savefig('Upwelling_summer_event_strong', format='png', dpi=600,bbox_inches='tight')
#############################
# Winter
##############################

# Strong 22
# Weak 15
my_ind = 22
ind_YY = IND_WIN[my_ind,:].flatten()
windplot = WIND_WIN[my_ind,:].flatten()

fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(3,7)
for n in range(0,7):
    # add a new subplot iteratively
    ax1 = fig.add_subplot(gs[0, n])
    map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
                  llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=ax1, resolution='l')
    map.fillcontinents(color='grey')
    map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
    map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
    x,y = map(lon_CROCO,lat_CROCO)
    X, Y = np.meshgrid(x, y)
    im1 = ax1.pcolor(X,Y,GRADIENT_WIN[:,:,ind_YY[n]], cmap=cmo.thermal, vmin = 0, vmax = 0.2)
    ax1.contour(X,Y,EDGE_WIN[:,:,ind_YY[n]], linewidths=0.2, linestyles='solid', colors='white')
    ax1.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='lightgreen')
    map.fillcontinents(color='grey')
    
    ax2 = fig.add_subplot(gs[1, n])
    map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
                  llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], ax=ax2, resolution='l')
    map.fillcontinents(color='grey')
    map.drawparallels(np.arange(np.round(min(lat_CROCO)),np.round(max(lat_CROCO))+1,3.),labels=[True,True,False,False],dashes=[3,3])
    map.drawmeridians(np.arange(np.round(min(lon_CROCO)),np.round(max(lon_CROCO))+1,1.),labels=[False,False,False,True],dashes=[3,3])
    x,y = map(lon_CROCO,lat_CROCO)
    X, Y = np.meshgrid(x, y)
    im2 = ax2.pcolor(X,Y,W_WIN[:,:,ind_YY[n]], cmap=cmo.balance, vmin = -5e-5, vmax = 5e-5)
#    ax2.contour(X,Y,EDGE_WIN[:,:,ind_YY[n]], linewidths=0.2, linestyles='solid', colors='black')
    ax2.contour(X,Y,h_CROCO,levels=[200,300,500], linewidths=1.0, linestyles='dashed', colors='green')
    map.fillcontinents(color='grey')
    
ax3 = fig.add_subplot(gs[2, 2:5])
ax3.plot(np.arange(1,7+1, 1, dtype=int),windplot,'--')
ax3.set_xlabel('Days', fontsize = 12)
ax3.set_ylabel('Wind-stress [N/$m^{2}$]',fontsize=12)
ax3.set_ylim([0, 0.2])

ax4 = ax3.twinx()
ax4.spines["right"].set_position(("axes", 1.07))
ax4.plot(np.arange(1,7+1, 1, dtype=int), BVF_WIN_shelf[my_ind,:],color="red",linestyle='-',marker = "o")
ax4.plot(np.arange(1,7+1, 1, dtype=int),BVF_WIN_up[my_ind,:],color="blue",linestyle='-',marker = "o")
ax4.set_ylabel("log10 ΔBVF [$s^{-2}$]",color='red',fontsize=13)
ax4.yaxis.set_label_position('right')

cbar_ax = fig.add_axes([0.94, 0.69, 0.02, 0.15])
cb = fig.colorbar(im1, cax=cbar_ax,pad=0.2)
cb.set_label(label="Gradient "+"[" + u"\N{DEGREE SIGN}"+"C"+"/km"+"]", size=12)
cb.ax.tick_params(labelsize=11)

cbar_ax = fig.add_axes([0.94, 0.42, 0.02, 0.15])
cb = fig.colorbar(im2, cax=cbar_ax,pad=0.2)
cb.set_label(label="Vertical velocities [m/s]", size=12)
cb.ax.tick_params(labelsize=11)

plt.show()

fig.savefig('Upwelling_winter_event_strong', format='png', dpi=600,bbox_inches='tight')
#%%