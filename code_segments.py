#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 10:33:56 2022

@author: jono
"""

############################
# Old segements of code
############################

#%%%%%%%%%%%%#############################
# Lag correlation of summer wind vs grad
##########################################

def df_derived_by_shift(df,lag=0,NON_DER=[]):
    df = df.copy()
    if not lag:
        return df
    cols ={}
    for i in range(1,lag+1):
        for x in list(df.columns):
            if x not in NON_DER:
                if not x in cols:
                    cols[x] = ['{}_{}'.format(x, i)]
                else:
                    cols[x].append('{}_{}'.format(x, i))
    for k,v in cols.items():
        columns = v
        dfn = pd.DataFrame(data=None, columns=columns, index=df.index)    
        i = 1
        for c in columns:
            dfn[c] = df[k].shift(periods=i)
            i+=1
        df = pd.concat([df, dfn], axis=1)
        df = df.reindex(df.index)
    return df


####################
# Add Date variabe
####################

fields = ['WIND','GRAD'] # mdct is datetime
x = ST_HELENA['Shelf_break'][fields]
x = x.iloc[ind_summer.astype(int)]
NON_DER= []
df_new = df_derived_by_shift(x, 6,NON_DER)
df_new = df_new.dropna()

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
colormap = plt.cm.RdBu
plt.figure(figsize=(15,10))
plt.title(u'St Helena Shelf-break Region 6-Day Lag Summer', y=1.05, size=16)

mask = np.zeros_like(df_new.corr())
mask[np.triu_indices_from(mask)] = True

svm = sns.heatmap(df_new.corr(), mask=mask, linewidths=0.1,vmax=1.0, 
            square=True, cmap=colormap, linecolor='white', annot=True)
 
#fig.savefig('LAG_CORR_Summer_SH_Shelf', format='png', dpi=600)

#%%

# Use a threshold to count extreme wind events per month, per region, per year to see the interannual variability
# From there, how do fronts behave in extreme wind events relative to the mean state and what drives these differences?
# Give an exact example of a particular date

# Looking at the distribution, will use 0.05 N/m^2

# Start by just plotting the wind

fig = plt.figure(figsize=(12,12))

labels = ['St Helena', 'Mid-Shelf', 'Namaqua'] 

plt.title(u'Alongshore Wind-stress [Daily]', y=1.05, size=16)
for frame in [ST_HELENA['Shelf_break'], MID_SHELF['Shelf_break'], NAMAQUA['Shelf_break']]:
    plt.plot(frame.index, frame['WIND'], linestyle='dashed')
plt.hlines(y=0.057, xmin=frame.index[0], xmax=frame.index[-1], linewidth=2, color='r')
plt.hlines(y=0.022, xmin=frame.index[0], xmax=frame.index[-1], linewidth=2, color='b')
plt.xlabel('Date', size=14)
plt.ylabel('Wind-stress [N/m^2]', size=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(labels,borderpad=1)
plt.grid()
plt.show()

#fig.savefig('WIND_DAILY', format='png', dpi=600)
plt.close()

#%%

#%%

# As a proxy for the front deepening or shallowing, we will use the BVF @ 30 m
# Need to extract the data for the BVF as done for the other variables for the frontal diagnostics

#######################
# Get the BVF
#######################

# Edge map 

f = open('EDGE.pckl', 'rb')
EDGE = pickle.load(f)
f.close()

# BVF map

f = open('BVF.pckl', 'rb')
BVF = pickle.load(f)
f.close()

########################################

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
        llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)

# Shelf-break region
BVF_shelf = []
for i in range(0,len(TIME)): 
    loc_data = np.where(EDGE[:,:,i]*sb_mask*NM_mask == 1)
    BVF_shelf.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF 
    
###########################################

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='cyl',llcrnrlon=lon_CROCO[0], urcrnrlon=lon_CROCO[-1], 
        llcrnrlat=lat_CROCO[0], urcrnrlat=lat_CROCO[-1], resolution='l')
map.fillcontinents(color='grey')
x,y = map(lon_CROCO,lat_CROCO)
X, Y = np.meshgrid(x, y)

# Upwelling region
BVF_up = []
for i in range(0,len(TIME)): 
    loc_data = np.where(EDGE[:,:,i]*up_mask*NM_mask == 1)
    BVF_up.append(np.nanmean(BVF[loc_data[0],loc_data[1],i]))  # BVF     

#############################################

#ST_HELENA['Shelf_break']['BVF'] = BVF_shelf 
#ST_HELENA['Upwelling']['BVF'] = BVF_up

#MID_SHELF['Shelf_break']['BVF'] = BVF_shelf
#MID_SHELF['Upwelling']['BVF'] = BVF_up

#NAMAQUA['Shelf_break']['BVF'] = BVF_shelf 
#NAMAQUA['Upwelling']['BVF'] = BVF_up

#%%

#%%
###################################################################
# Evenets above threshold
###################################################################

import matplotlib.dates as mdate
import datetime as dt
from matplotlib.dates import DateFormatter

# Number of high events per season month
# fields = ['WIND','W','DIS','LEN','GRAD','NUM','BVF']

fields = ['WIND','W','DIS','LEN','GRAD','NUM']

# Summer

Threshold = 0.057

gr1 = ST_HELENA['Shelf_break'][fields].iloc[ind_summer.astype(int),:]
gr1 = gr1[fields][(gr1['WIND'] > Threshold)]
gr1 = gr1.groupby([gr1.index.month]).mean() # CLIM

gr1 = gr1.groupby(gr1.index.strftime('%Y-%m')).mean()

gr2 = MID_SHELF['Shelf_break']['WIND'][(MID_SHELF['Shelf_break']['WIND'] > Threshold)]
gr2 = gr2.groupby(gr2.index.strftime('%Y-%m')).count()

gr3 = NAMAQUA['Shelf_break']['WIND'][(NAMAQUA['Shelf_break']['WIND'] > Threshold)]
gr3 = gr3.groupby(gr3.index.strftime('%Y-%m')).count()

gr = pd.concat([gr1,gr2,gr3], axis=1, sort = True)

fig = plt.figure(figsize=(15,12))

labels = ['St Helena', 'Mid-Shelf', 'Namaqua'] 

#plt.title(u'Alongshore Wind-stress [Monthly]', y=1.05, size=16)

for i in range(0,3):
    plt.plot(gr.index, gr.iloc[:,i])
plt.xlabel('Date', size=14)
plt.ylabel('Number of extreme wind events', size=14)
plt.xticks(fontsize=12)
plt.xticks(np.arange(0, len(gr.index)+1, 11))
plt.yticks(fontsize=12)
plt.legend(labels,borderpad=1)
plt.grid()
plt.show()

#fig.savefig('WIND_EVENTS_MONTHLY', format='png', dpi=600)
plt.close()
