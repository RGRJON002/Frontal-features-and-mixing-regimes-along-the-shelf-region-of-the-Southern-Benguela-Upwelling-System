#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 12:19:22 2022

@author: jono
"""

# This script is intended to calculate the Bunt-Vasallis Frequency for the near coastal
# region of the SBUS

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

# Read in a test file

file_base =  '/media/data/CHPC_SBUS_3km/'
file_name = 'croco_avg_Y2005M6.nc'   

nc_file = file_base + file_name
root = Dataset(nc_file)
zeta = root.variables['zeta'][:]  # Free surface
h = root.variables['h'][:]    # Topography

#%%

N = 60;
theta_s    =  5.;
theta_b    =  0.;
hc         = 10.;
vtransform =  2.;
vartype = 'r';

#%%

zeta = zeta[0,:,:]

def csf(sc,theta_s,theta_b):
    import math
    if theta_s>0:
        csrf = np.divide(1-np.cosh(np.multiply(sc,theta_s)),(np.cosh(theta_s)-1))
    else:
        csrf = -sc**2
    if theta_b > 0:
        h = np.divide(math.exp(np.multiply(theta_b,csrf))-1,1-math.exp(-theta_b))
    else:
        h = csrf     
    return h

def zlevels(h,zeta,theta_s,theta_b,hc,N,type):

    M = np.size(h,axis=0)
    L  = np.size(h,axis=1)

    sc_r=np.zeros((N,1))
    Cs_r=np.zeros((N,1))
    sc_w = np.zeros((N+1,1))
    Cs_w = np.zeros((N+1,1))
        
    ds = 1./N
    if type == 'w':
        sc_w[0] = -1.0
        sc_w[N+1] = 0
        Cs_w[0] = -1.0
        Cs_w[N+1] = 0
        
        sc_w[1:N] = ds*(np.arange(1,N)-N)        
        
    sc = ds*(np.arange(1,N+1)-N-0.5)

    
     
    Cs_r = csf(sc,theta_s,theta_b)

    sc_r = sc

    h[h==0] = 1.e-2
    Dcrit = 0.01
    zeta[zeta<(Dcrit-h)] = Dcrit -h[zeta<(Dcrit-h)]

    z = np.zeros((N,M,L))

    cff1 = Cs_r
    sc = sc_r

    h2 = (h+hc)
    cff = hc*sc
    h2inv = np.divide(1,h2)

    for k in range(0, N):
        z0 = cff[k] + (cff1[k]*h)
        z[k,:,:] = np.multiply(z0,np.divide(h,h2)) + np.multiply(zeta,(1+np.multiply(z0,h2inv)))

    dz = np.abs(z[0:N-1,:,:]) - np.abs(z[1:N,:,:]) 
    dz_top = h - np.abs(z[0,:,:])
    dz_top = dz_top.data
    dz = np.append(np.transpose(np.atleast_3d(dz_top),(2,0,1)),dz, axis=0)
    
    return z, dz

z = zlevels(h,zeta,theta_s,theta_b,hc,N)



#%%

# Equation

# BV = sqrt(-g/rho_0 d rho/d z)


# Red in the value of density










