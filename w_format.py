#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 11:05:37 2022

@author: jono
"""
##### This script will just format the .mat files for the vertical velocity
##### into a python variable that can be used 

#%%
import scipy.io
import numpy as np
import pickle

#%% Get dimensional variables

# Load the first file

file = file = 'W_'+ str(1) + '.mat'
mat = scipy.io.loadmat(file)

# Get the dimensions and format

mylat = mat['CROCO_lat']
mylat = mylat.transpose()[:,0]
mylat = mylat[0:-1]

mylon = mat['CROCO_lon']
mylon = mylon.transpose()[0,:]
mylon = mylon[0:-1]

##############
# LOOP
##############

W = np.empty((len(mylat),len(mylon),0))
for i in range(1,4+1):
    file = 'W_'+ str(i) + '.mat'
    mat = scipy.io.loadmat(file)
    W_mat = mat['W_'+str(i)]
    W_mat = W_mat[0:-1,0:-1,:]
    W = np.concatenate((W, W_mat), axis=2)

##############
# SAVE
##############

f = open('W.pckl', 'wb')
pickle.dump(W, f)
f.close()
