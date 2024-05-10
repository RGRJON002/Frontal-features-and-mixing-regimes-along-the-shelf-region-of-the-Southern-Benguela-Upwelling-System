#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 10:11:55 2022

@author: jono
"""

import scipy.io
import numpy as np
import pickle

#%% Get dimensional variables

# Load the first file

file = file = 'BVF_'+ str(1) + '.mat'
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

BVF = np.empty((len(mylat),len(mylon),0))
for i in range(1,4+1):
    file = 'BVF_'+ str(i) + '.mat'
    mat = scipy.io.loadmat(file)
    BVF_mat = mat['BVF_'+str(i)]
    BVF_mat = np.transpose(BVF_mat,(1,0,2)) 
    BVF_mat = BVF_mat[0:-1,0:-1,:]
    BVF = np.concatenate((BVF, BVF_mat), axis=2)

##############
# SAVE
##############

f = open('BVF.pckl', 'wb')
pickle.dump(BVF, f)
f.close()
