#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 10:44:25 2022

@author: jono
"""

import pandas as pd
from sklearn import datasets
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler

iris = datasets.load_iris()

target_names = {
    0:'setosa',
    1:'versicolor',
    2:'virginica'
}

df = pd.DataFrame(iris.data,columns=iris.feature_names)
df['target'] = iris.target
df['target_names'] = df['target'].map(target_names)

sns.countplot(x='target_names',data = df)
plt.title('Iris targets value count')
plt.show()

iris = datasets.load_iris()
X = iris.data
Y = iris.target

# data scaling
x_scaled = StandardScaler().fit_transform(X)

from sklearn.decomposition import PCA
 
pca = PCA(n_components=3)
 
pca_features = pca.fit_transform(x_scaled)
 
print('Shape before PCA: ', x_scaled.shape)
print('Shape after PCA: ', pca_features.shape)
 
pca_df = pd.DataFrame(data=pca_features, columns=['PC1', 'PC2', 'PC3'])

import matplotlib.pyplot as plt 
 
from sklearn.decomposition import PCA
sns.set()
 
# Reduce from 4 to 3 features with PCA
pca = PCA(n_components=3)
 
# Fit and transform data
pca.fit_transform(x_scaled)
 
# Bar plot of explained_variance
plt.bar(range(1,len(pca.explained_variance_)+1),pca.explained_variance_)
 
 
plt.xlabel('PCA Feature')
plt.ylabel('Explained variance')
plt.title('Feature Explained Variance')
plt.show()

import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
 
from sklearn.decomposition import PCA
sns.set()
 
# Reduce from 4 to 3 features with PCA
pca = PCA(n_components=3)
 
# Fit and transform data
reduced_features = pca.fit_transform(x_scaled)
 
# Bar plot of explained_variance
plt.bar(range(1,len(pca.explained_variance_)+1),pca.explained_variance_)
 
plt.plot(range(1,len(pca.explained_variance_ )+1),np.cumsum(pca.explained_variance_),c='red',
    label='Cumulative Explained Variance')
 
plt.legend(loc='upper left')
plt.xlabel('Number of components')
plt.ylabel('Explained variance (eignenvalues)')
plt.title('Scree plot')
 
plt.show()