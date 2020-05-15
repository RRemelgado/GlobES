# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 20:34:15 2020

@author: rr70wedu
"""

#----------------------------------------------------------------------------------------------------------------------------#
# 0 load modules and arguments
#----------------------------------------------------------------------------------------------------------------------------#

from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from argparse import ArgumentParser
import rasterio as rt
import pandas as pd
import numpy as np

# read path to output directory
parser = ArgumentParser(description = 'fill data gaps in the Global Surface Water dataset')
parser.add_argument("variables", help = "input csv with listing input variables")
parser.add_argument("output", help = "path to output directory")

options = parser.parse_args()
variables = pd.read_csv(options.variables)
odir = options.output

#----------------------------------------------------------------------------------------------------------------------------#
# 1. acces input variables
#----------------------------------------------------------------------------------------------------------------------------#

mwater = rt.open(variables['variable'][variables['variable'] == 'missing data'].values[0]) # 300 m
pwater = rt.open(variables['variable'][variables['variable'] == 'permanent water (stack)'].values[0]) # 300 m
swater = rt.open(variables['variable'][variables['variable'] == 'seasonal water (stack)'].values[0]) # 1 km
ptotal = rt.open(variables['variable'][variables['variable'] == 'precipitation (year total stack)'].values[0]) # 1km
sflow = rt.open(variables['variable'][variables['variable'] == 'streamflow (year mean stack)'].values[0]) # 1 km
mnPDSI = rt.open(variables['variable'][variables['variable'] == 'PDSI (min stack)'].values[0]) # 4 km
mxPDSI = rt.open(variables['variable'][variables['variable'] == 'PDSI (max stack)'].values[0]) # 4 km
rivers = rt.open(variables['variable'][variables['variable'] == 'PDSI (max stack)'].values[0]) # 300 m
waterCover = rt.open(variables['variable'][variables['variable'] == 'PDSI (max stack)'].values[0]) # 300 m

#----------------------------------------------------------------------------------------------------------------------------#
# 2. extract relevant metadata
#----------------------------------------------------------------------------------------------------------------------------#

dims1 = mwater.shape # dimensions for 300 m res. data
dims2 = sflow[0].shape # dimensions for 1 km res. data
dims3 = mnPDSI[0].shape # dimensions for 4 km res. data

time = list(range(1, dims1[2]+1)) # time vector (neeed for modelling)

# extract and update 
meta = mwater.meta.copy()
meta.updata(dtype='float32', compress='deflate', predict=2, zlevel=9, bands=1)

#----------------------------------------------------------------------------------------------------------------------------#
# 3. read reference data and predictive variables
#----------------------------------------------------------------------------------------------------------------------------#

ia0 = np.zeros(dims1, dtype='uint8') # missing data
ia1 = np.zeros(dims2, dtype='uint8') # permanent water
ia2 = np.zeros(dims1, dtype='uint8') # seasonal water
sfa = np.zeros(dims2, dtype='float32') # streamflow
d1a = np.zeros(dims2, dtype='float32') # PDSI (min)
d2a = np.zeros(dims3, dtype='float32') # PDSI (max)
tpa = np.zeros(dims3, dtype='float32') # precipitation

for y in range(0, dims2[2]):
    ia0[:,:,y] = mwater.read(y)
    ia1[:,:,y] = swater.read(y)
    ia2[:,:,y] = pwater.read(y)
    sfa[:,:,y] = sflow.read(y)
    d1a[:,:,y] = mnPDSI.read(y)
    d2a[:,:,y] = mxPDSI.read(y)
    tpa[:,:,y] = ptotal.read(y)

#----------------------------------------------------------------------------------------------------------------------------#
# 4. find pixels to interpolate
#----------------------------------------------------------------------------------------------------------------------------#

mmv = np.max(ia0, axis=-1) # maximum ratio of missing data
px = np.where(mmv > 0) # find reference pixels
n = len(px[0]) # number of pixels

mmv = None

#----------------------------------------------------------------------------------------------------------------------------#
# 5. interpolate missing values
#----------------------------------------------------------------------------------------------------------------------------#

# arrays with correlations
r2s = np.zeros((dims1[0],dims1[1]), dtype='float32')
r2p = np.zeros((dims1[0],dims1[1]), dtype='float32')

for p in range(0, n):
    
    #=====================================================================#
    # 5.1 find misisng data in time and read water percent data
    #=====================================================================#
    
    m = ia0[px[0][p],px[1][p],]
    i = np.where(m == 0)
    o = np.where(m > 0)
    
    y1 = ia1[px[0][p],px[1][p],] # permanent water
    y2 = ia2[px[0][p],px[1][p],] # seasonal water
        
    #=====================================================================#
    # 5.2. read predictive variables
    #=====================================================================#
    
    # determine coordinate of current pixel
    xy = mwater.xy(px[0][p],px[1][p])
    
    # find matching pixel for 1km data
    po = sflow[0].index(xy[0],xy[1])
    v1 = sfa[po[0],po[1],] # streamflow
    v2 = tpa[po[0],po[1],] # precipitation
    
    # find matching pixel for 4km data
    po = mnPDSI.index(xy[0],xy[1])
    v3 = d1a[po[0],po[1],] # PDSI (min)
    v4 = d2a[po[0],po[1],] # PDSI (max)
    
    # prepare reference table
    train_x = pd.DataFrame({'y':time[i], 's':v1[i], 'p':v2[i], 'd1':v3[i], 'd2':v4[i]})
    pred_x = pd.DataFrame({'y':time[o], 's':v1[i], 'p':v2[o], 'd1':v3[o], 'd2':v4[o]})
    
    #=====================================================================#
    # 5.3. build predictive models and fill missing data
    #=====================================================================#
    
    # build model and predict permanent water
    rf = RandomForestRegressor().fit(train_x,y1[i]) # build model
    y1[o] = rf.predict(pred_x) # predict missing values
    ia1[px[0][p],px[1][p],] = y1 # update array
    
    # correlate "good" data and corresponding predicted values (info on quality)
    r2s[px[0][p],px[1][p]] = r2_score(y1[i],rf.predict(train_x))
    
    # build model and predict seasonal water
    rf = RandomForestRegressor().fit(train_x,y2[i])
    y2[o] = rf.predict(pred_x)
    ia2[px[0][p],px[1][p],] = y2
    
    r2p[px[0][p],px[1][p]] = r2_score(y2[i],rf.predict(train_x))
    
    # update missing data
    ia0[px[0][p],px[1][p],o] = 0
    
    print(p)

#=========================================================================#
# 6. write layers with per-pixel correlation values
#=========================================================================#

ods = rt.open(odir + 'permanent_r2.tif', 'w', **meta)
ods.write(r2s, indexes=1)
ods.close()
r2s = None

ods = rt.open(odir + 'seasonal_r2.tif', 'w', **meta)
ods.write(r2p, indexes=1)
ods.close()
r2p = None
