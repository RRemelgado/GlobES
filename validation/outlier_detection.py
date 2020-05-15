# -*- coding: utf-8 -*-
"""
Created on Fri May 15 18:14:46 2020

@author: rr70wedu
"""

#----------------------------------------------------------------------------------------------------------------------------#
# 0. load modules and arguments
#----------------------------------------------------------------------------------------------------------------------------#

from argparse import ArgumentParser
import rasterio as rt
import numpy as np

# read path to output directory
parser = ArgumentParser(description = 'quantify ecossystem change')
parser.add_argument("data", help = "stack with ecosystem time-series")
parser.add_argument("output", help = "path to output directory")

options = parser.parse_args()
target = options.data
odir = options.output

#----------------------------------------------------------------------------------------------------------------------------#
# 1. prepare function to detect outliers
#----------------------------------------------------------------------------------------------------------------------------#

def main(x):
    
    # estimater inter-annual percent changes
    v = [0]
    v = np.array(v + list(np.diff(x) / x[0:(len(x)-1)]))
    
    o = v * 0
    o[x > 0] = 1 # start from the premise that there are no outliers
    
    q2 = np.median(v) # 50% quantile
    q1 = np.median(v[v < q2]) # 25% quantile
    q3 = np.median(v[v > q2]) # 75% quantile
    
    if not np.isnan(q1) and not np.isnan(q3):
        
        iqr = 1.5*(q3-q1) # inter-quantile range
        o[((v < (q1-iqr)) | (v > (q3+iqr))) & (x > 0)] = 2 # assign outlier ID
        
    return(o)

#----------------------------------------------------------------------------------------------------------------------------#
# 2. find outliers
#----------------------------------------------------------------------------------------------------------------------------#

# access ecosystem stack
data = rt.open(target)
dims = data.shape

# read input data
ia = np.zeros(dims, dtype=data.dtype)
for i in range(0, dims[2]):
    ia[:,:,i] = data.read(i)

# output array
oa = np.zeros(dims, dtype=data.dtype)

# identify outliers
for r in range(0, dims[0]):
    for c in range(0, dims[1]):
        v = ia[r,c,] # target vector
        oa[r,c,] = main(v) # confidence classes

#----------------------------------------------------------------------------------------------------------------------------#
# 3. write output
#----------------------------------------------------------------------------------------------------------------------------#

meta = data.meta.copy()
meta.update(driver='GTiff', dtype='uint8', compress='deflate', predict=2, zlevel=9, nodata=None)
ods = rt.open(odir + target.split('.') + '_qa.tif', 'w', **meta)
for i in range(0,dims[2]):
    ods.write(oa[:,:,i], indexes=i)
ods.close()
