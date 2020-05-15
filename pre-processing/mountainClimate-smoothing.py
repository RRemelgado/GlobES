# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 20:49:22 2019

@author: rr70wedu
"""

#----------------------------------------------------------------------------------------------------------------------------#
# 0. load required modules
#----------------------------------------------------------------------------------------------------------------------------#

from rasterio.mask import mask
import rasterio as rt
import fiona as fn
import numpy as np

# climate codes in the BECK classification map
brc = [19, 20, 23, 24, 27, 28] # boreal
tpc = [7, 14, 15, 16, 17, 18, 21, 22, 25, 26] # temperate
tdc = [3, 6, 11, 12, 13] # tropical (dry)
twc = [1, 2] # tropical (wet)
mdc = [8, 9, 10] # mediterranean
dch = 4 # desert (hot)
dcc = 5 # desert (cold)

#----------------------------------------------------------------------------------------------------------------------------#
# 1. load climate map
#----------------------------------------------------------------------------------------------------------------------------#

i = '/data/idiv_meyer/temp/GMBA/' # GMBA data directory
ids = rt.open(i + 'BECK_climateMap-presentClass_19860101-20160101_30arcSec.tif', 'r+') # read climate map

p = ids.meta.copy()
p.update(driver='GTiff', compress='deflate', predict=2, zlevel=9, nodata=0)

#----------------------------------------------------------------------------------------------------------------------------#
# 1. load climate map
#----------------------------------------------------------------------------------------------------------------------------#

i = '/data/idiv_meyer/00_data/processed/utilityLayers/'
latMap = rt.open(i + 'utilityLayers-latitude_NA_10arcSec.tif')

i = '/data/idiv_meyer/temp/GMBA/' # GMBA data directory
sp = fn.open(i + 'GMBA-merge.shp') # read mountain shapefile

#----------------------------------------------------------------------------------------------------------------------------#
# 1. load climate map
#----------------------------------------------------------------------------------------------------------------------------#

# map of mountains with climate info
ods = rt.open(i + 'tropicalMountains.tif', 'w+', **p)

#----------------------------------------------------------------------------------------------------------------------------#
# 1. load climate map
#----------------------------------------------------------------------------------------------------------------------------#

for s in sp:
    
    #=========================================================================#
    # 5.1. load array subset(s)
    #=========================================================================#
    
    f = [s['geometry']] # extract geometry from feature
    o = mask(ids, f, crop=True, indexes=1, all_touched=True, pad=True)
    t = o[1] # transform (used to determine I/) window)
    o = o[0] # extract array
    i = (o > 0) # derive valid pixel mask 
    lm = mask(latMap, f, crop=True, all_touched=True, pad=True, indexes=1)[0]
    lm = np.max(lm) # maximum latitude
    
    #=========================================================================#
    # 5.2. classify climates into major groups
    #=========================================================================#
    
    px = np.isin(o, [0, 29, 30], invert=True)
    r = o[px] # subet array to non-articl/boreal climates (reference)
    a = o[px] # subet array to non-articl/boreal climates (to classify)
    a[np.isin(r, brc)] = 2 # boreal
    a[np.isin(r, mdc)] = 3 # mediterranean
    a[np.isin(r, tpc)] = 4 # temperate
    a[np.isin(r, tdc)] = 5 # tropical (dry)
    a[np.isin(r, twc)] = 6 # tropical (wet)
    a[np.isin(r, dch)] = 7 # desert (hot)
    a[np.isin(r, dcc)] = 8 # desert (cold)
    
    m = None
    o = None
    
    #=========================================================================#
    # 5.3. infer dominant climate within mountain region
    #=========================================================================#
    
    # count climate contributions in terms or nr. of pixels
    dc = np.unique(a, return_counts=True)
    ct = dc[0] # climate types
    cc = dc[1] # climate types count
    
    dc = None
    a = None
    
    # find dominant climate
    if len(ct > 0):
        
        dc = ct[sorted(range(len(ct)), key=cc.__getitem__)[len(ct)-1]]
        
        # are tropical climates present (priority class)
        if (dc != 6) & np.isin(ct, 6).any():
            
            px1 = np.where(ct == 5) # locate dry tropical pixel count
            px2 = np.where(ct == 6) # locate moist tropical pixel count
            if ((cc[px1]+cc[px2]) / np.sum(cc)) > 0.3:
                dc = 6
                px1 = None
                px2 = None
        
        # accounts for the atacama region
        if (dc == 8) & (lm <= 0):
            dc = 6
        
    else:
        
        dc = 1 # default to arctic/antarctic when no other climate is found
    
    ct = None
    cc = None
    
    #=========================================================================#
    # 5.4. update climate-mountain map
    #=========================================================================#
    
    # determine I/O window
    px = ids.index(t.c, t.f)
    ad = i.shape
    w = rt.windows.Window(px[1], px[0], ad[1], ad[0])
    
    px = None
    ad = None
    
    # write array to mountain mask (if tropical moist or cold desert)
    #if (dc == 5) | (dc == 6):
    
    o = ods.read(1, window=w)
    o[(o == 0) & i] = dc
    ods.write(o, window=w, indexes=1)
    
    #=========================================================================#
    # 5.5. update climate map
    #=========================================================================#
    
    # count climate contributions in terms or nr. of pixels
    cc = np.unique(r, return_counts=True)
    ct = cc[0] # climate types
    cc = cc[1] # climate types count
    
    r = None
    
    # "fix" cold desert pixels related to mountains
    if (dc == 7):
        if (np.min(ct) == 4) & (np.min(ct) == 5):
            ca = ids.read(1, window=w)
            ca[i & (ca == 5)] = 4
            ca[i & np.isin(ca, brc)] = 4
            ids.write(ca, window=w, indexes=1)
    
    # "fix" boreal pixels related to mountains
    if (dc == 8):
        if (np.min(ct) == 4) & (np.min(ct) == 5):
            ca = ids.read(1, window=w)
            ca[i & np.isin(ca, brc)] = 5
            ids.write(ca, window=w, indexes=1)
    
    # "fix" tundra/artic climates in mountain areas
    if (dc > 1) & (dc < 7):
        
        # read climate map subset
        ca = ids.read(1, window=w)
        
        # update climate data
        if dc == 2: px = np.isin(ct, brc)
        if dc == 3: px = np.isin(ct, mdc)
        if dc == 4: px = np.isin(ct, tpc)
        if dc == 5: px = np.isin(ct, tdc)
        if dc == 6: px = np.isin(ct, twc)
        
        # determine dominant sub-climate type
        ac = ct[px][sorted(range(len(ct[px])), key=cc.__getitem__)[len(ct[px])-1]]
        
      
    # write data for target mountain
    ids.write(ca, window=w, indexes=1)
    ca = None
    
    f = None
    i = None

# close output array (needed to finish writting)
ods.close()
ids.close()
