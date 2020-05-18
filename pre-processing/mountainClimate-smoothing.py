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

# access reference climate map
ids = rt.open('path to original BECK climate map')

# initiate output
p = ids.meta.copy()
p.update(compress='deflate', predict=2, zlevel=9)
ods = rt.open('path to smoothed BECK climate map', 'w', **p)

# access reference mountain polygons
sp = fn.open('path to mountain vector data')

#----------------------------------------------------------------------------------------------------------------------------#
# 2. iterate through each mountain and re-classify boreal/artic/tundra climate when applicable
#----------------------------------------------------------------------------------------------------------------------------#

for s in sp:
    
    #=========================================================================#
    # 2.1. load array subset(s)
    #=========================================================================#
    
    f = [s['geometry']] # extract geometry from feature
    o = mask(ids, f, crop=True, indexes=1, all_touched=True, pad=True)
    t = o[1] # transform (used to determine I/) window)
    o = o[0] # extract array
    i = (o > 0) # valid pixel mask  (no data is 0)
    
    f = None
    
    #=========================================================================#
    # 2.2. classify climates into major groups
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
    # 2.3. infer dominant climate within mountain region
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
        
        px1 = np.where(ct == 5) # locate dry tropical pixel count
        px2 = np.where(ct == 6) # locate moist tropical pixel count
        if ((cc[px1]+cc[px2]) / np.sum(cc)) > 0.3:
            dc = 6
            px1 = None
            px2 = None
        
    else:
        
        dc = 1 # default to arctic/antarctic when no other climate is found
    
    ct = None
    cc = None
    
    #=========================================================================#
    # 5.4. update climate map
    #=========================================================================#
    
    # determine I/O window
    px = ids.index(t.c, t.f)
    ad = i.shape
    w = rt.windows.Window(px[1], px[0], ad[1], ad[0])
    
    px = None
    ad = None
    
    # count pixel contributions of each climate type
    cc = np.unique(r, return_counts=True)
    ct = cc[0] # climate types
    cc = cc[1] # pixel count per climate type
    
    r = None
    
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
    ca[i & (np.isin(ca, brc) | np.isin(ca, [29,30]))] = ac
    ids.write(ca, window=w, indexes=1)
    
    ca = None
    i = None

# close output array (needed to finish writting)
ods.close()
ids.close()
