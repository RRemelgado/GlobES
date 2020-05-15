# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 16:34:17 2020

@author: rr70wedu
"""

#-----------------------------------------------------------------------------#
#
#-----------------------------------------------------------------------------#

import rasterio as rt
import pandas as pd
import numpy as np
import glob2 as g
import os
import re

odir = '/work/remelgad/IUCN_ecoMap_val/iucn/'
rdir = '/work/remelgad/IUCN_ecoMap_val/iucn/IUCN_habRange_2dg/'
hdir = '/work/remelgad/IUCN_ecoMap_val/iucn/IUCN_ecoMapME_2dg/'

#-----------------------------------------------------------------------------#
#
#-----------------------------------------------------------------------------#

hab = pd.read_csv(odir + 'IUCN_ecoMap_val-iucn_selectedSpecies_subset.csv', dtype={'taxonID':str,'cid': str})
files = g.glob(rdir + '*.tif')

#-----------------------------------------------------------------------------#
#
#-----------------------------------------------------------------------------#

# build reference layers
ref = rt.open(files[0])

accuracy = np.zeros(ref.shape, dtype='float32')
constant = np.zeros(ref.shape, dtype='float32')

for f in range(0, len(files)):
    
    tid = re.split('[-|_]', os.path.basename(files[f]))[3]
    
    try:
        
        # read species range map
        rds = rt.open(rdir + 'IUCN_habRange_2dg-' + tid + '_NA_2degree.tif').read(1)
        
        # nr. of ecosystems
        eco = hab[hab['taxonID'] == tid]
        
        # extract accuracy
        h = c = np.zeros(rds.shape, dtype='float32')
        for i in range(0, eco.shape[0]):
            try:
                ids = hdir + 'IUCN_ecoMapME_2dg-' + eco['cid'].values[i] + '_19920101-20180101_2degree.tif'
                h = h + (rt.open(ids).read(1) > 0).astype('uint8') * rds
                c = rds
            except:
                h = h + 0
        
        accuracy = accuracy + (h > 0).astype('uint8')
        constant = constant + c
        
        eco = None
        n = None
        h = None
        rds = None
        
        print(str(f) + ' done')
    
    except:
        
        print(str(f) + ' skipped')

#-----------------------------------------------------------------------------#
#
#-----------------------------------------------------------------------------#

oa0 = np.zeros(accuracy.shape, dtype='float32')
oa0[:] = -1
px = np.where(accuracy > 0)
oa0[px] = accuracy[px] / constant[px]
px = np.where((accuracy == 0) & (constant > 0))
oa0[px] = 0

#-----------------------------------------------------------------------------#
#
#-----------------------------------------------------------------------------#

p = ref.meta.copy()
p.update(dtype='float32', compress='deflate', predict=2, zlevel=9, nodata=-1)
ods = rt.open(odir + 'IUCN_ecoMap_val_2dg-speciesAcc_19920101-20180101_2degree.tif', 'w', **p)
ods.write(oa0, indexes=1)
ods.close()

p.update(dtype='int16', compress='deflate', predict=2, zlevel=9, nodata=-1)
ods = rt.open(odir + 'IUCN_ecoMap_val_2dg-speciesCount_19920101-20180101_2degree.tif', 'w', **p)
ods.write(constant.astype('int16'), indexes=1)
ods.close()
