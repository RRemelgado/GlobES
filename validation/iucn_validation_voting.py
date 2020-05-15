# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 09:34:29 2020

@author: rr70wedu
"""

#------------------------------------------------------------------------------------------------------------------#
# 0. import modules
#------------------------------------------------------------------------------------------------------------------#

from argparse import ArgumentParser
from rasterio.mask import mask
import rasterio as rt
import pandas as pd
import numpy as np
import fiona as fn
import glob2 as g

#-----------------------------------------------------------------------------#
#
#-----------------------------------------------------------------------------#

hab = pd.read_csv(idir + 'IUCN_ecoMap_val-iucn_selectedSpecies_subset.csv', dtype={'taxonID': str, 'cid': str})

fileCode = ('{0:0'+str(len(str(hab.size))) + 'd}')

tid = list(np.unique(hab['taxonID'].values))[index]
ref = rt.open(g.glob(hdir + '*.tif')[0]) # reference array

#==============================================================================================================#
# 4. build function to validat ecossystem map
#==============================================================================================================#

eco = list(np.unique(hab['cid'][hab['taxonID'] == tid]))
n = len(eco)

# access range map
sp = fn.open(rdir + 'IUCN_habRange-' + tid + '_NA_NA.json')
f = [s['geometry'] for s in sp]

val = np.zeros(n, dtype='uint8')
npx = np.zeros(n, dtype='float32')

for c in range(0, n):
    
    # access ecosystem map
    ecoMap = rt.open(hdir + 'IUCN_ecoMapME_5am-' + eco[c] + '_19920101-20180101_5arcMin.vrt')
    
    # validate
    ia = mask(ecoMap, f, crop=True, all_touched=True, nodata=-1, indexes=1)[0]
    val[c] = (np.max(ia[ia > -1]) > 0).astype('uint8')
    npx[c] = np.sum(ia > -1)
    
    ia = None

# return data frame (with validation)
df = pd.DataFrame({'taxonID':tid, 'check':list(val), 'code':eco, 'nr_px':npx, 'nr_eco':n})
df.to_csv(odir + 'ecoVal-' + fileCode.format(index) + '_19920101-20180101_5arcMin.csv', index=False)

print(str(index) + ' done!')
