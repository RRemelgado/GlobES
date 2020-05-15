# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 09:34:29 2020

@author: rr70wedu
"""

#------------------------------------------------------------------------------------------------------------------#
# 0. import modules and read arguments
#------------------------------------------------------------------------------------------------------------------#

from argparse import ArgumentParser
from rasterio.mask import mask
import rasterio as rt
import pandas as pd
import numpy as np
import fiona as fn

parser = ArgumentParser(description = 'IUCN ecosystem validation (example using one ecosystem and one species)')
parser.add_argument("feature", help = "vector file (e.g. shapefile) with species range map")
parser.add_argument("ecosystem", help = "index")
parser.add_argument("output", help = "output directory")

options = parser.parse_args()
feature = options.feature
ecosystem = options.ecosystem
odir = options.output

#------------------------------------------------------------------------------------------------------------------#
# 1. read reference data
#------------------------------------------------------------------------------------------------------------------#

# read range map geometry
sp = fn.open(feature)
f = [s['geometry'] for s in sp]

# extract subset of ecosystem map
ecoMap = rt.open(ecosystem)
ia = mask(ecoMap, f, crop=True, all_touched=True, nodata=-1, indexes=1)[0]

#------------------------------------------------------------------------------------------------------------------#
# 2. validate
#------------------------------------------------------------------------------------------------------------------#

# validate and extract relavent information
val = (np.max(ia[ia > -1]) > 0).astype('uint8') # check if ecosystem occurs in range map (1 if TRUE, 0 if FALSE)
npx = np.sum(ia > -1) # size of range map in pixels

ia = None
f = None
sp = None

#------------------------------------------------------------------------------------------------------------------#
# 3. export validation results
#------------------------------------------------------------------------------------------------------------------#

# return data frame (with validation)
df = pd.DataFrame({'vote':list(val), 'nr_px':npx})
df.to_csv(odir + 'validation.csv', index=False)
