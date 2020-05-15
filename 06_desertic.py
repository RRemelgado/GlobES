# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 09:33:13 2019

@author: rr70wedu
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 18:29:01 2019

@author: rr70wedu
"""

#----------------------------------------------------------------------------------------------------------------------------#
# 0 load modules and arguments
#----------------------------------------------------------------------------------------------------------------------------#

from argparse import ArgumentParser
from scipy import ndimage as nd
import rasterio as rt
import pandas as pd
import numpy as np

parser = ArgumentParser(description = 'Classification of desertic ecosystem types')
parser.add_argument("variables", help = "input csv with listing input variables")
parser.add_argument("reference", help = "potential ecosystem layer")
parser.add_argument("output", help = "path to output directory")

options = parser.parse_args()
reference = options.reference
variables = pd.read_csv(options.variables)
odir = options.output

factor = 3 # aggregation factor for total area estimation

# create structure (3x3 window) used for focal operations
s = [[1,1,1], [1,1,1], [1,1,1]]

#----------------------------------------------------------------------------------------------------------------------------#
# 1. load data for input variables
#----------------------------------------------------------------------------------------------------------------------------#

climate_artic = rt.open(variables['variable'][variables['variable'] == 'Artic climate'].values[0]).read(1)
climate_tundra = rt.open(variables['variable'][variables['variable'] == 'Tundra climate'].values[0]).read(1)
snowCover = rt.open(variables['variable'][variables['variable'] == 'snow cover'].values[0]).read(1)
snowRecurrence = rt.open(variables['variable'][variables['variable'] == 'snow cover recurrence'].values[0]).read(1)
soilDepth = rt.open(variables['variable'][variables['variable'] == 'soil depth'].values[0]).read(1)
climate_hdesert = rt.open(variables['variable'][variables['variable'] == 'hot desert'].values[0]).read(1)
climate_cdesert = rt.open(variables['variable'][variables['variable'] == 'cold desert'].values[0]).read(1)
climate_hsteppe = rt.open(variables['variable'][variables['variable'] == 'hot steppe climate'].values[0]).read(1)
climate_csteppe = rt.open(variables['variable'][variables['variable'] == 'cold steppe climate'].values[0]).read(1)
oceanDistance = rt.open(variables['variable'][variables['variable'] == 'ocean distance'].values[0]).read(1)
pWater = rt.open(variables['variable'][variables['variable'] == 'permanent water'].values[0]).read(1)

# make combined desert climate mask
climate_desert = ((climate_hdesert > 0) | (climate_hdesert > 0) | (climate_hsteppe > 0) | (climate_csteppe > 0)).astype('uint8')

climate_hdesert = None
climate_cdesert = None
climate_hsteppe = None
climate_csteppe = None

#----------------------------------------------------------------------------------------------------------------------------#
# 2. read reference data and metadata
#----------------------------------------------------------------------------------------------------------------------------#

# extract/update metadata (needed to write output)
potentialMap = rt.open(reference)
meta = potentialMap.meta.copy()
nodata = potentialMap.nodata # no data value
nc = potentialMap.height # number of rows (original)
onc = nc/factor # output number of rows (output)
nr = potentialMap.width # number of columns
onr = nr/factor # output number of columns (output)
meta.update(driver='GTiff', dtype='float32', compress='deflate', predict=2, zlevel=9, nc=onc, nr=onr, nodata=255)

# read reference map data
potentialMap = potentialMap.read(1)

#----------------------------------------------------------------------------------------------------------------------------#
# 3. read and update data on pixel areas
#----------------------------------------------------------------------------------------------------------------------------#

pixelArea = rt.open(variables['variable'][variables['variable'] == 'pixel area'].values[0]).read(1)

# remove proportion of area occupied by permanent water (e.g. rivers, lakes)
pixelArea = pixelArea - (pixelArea*pWater)

# aggregate pixel area to 1km
pixelArea = pixelArea.reshape(onr, factor, onc, factor).sum(axis=(1, 3))

#----------------------------------------------------------------------------------------------------------------------------#
# 4. classify ecosystem types at 300m (1)
# NOTE: focuses on rocky areas (driven by soil depth) and cold ecosystems (driven by climate and snow recurrence)
#----------------------------------------------------------------------------------------------------------------------------#

# make water mask
oa = np.zeros(potentialMap.shape, dtype='uint8')

# inland rocky areas
px = np.where((potentialMap > 0) & (soilDepth < 1))
oa[px] = 1

# cold desert driven by climate
px = np.where((oa == 0) & (potentialMap > 0) & ((climate_artic > 0) | (climate_tundra > 0)))
oa[px] = 4

# cold desert driven by snow
px = np.where((oa == 0) & (potentialMap > 0) & ((climate_desert > 0) & (snowRecurrence >= 0.5)) | (snowCover > 0))
oa[px] = 4

#----------------------------------------------------------------------------------------------------------------------------#
# 5. find individual pixel regions to classify and derive region statistics
#----------------------------------------------------------------------------------------------------------------------------#

# label connected pixel regions
larr, nf = nd.measurements.label(((potentialMap > 0) & (oa == 0)).astype('uint8'), s)
region_id = np.array(range(1, nf+1)) # unique region ID's

# count of pixels per region
pixel_count = np.unique(larr[larr > 0], return_counts=True)[1]

# auxiliary variables
coastal = nd.sum((oceanDistance < 333), labels=larr, index=region_id) / pixel_count # % of pixels within 100 km from the ocean
hot = nd.sum(climate_hdesert, labels=larr, index=region_id) / pixel_count # proportion of cold desert climate pixels
cold = nd.sum(climate_cdesert, labels=larr, index=region_id) / pixel_count # proportion of cold desert climate pixels

#----------------------------------------------------------------------------------------------------------------------------#
# 6. lassify ecosystem types at 300m (1)
# NOTE: object-based classification
#----------------------------------------------------------------------------------------------------------------------------#

# cold (coastal) deserts
suitable_region = region_id[(coastal >= 0.5)] # relevant regions
oa[(oa == 0) & np.isin(larr, list(suitable_region))] = 4

# hot desert
suitable_region = region_id[(hot >= 0.5)] # relevant regions
oa[(oa == 0) & np.isin(larr, list(suitable_region))] = 2

# temperate desert
suitable_region = region_id[(cold >= 0.5)] # relevant regions
oa[(oa == 0) & np.isin(larr, list(suitable_region))] = 3

#----------------------------------------------------------------------------------------------------------------------------#
# 7. derive and write 1km ecosystem maps
#----------------------------------------------------------------------------------------------------------------------------#

# rocky areas
oa = (oa == 1).astype('uint8') * potentialMap # find ecosystem pixels
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3)) # aggregate selected pixels
oa = oa * pixelArea # extract area prortion (note: proportion is based on "selected pixels / "non-nodata potential pixels")
ods = rt.open(odir + 'GlobES-0600_1km.tif', 'w', **meta) # initiate area map file
ods.write(oa, indexes=1) # write ecosystem array
ods.close() # close conenction

# hot deserts
oa = (oa == 2).astype('uint8') * potentialMap
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea
ods = rt.open(odir + 'GlobES-0801_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# temperate deserts
oa = (oa == 3).astype('uint8') * potentialMap
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea
ods = rt.open(odir + 'GlobES-0802_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# cold deserts
oa = (oa == 4).astype('uint8') * potentialMap
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea
ods = rt.open(odir + 'GlobES-0803_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()
