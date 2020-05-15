# -*- coding: utf-8 -*-
"""
Created on Thu May 14 16:52:00 2020

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

parser = ArgumentParser(description = 'Classification of open-water \
                        ecosystem types (NOTE: unlike in GlobES, which uses historical \
                        data on water extent to classify water bodies, this script \
                        ignores the historical aspect for the sake of simplicity)')
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

oceanDistance = rt.open(variables['variable'][variables['variable'] == 'ocean distance'].values[0]).read(1)
coralDistance = rt.open(variables['variable'][variables['variable'] == 'coral distance'].values[0]).read(1)
riverDistance = rt.open(variables['variable'][variables['variable'] == 'river distance'].values[0]).read(1)
izDistance = rt.open(variables['variable'][variables['variable'] == 'intertidal zone distance'].values[0]).read(1)
estuaries = rt.open(variables['variable'][variables['variable'] == 'estuaries'].values[0]).read(1) # includes tidal rivers
deltas = rt.open(variables['variable'][variables['variable'] == 'deltas'].values[0]).read(1)
inlandRivers = rt.open(variables['variable'][variables['variable'] == 'inland rivers'].values[0]).read(1)
fwLakes = rt.open(variables['variable'][variables['variable'] == 'freshwater lakes'].values[0]).read(1)
swLakes = rt.open(variables['variable'][variables['variable'] == 'saline lakes'].values[0]).read(1)
bareLand = rt.open(variables['variable'][variables['variable'] == 'bare land and sparse vegetation cover'].values[0]).read(1)
climate_hdesert = rt.open(variables['variable'][variables['variable'] == 'hot desert climate'].values[0]).read(1)
climate_cdesert = rt.open(variables['variable'][variables['variable'] == 'cold desert climate'].values[0]).read(1)
climate_hsteppe = rt.open(variables['variable'][variables['variable'] == 'hot steppe climate'].values[0]).read(1)
climate_csteppe = rt.open(variables['variable'][variables['variable'] == 'cold steppe climate'].values[0]).read(1)
sWater = rt.open(variables['variable'][variables['variable'] == 'seasonal water'].values[0]).read(1)
pWater = rt.open(variables['variable'][variables['variable'] == 'permanent water'].values[0]).read(1)

# make combined desert climate mask
climate_desert = ((climate_hdesert > 0) | (climate_cdesert > 0) | (climate_hsteppe > 0) | (climate_csteppe > 0)).astype('uint8')

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

pixelArea_hr = rt.open(variables['variable'][variables['variable'] == 'pixel area'].values[0]).read(1)

# aggregate pixel area to 1km
pixelArea_agg = pixelArea_hr.reshape(onr, factor, onc, factor).sum(axis=(1, 3))

#----------------------------------------------------------------------------------------------------------------------------#
# 4. classify ecosystem types at 300m, aggregate to 1km and write output (focus on rivers)
#----------------------------------------------------------------------------------------------------------------------------#

# permanent rivers
oa = potentialMap * inlandRivers * pWater # find ecosystem pixels
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3)) # aggregate selected pixels
ods = rt.open(odir + 'GlobES-0501_1km.tif', 'w', **meta) # initiate area map file
ods.write(oa, indexes=1) # write ecosystem array
ods.close() # close conenction

# seasonal rivers
oa = potentialMap * inlandRivers * sWater
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
ods = rt.open(odir + 'GlobES-0502_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# filter classified pixels
potentialMap[inlandRivers > 0] = 0

inlandRivers = None

#----------------------------------------------------------------------------------------------------------------------------#
# 5. find individual pixel regions to classify and derive region statistics
#----------------------------------------------------------------------------------------------------------------------------#

# label connected pixel regions
larr, nf = nd.measurements.label((potentialMap > 0).astype('uint8'), s)
region_id = np.array(range(1, nf+1)) # unique region ID's

# count of pixels per region
pixel_count = np.unique(larr[larr > 0], return_counts=True)[1]

# check for basic statistic applicable to all regions
min_od = nd.minimum(oceanDistance, labels=larr, index=region_id) # min. ocean distance
max_od = nd.maximum(oceanDistance, labels=larr, index=region_id) # max. ocean distance
min_izd = nd.minimum(izDistance, labels=larr, index=region_id) # min. intertidal zone distance
min_cd = nd.minimum(coralDistance, labels=larr, index=region_id) # min. coral distance
min_rd = nd.minimum(riverDistance, labels=larr, index=region_id) # min. river distance
region_area = nd.sum(pixelArea_hr * ((sWater+pWater)*0.01), labels=larr, index=region_id) # region area
bare_inside = nd.sum(bareLand, labels=larr, index=region_id) / pixel_count # percent of pixels covered by bare land
desertic = nd.sum(climate_desert, labels=larr, index=region_id) / pixel_count # percent of pixels in desertic climates

oceanDistance = None
coralDistance = None
riverDistance = None
izDistance = None
desertClimate = None

# check for the presence of known water bodies
fwl = nd.maximum(fwLakes, labels=larr, index=region_id) # freshwater lake occurrence
swl = nd.maximum(swLakes, labels=larr, index=region_id) # saline lake occurrence
est = nd.maximum(estuaries, labels=larr, index=region_id) # estuary occurrence
dlt = nd.maximum(deltas, labels=larr, index=region_id) # deltas occurrence

estuaries = None
corals = None
deltas = None
fwLakes = None
swLakes = None

#----------------------------------------------------------------------------------------------------------------------------#
# 6. classify known water bodies
#----------------------------------------------------------------------------------------------------------------------------#

# base classification layer
oa = np.zeros(potentialMap.shape, dtype='uint8')

# freshwater lakes
suitable_region = region_id[np.where(swl > 0)]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 1
larr[px] = 0

# salt lakes
suitable_region = region_id[np.where(swl > 0)]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 2
larr[px] = 0

# estuaries
suitable_region = region_id[np.where(est > 0)]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 3
larr[px] = 0

# coastal freshwater lakes (including deltas)
suitable_region = region_id[np.where(dlt > 0)]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 4
larr[px] = 0

fwl = None
swl = None
est = None
dlt = None

#----------------------------------------------------------------------------------------------------------------------------#
# 6. classify water bodies based on the distinction of pixel regions (1)
#----------------------------------------------------------------------------------------------------------------------------#

# estuaries (2)
suitable_region = region_id[np.where((min_od <= 1) & (min_rd <= 2))]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 3
larr[px] = 0

# coral lagoons
suitable_region = region_id[np.where((min_cd == 0) & (min_rd > 2))]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 5
larr[px] = 0

# saline lakes (2)
suitable_region = region_id[np.where((region_area >= 8) & (min_izd >= 4) & ((bare_inside >= 0.5) & (desertic >= 0.5)))]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 2
larr[px] = 0

# frehwater lakes (2)
suitable_region = region_id[np.where((region_area >= 8) & (min_izd >= 4) & (desertic < 0.5))]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 1
larr[px] = 0

# coastal freshwater lakes
suitable_region = region_id[np.where((region_area >= 8) & (min_izd > 0) & (min_izd < 4))]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 4
larr[px] = 0

# lagoon
suitable_region = region_id[np.where(((region_area >= 8) & (min_izd == 0)) | ((min_od == 0) & (min_rd > 1)))]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 6
larr[px] = 0

# find freshwater pools
suitable_region = region_id[np.where((region_area < 8) & (min_izd >= 4) & (desertic < 0.5))]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 7
larr[px] = 0

# find saline pools
suitable_region = region_id[np.where((region_area < 8) & (min_izd >= 4) & (desertic >= 0.5))]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 8
larr[px] = 0

# find tidepools
suitable_region = region_id[np.where((region_area < 8) & (min_izd == 0))]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 9
larr[px] = 0

# find floodplains
suitable_region = region_id[np.where((region_area < 8) & (min_izd == 0))]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 10
larr[px] = 0

# isolate floodplain pixels with less than 50% permanent water
# NOTE: floodplain pixels with a high proportion of permanent water are likely lakes
larr[(oa == 10) & (oa < 0.5)] = 0
oa[(oa == 10) & (pWater < 0.5)] = 0
larr[px] = 0

min_od = None
max_od = None
min_izd = None
min_rd = None
min_cd = None

#----------------------------------------------------------------------------------------------------------------------------#
# 7. classify water bodies based on the distinction of pixel regions (2)
# NOTE: focus on potential saltine wetlands/pools
#----------------------------------------------------------------------------------------------------------------------------#

if np.max(larr) > 0:
    
    region_id = np.unique(larr[larr > 0]) # isolate missed region ID's
    region_area = nd.sum(pixelArea_hr * ((sWater+pWater)*0.01), labels=larr, index=region_id) # region area
    bare_inside = nd.sum(bareLand, labels=larr, index=region_id) / pixel_count # percent of pixels covered by bare land
    desertic = nd.sum(climate_desert, labels=larr, index=region_id) / pixel_count # percent of pixels in desertic climates
    bare_outside = np.zeros(len(region_id), dtype='float32')
    
    # iterate through each region
    for l in range(0, len(region_id)):
        m = nd.morphology.binary_dilation(larr == region_id[l], s).astype('uint8') # dilate region
        px = np.where((m == 1) & (larr != region_id[l])) # find pixels circumventing region
        bare_outside[l] = np.sum(bareLand[px]) / np.sum(m[px]) # percentage of desertic pixels circumventing the region
    
    # freshwater lakes
    suitable_region = region_id[np.where((region_area >= 8) & (bare_outside >= 0.5))]
    oa[np.isin(larr, list(suitable_region))] = 1
    
    # saline lakes
    suitable_region = region_id[np.where((region_area >= 8) & (bare_outside >= 0.5))]
    oa[np.isin(larr, list(suitable_region))] = 2
    
    # freshwater lakes
    suitable_region = region_id[np.where((region_area < 8) & (bare_outside < 0.5))]
    oa[np.isin(larr, list(suitable_region))] = 7
    
    # freshwater lakes
    suitable_region = region_id[np.where((region_area < 8) & (bare_outside >= 0.5))]
    oa[np.isin(larr, list(suitable_region))] = 8


pixelArea_hr = None
bareLand = None
larr = None

#----------------------------------------------------------------------------------------------------------------------------#
# 8. classify freshwater, coastal and wetland ecosystems
#----------------------------------------------------------------------------------------------------------------------------#

# permanent freshwater lakes
oa = (oa == 1).astype('uint8') * pWater # find ecosystem pixels
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3)) # aggregate selected pixels
oa = oa * pixelArea_agg # extract area prortion (note: proportion is based on "selected pixels / "non-nodata potential pixels")
ods = rt.open(odir + 'GlobES-0505_1km.tif', 'w', **meta) # initiate area map file
ods.write(oa, indexes=1) # write ecosystem array
ods.close() # close conenction

# seasonal freshwater lakes
oa = (oa == 1).astype('uint8') * sWater
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea_agg
ods = rt.open(odir + 'GlobES-0506_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# permanent saline lakes
oa = (oa == 2).astype('uint8') * pWater
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea_agg
ods = rt.open(odir + 'GlobES-0514_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# seasonal saline lakes
oa = (oa == 2).astype('uint8') * sWater
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea_agg
ods = rt.open(odir + 'GlobES-0515_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# estuaries
oa = (oa == 3).astype('uint8') * pWater
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea_agg
ods = rt.open(odir + 'GlobES-0910_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# coastal freshwater lakes
oa = (oa == 4).astype('uint8') * pWater
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea_agg
ods = rt.open(odir + 'GlobES-1304_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# coral lagoons
oa = (oa == 5).astype('uint8') * pWater
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea_agg
ods = rt.open(odir + 'GlobES-1305_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# lagoons (not related to corals)
oa = (oa == 6).astype('uint8') * pWater
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea_agg
ods = rt.open(odir + 'GlobES-0910_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# tidepools
oa = (oa == 9).astype('uint8') * pWater
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea_agg
ods = rt.open(odir + 'GlobES-1206_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

#----------------------------------------------------------------------------------------------------------------------------#
# 9. write layers needed to classify wetland ecosystems
#----------------------------------------------------------------------------------------------------------------------------#

# update metadata, since this output is in byte format
meta.update(dtype='uint8')

# export saline pools (used in the classification of wetlands)
oa = (oa == 2).astype('uint8')
ods = rt.open(odir + 'saline_lakes.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# export freshwater pools
oa = (oa == 7).astype('uint8')
ods = rt.open(odir + 'freshwater_pools.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# export mask of saline pools
oa = (oa == 8).astype('uint8')
ods = rt.open(odir + 'saline_pools.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()
