# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:24:06 2020

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

parser = ArgumentParser(description = 'Classification of wetland ecosystems')
parser.add_argument("variables", help = "input csv with listing input variables")
parser.add_argument("reference", help = "potential ecosystem layer with per-pixel area ratios")
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

swLakes = rt.open(variables['variable'][variables['variable'] == 'saline lakes'].values[0]).read(1) # classified in 04_openWater.py
alpine = rt.open(variables['variable'][variables['variable'] == 'alpine mountains'].values[0]).read(1)
intertidalZone = rt.open(variables['variable'][variables['variable'] == 'intertidal zone'].values[0]).read(1)
bareLand = rt.open(variables['variable'][variables['variable'] == 'bare land and sparse vegetation cover'].values[0]).read(1)
waterCover = rt.open(variables['variable'][variables['variable'] == 'water (land) cover'].values[0]).read(1)
climate_artic = rt.open(variables['variable'][variables['variable'] == 'artic climate'].values[0]).read(1)
climate_tundra = rt.open(variables['variable'][variables['variable'] == 'tundra climate'].values[0]).read(1)
climate_hdesert = rt.open(variables['variable'][variables['variable'] == 'hot desert climate'].values[0]).read(1)
climate_cdesert = rt.open(variables['variable'][variables['variable'] == 'cold desert climate'].values[0]).read(1)
climate_hsteppe = rt.open(variables['variable'][variables['variable'] == 'hot steppe climate'].values[0]).read(1)
climate_csteppe = rt.open(variables['variable'][variables['variable'] == 'cold steppe climate'].values[0]).read(1)
sWater = rt.open(variables['variable'][variables['variable'] == 'seasonal water'].values[0]).read(1)
pWater = rt.open(variables['variable'][variables['variable'] == 'permanent water'].values[0]).read(1)

# make mask of desertic climates
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

pixelArea = rt.open(variables['variable'][variables['variable'] == 'pixel area'].values[0]).read(1)

# aggregate pixel area to 1km
pixelArea = pixelArea.reshape(onr, factor, onc, factor).sum(axis=(1, 3))

#----------------------------------------------------------------------------------------------------------------------------#
# 4. classify ecosystem types at 300m, aggregate to 1km and write output (focus on rivers)
#----------------------------------------------------------------------------------------------------------------------------#

# mud flats
oa = waterCover * intertidalZone * pWater * pixelArea
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
ods = rt.open(odir + 'IUCN_Habitat-1204_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()
potentialMap[oa > 0] = 0

# salt marshes
oa = potentialMap * intertidalZone
ods = rt.open(odir + 'GlobES-1205_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()
potentialMap[oa > 0] = 0

intertidalZone = None

# alpine wetlands
oa = potentialMap * alpine
ods = rt.open(odir + 'GlobES-0511_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()
potentialMap[oa > 0] = 0

alpine = None

# tundra
oa = potentialMap * ((climate_artic > 0) | (climate_tundra > 0)).astype('uint8')
ods = rt.open(odir + 'GlobES-0510_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()
potentialMap[oa > 0] = 0

climate_artic = None
climate_tundra = None

#----------------------------------------------------------------------------------------------------------------------------#
# 5. find individual pixel regions to classify and derive region statistics
#----------------------------------------------------------------------------------------------------------------------------#

# label connected pixel regions
larr, nf = nd.measurements.label((potentialMap > 0).astype('uint8'), s)
region_id = np.array(range(1, nf+1)) # unique region ID's

# count of pixels per region
pixel_count = np.unique(larr[larr > 0], return_counts=True)[1]

# check for basic statistic applicable to all regions
region_area = nd.sum(pixelArea * potentialMap, labels=larr, index=region_id) # region area
bare_inside = nd.sum(bareLand, labels=larr, index=region_id) / pixel_count # percent of pixels covered by bare land
desertic = nd.sum(climate_desert, labels=larr, index=region_id) / pixel_count # percent of pixels in desertic climates
saltLake = nd.sum(swLakes, labells=larr, index=region_id) / pixel_count # occurrence of salt lakes

swLakes = None

#----------------------------------------------------------------------------------------------------------------------------#
# 6. classify known water bodies
#----------------------------------------------------------------------------------------------------------------------------#

# base classification layer
oa = np.zeros(potentialMap.shape, dtype='uint8')

# find saline wetlands over known saline lakes
suitable_region = region_id[np.where(saltLake > 0)]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 1
larr[px] = 0

# find saline wetlands were the area is dominated by bare land and sparse vegetation (likely related to pools)
suitable_region = region_id[np.where(bare_inside > 0.5)]
px = np.where(np.isin(larr, list(suitable_region)))
oa[px] = 1
larr[px] = 0

#----------------------------------------------------------------------------------------------------------------------------#
# 7. classify water bodies based on the distinction of pixel regions (2)
# NOTE: focus on potential saltine wetlands/pools
#----------------------------------------------------------------------------------------------------------------------------#

if np.max(larr) > 0:
    
    region_id = np.unique(larr[larr > 0]) # isolate missed region ID's
    region_area = nd.sum(pixelArea * ((sWater+pWater)*0.01), labels=larr, index=region_id) # region area
    desertic = nd.sum(climate_desert, labels=larr, index=region_id) / pixel_count # percent of pixels in desertic climates
    bare_outside = np.zeros(len(region_id), dtype='float32')
    
    # iterate through each region
    for l in range(0, len(region_id)):
        m = nd.morphology.binary_dilation(larr == region_id[l], s).astype('uint8') # dilate region
        px = np.where((m == 1) & (larr != region_id[l])) # find pixels circumventing region
        bare_outside[l] = np.sum(bareLand[px]) / np.sum(m[px]) # percentage of desertic pixels circumventing the region
    
    # saline wetlands
    suitable_region = region_id[np.where((bare_outside >= 0.5) & (desertic >= 0.5))]
    oa[np.isin(larr, list(suitable_region))] = 1
    
    # freshwater wetlands (>= 8ha)
    suitable_region = region_id[np.where(region_area >= 8)]
    oa[(oa == 0) & np.isin(larr, list(suitable_region))] = 2
    
    # freshwater wetlands (< 8ha)
    suitable_region = region_id[np.where(region_area < 8)]
    oa[(oa == 0) & np.isin(larr, list(suitable_region))] = 3

larr = None    
bareLand = None
region_area = None
desertic = None

#----------------------------------------------------------------------------------------------------------------------------#
# 8. classify freshwater, coastal and wetland ecosystems
#----------------------------------------------------------------------------------------------------------------------------#

# freshwater wetlands (>= 8 ha)
oa = (oa == 2).astype('uint8') # find ecosystem pixels
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3)) # aggregate selected pixels
oa = oa * pixelArea # extract area prortion (note: proportion is based on "selected pixels / "non-nodata potential pixels")
ods = rt.open(odir + 'GlobES-0504_1km.tif', 'w', **meta) # initiate area map file
ods.write(oa, indexes=1) # write ecosystem array
ods.close() # close conenction

# permanent freshwater pools
oa = (oa == 3).astype('uint8') * pWater
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea
ods = rt.open(odir + 'GlobES-0507_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# seasonal freshwater pools
oa = (oa == 3).astype('uint8') * sWater
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea
ods = rt.open(odir + 'GlobES-0508_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# permanent saline wetlands
oa = (oa == 1).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * (pixelArea*sWater)
ods = rt.open(odir + 'GlobES-0516_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# seasonal freshwater lakes
oa = (oa == 1).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).mean(axis=(1, 3))
oa = oa * pixelArea
ods = rt.open(odir + 'GlobES-0517_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()
