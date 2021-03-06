# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 12:57:17 2019

@author: rr70wedu
"""

#----------------------------------------------------------------------------------------------------------------------------#
# 0 load modules and arguments
#----------------------------------------------------------------------------------------------------------------------------#

from argparse import ArgumentParser
import rasterio as rt
import pandas as pd

parser = ArgumentParser(description = 'Classification of (dry) forest ecosystem types')
parser.add_argument("variables", help = "input csv with listing input variables")
parser.add_argument("reference", help = "potential ecosystem layer with per-pixel area ratios")
parser.add_argument("output", help = "path to output directory")

options = parser.parse_args()
reference = options.reference
variables = pd.read_csv(options.variables)
odir = options.output

factor = 3 # aggregation factor for total area estimation

#----------------------------------------------------------------------------------------------------------------------------#
# 1. load data for input variables
#----------------------------------------------------------------------------------------------------------------------------#

climate_boreal = rt.open(variables['variable'][variables['variable'] == 'boreal climate'].values[0]).read(1)
climate_artic = rt.open(variables['variable'][variables['variable'] == 'artic climate'].values[0]).read(1)
climate_tundra = rt.open(variables['variable'][variables['variable'] == 'tundra climate'].values[0]).read(1)
climate_temperate = rt.open(variables['variable'][variables['variable'] == 'temperate climate'].values[0]).read(1)
climate_mediterranean = rt.open(variables['variable'][variables['variable'] == 'mediterranean climate'].values[0]).read(1)
climate_sdry = rt.open(variables['variable'][variables['variable'] == 'sub-/tropical dry climate'].values[0]).read(1)
climate_smoist = rt.open(variables['variable'][variables['variable'] == 'sub-/tropical moist climate'].values[0]).read(1)
climate_savanna = rt.open(variables['variable'][variables['variable'] == 'savanna climate'].values[0]).read(1)
climate_hdesert = rt.open(variables['variable'][variables['variable'] == 'hot desert'].values[0]).read(1)
climate_cdesert = rt.open(variables['variable'][variables['variable'] == 'cold desert'].values[0]).read(1)
climate_hsteppe = rt.open(variables['variable'][variables['variable'] == 'hot steppe climate'].values[0]).read(1)
climate_csteppe = rt.open(variables['variable'][variables['variable'] == 'cold steppe climate'].values[0]).read(1)
mountains = rt.open(variables['variable'][variables['variable'] == 'mountains'].values[0]).read(1)
islands = rt.open(variables['variable'][variables['variable'] == 'artic islands'].values[0]).read(1)
latitude  = rt.open(variables['variable'][variables['variable'] == 'latitude'].values[0]).read(1)
sWater = rt.open(variables['variable'][variables['variable'] == 'seasonal water'].values[0]).read(1)
pWater = rt.open(variables['variable'][variables['variable'] == 'permanent water'].values[0]).read(1)

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
# 4. classify ecosystem types at 300m, aggregate to 1km and write output
#----------------------------------------------------------------------------------------------------------------------------#

# maximum (possible) proportion of available area, discounting missing data pixels (used to adjust pixel areas at 1km)
maxProp = ((potentialMap > 0) & (potentialMap != nodata)).astype('float32').reshape(onr, factor, onc, factor).sum(axis=(1, 3))
potentialMap[potentialMap == nodata] = 0

# boreal
oa = potentialMap * climate_boreal * (latitude > 0).astype('uint8') # find ecosystem pixels
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3)) # aggregate selected pixels
oa = (oa/maxProp) * pixelArea # extract area prortion (note: proportion is based on "selected pixels / "non-nodata potential pixels")
ods = rt.open(odir + 'GlobES-0102_1km.tif', 'w', **meta) # initiate area map file
ods.write(oa, indexes=1) # write ecosystem array
ods.close() # close conenction

oa = None
climate_boreal = None

# subarctic
oa = potentialMap * ((climate_artic > 0) | (climate_tundra > 0)).astype('uint8') * (latitude > 0).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-0102_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

oa = None
latitude = None

# subantarctic
oa = potentialMap * islands * ((climate_artic > 0) | (climate_tundra > 0)).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3)) 
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-0103_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

oa = None
islands = None
climate_artic = None
climate_tundra = None

# temperate
oa = potentialMap * ((climate_temperate > 0) | (climate_mediterranean > 0) | (climate_csteppe > 0) | ((climate_cdesert > 0) & (mountains > 0))).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-0104_1km.tif', 'w', **meta) 
ods.write(oa, indexes=1)
ods.close()

oa = None
climate_temperate = None
climate_mediterranean = None

# subtropical/tropical dry
oa = potentialMap * ((climate_savanna > 0) | (climate_cdesert > 0) | (climate_hsteppe > 0) | (climate_sdry > 0)).astype('uint8') * (pWater == 0).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-0105_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

oa = None

# subtropical/tropical moist (lowland)
oa = potentialMap * ((climate_smoist > 0) | (((climate_sdry > 0) | (climate_cdesert > 0) | (climate_hsteppe > 0) | (climate_savanna > 0)) & (pWater > 0))).astype('uint8') * (mountains == 0).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-0106_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

oa = None

# subtropical/tropocal moist (mountain)
oa = potentialMap * ((climate_smoist > 0) | (((climate_sdry > 0) | (climate_cdesert > 0) | (climate_hsteppe > 0) | (climate_savanna > 0)) & (pWater > 0)) | (climate_cdesert > 0)).astype('uint8') * (mountains > 0).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-0109_1km.tif', 'w', **meta)
ods.write(oa.astype('uint8'), indexes=1)
ods.close()

oa = None
pw = None
