# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 21:22:56 2019

@author: rr70wedu
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:22:02 2019

@author: rr70wedu
"""

#----------------------------------------------------------------------------------------------------------------------------#
# 0 load modules and arguments
#----------------------------------------------------------------------------------------------------------------------------#

from argparse import ArgumentParser
import rasterio as rt
import pandas as pd

parser = ArgumentParser(description = 'Classification of oceanic ecosystem types')
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

oceanDistance = rt.open(variables['variable'][variables['variable'] == 'ocean distance'].values[0]).read(1)
elevation = rt.open(variables['variable'][variables['variable'] == 'elevation'].values[0]).read(1)
slope = rt.open(variables['variable'][variables['variable'] == 'slope'].values[0]).read(1)
abyssalPlains = rt.open(variables['variable'][variables['variable'] == 'abyssal plains'].values[0]).read(1)
abyssalMountains = rt.open(variables['variable'][variables['variable'] == 'abyssal mountains'].values[0]).read(1)
seamounts = rt.open(variables['variable'][variables['variable'] == 'Seamounts'].values[0]).read(1)
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

# epipelagic
oa = potentialMap * (elevation < 0).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3)) # aggregate selected pixels
oa = (oa/maxProp) * pixelArea # extract area prortion (note: proportion is based on "selected pixels / "non-nodata potential pixels")
ods = rt.open(odir + 'GlobES-1001_1km.tif', 'w', **meta) # initiate area map file
ods.write(oa, indexes=1) # write ecosystem array
ods.close() # close conenction

# mesopelagic
oa = potentialMap * (elevation < -200).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-1002_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# bathypelagic
oa = potentialMap * (elevation < -1000).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-1003_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# abyssopelagic
oa = potentialMap * (elevation < -4000).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-1004_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# bathyl zone
oa = potentialMap * ((elevation <= -200) & (elevation >= -2000) & (slope > 0)).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-1101_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# abyssal plains
oa = potentialMap * ((abyssalPlains > 0) & (slope == 0) & (elevation <= -4000)).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-1102_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# abyssal mountains
oa = potentialMap * (abyssalMountains > 0).astype('uint8') * (slope > 0) * (elevation <= -4000)
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-1103_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# hadal
oa = potentialMap * (elevation < -6000).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-1104_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()

# seamounts
oa = potentialMap * (seamounts  > 0).astype('uint8') * (slope > 0)
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-1105_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()
