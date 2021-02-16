# -*- coding: utf-8 -*-
"""
Created on Thu May 14 17:40:29 2020

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

parser = ArgumentParser(description = 'Classification of artificial ecosystem types related to agriculture')
parser.add_argument("variables", help = "input csv with listing input variables")
parser.add_argument("output", help = "path to output directory")

options = parser.parse_args()
reference = options.reference
variables = pd.read_csv(options.variables)
odir = options.output

factor = 3 # aggregation factor for total area estimation

#----------------------------------------------------------------------------------------------------------------------------#
# 1. load data for input variables and metadata
#----------------------------------------------------------------------------------------------------------------------------#

arableLand = rt.open(variables['variable'][variables['variable'] == 'arable land'].values[0])
pastureLand = rt.open(variables['variable'][variables['variable'] == 'pastureland'].values[0])
sWater = rt.open(variables['variable'][variables['variable'] == 'seasonal water'].values[0]).read(1)
pWater = rt.open(variables['variable'][variables['variable'] == 'permanent water'].values[0]).read(1)

# extract/update metadata (needed to write output)
meta = arableLand.meta.copy()
nodata1 = arableLand.nodata # no data value for arable land
nodata2 = pastureLand.nodata # no data value for arable land
nc = arableLand.height # number of rows (original)
onc = nc/factor # output number of rows (output)
nr = arableLand.width # number of columns
onr = nr/factor # output number of columns (output)
meta.update(driver='GTiff', dtype='float32', compress='deflate', predict=2, zlevel=9, nc=onc, nr=onr, nodata=255)

# read arable/pasture land data
arableLand = arableLand.read(1)
pastureLand = pastureLand.read(1)

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
maxProp = (((arableLand > 0) | (pastureLand > 0)) & (arableLand != nodata1) & (pastureLand != nodata2)).astype('float32').reshape(onr, factor, onc, factor).sum(axis=(1, 3))
arableLand[arableLand == nodata1] = 0
pastureLand[pastureLand == nodata2] = 0

# Arable Land
oa = arableLand * (sWater == 0).astype('uint8') # find ecosystem pixels
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3)) # aggregate selected pixels
oa = (oa/maxProp) * pixelArea # extract area prortion (note: proportion is based on "selected pixels / "non-nodata potential pixels")
ods = rt.open(odir + 'GlobES-1401_1km.tif', 'w', **meta) # initiate area map file
ods.write(oa, indexes=1) # write ecosystem array
ods.close() # close conenction

a = None

# pastureland
oa = pastureLand * (sWater == 0).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-1402_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()


a = None
oa = None

# seasonally flooded agriculture
oa = (arableLand + pastureLand) * (sWater > 0).astype('uint8')
oa = oa.reshape(onr, factor, onc, factor).sum(axis=(1, 3))
oa = (oa/maxProp) * pixelArea
ods = rt.open(odir + 'GlobES-1508_1km.tif', 'w', **meta)
ods.write(oa, indexes=1)
ods.close()
