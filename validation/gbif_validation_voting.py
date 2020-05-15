# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 16:07:55 2020

@author: rr70wedu
"""

#------------------------------------------------------------------------------------------------------------------#
# 0. import modules and read arguments
#------------------------------------------------------------------------------------------------------------------#

from argparse import ArgumentParser
import rasterio as rt
import pandas as pd
import numpy as np

parser = ArgumentParser(description = 'GBIF ecosystem validation (example using one ecosystem and one species)')
parser.add_argument("samples", help = "csv file with species samples and columns 'x', 'y' and 'year'")
parser.add_argument("ecosystem", help = "index")
parser.add_argument("output", help = "output directory")

options = parser.parse_args()
samples = options.samples
ecosystem = options.ecosystem
odir = options.output

# search diameter
search_diameter = [0,1,2,3,4,5,6,7,8,9]

#------------------------------------------------------------------------------------------------------------------#
# 1. read reference data
#------------------------------------------------------------------------------------------------------------------#

eds = rt.open(ecosystem) # access ecosystem map
ids = pd.read_csv(samples, dtype={'year': str}) # read samples for target species

#------------------------------------------------------------------------------------------------------------------#
# 2. validate map for different search_diameters
#------------------------------------------------------------------------------------------------------------------#

val = np.zeros((ids.shape[0], len(search_diameter)), dtype='float32') # validation vector

pr = eds.res[0] # pixel search_diameter

for i in range(0, ids.shape[0]):
    
    #================================================================================#
    # iterate through each search_diameter
    # the "diameter"/2 defines the search buffer
    #================================================================================#
    
    for p in range(0, len(search_diameter)):
        
        r = (search_diameter[p]*pr)/2 # search buffer radius
        sxy = eds.index(ids['x'][i]-r, ids['y'][i]+r) # window start x
        exy = eds.index(ids['x'][i]+r, ids['y'][i]-r) # window starty y
        nc = (exy[0]-sxy[0])+1 # number oc columns
        nr = (exy[1]-sxy[1])+1 # number of rows
        
        # determine maximum area within searc hwindow
        a = np.max(eds.read(1, window=rt.windows.Window(sxy[1], sxy[0], nc, nr)))
        val[i,p] = val[i,p] + int(a > 0) # add area results to output vector
 
eds.close()
ids = None

#------------------------------------------------------------------------------------------------------------------#
# 3. export validation results
#------------------------------------------------------------------------------------------------------------------#

v = [np.sum(val[:,i] > 0) for i in range(0, val.shape[1])] # sum correct cases for each resolution
resolution = search_diameter+1 # convert search diameter to resolution
df = pd.DataFrame({'nr_votes':v, 'resolution':resolution, 'nr_pts':ids.shape[0]})
df.to_csv(odir + 'GBIF_validation.csv', index=False)
