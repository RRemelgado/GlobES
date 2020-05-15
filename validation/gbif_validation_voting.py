# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 16:07:55 2020

@author: rr70wedu
"""

#-----------------------------------------------------------------------------#
#
#-----------------------------------------------------------------------------#

from argparse import ArgumentParser
from os.path import basename as b
from os.path import isfile as fi
import rasterio as rt
import pandas as pd
import numpy as np
import glob2 as g

# read path to output directory
parser = ArgumentParser(description = 'quantify ecossystem area')
parser.add_argument("index", help = "file index")
options = parser.parse_args()
index = int(options.index)

bsize = [0,1,2,3,4,5,6,7,8,9]

idir = '/work/remelgad/IUCN_ecoMap_30as/'
odir1 = '/work/remelgad/IUCN_ecoMap_val/gbif/all/'
odir2 = '/work/remelgad/IUCN_ecoMap_val/gbif/summary/'

rds = rt.open(idir + 'IUCN_ecoMap_30as-0101_19920101_30arcSec.vrt')
pr = rds.res[0]

#-----------------------------------------------------------------------------#
#
#-----------------------------------------------------------------------------#

files = g.glob('/work/remelgad/IUCN_ecoMap_val/gbif/tmp/*.csv')
f = files[index]

fbn = b(f).split('.')[0]
oname1 = odir1 + fbn + '_all.csv'
oname2 = odir2 + fbn + '_summary.csv'

ifile = '/work/remelgad/IUCN_ecoMap_val/gbif/IUCN_ecoMap_val-gbif_selectedSpecies.csv'
ref = pd.read_csv(ifile, dtype={'cid': str}, low_memory=False)

#-----------------------------------------------------------------------------#
#
#-----------------------------------------------------------------------------#

if not fi(oname1) and not fi(oname2):
    
    #=========================================================================#
    # iterate through each range
    #=========================================================================#
    
    ids = pd.read_csv(f, dtype={'year': str})
    n = np.zeros(len(bsize), dtype='float32')
    n[:] = ids.shape[0] # number of observations
    
    eco = list(np.unique(ref[ref['binomial'] == ids['species'][0]]['cid']))
    
    val = np.zeros((ids.shape[0], len(bsize)), dtype='float32')
    
    if len(eco) > 0:
        
        for i in range(0, ids.shape[0]):
            
            #=====================================================================#
            # iterate through each range
            #=====================================================================#
            
            for p in range(0, len(bsize)):
                
                r = (bsize[p]*pr)/2
                sxy = rds.index(ids['longitude'][i]-r, ids['latitude'][i]+r)
                exy = rds.index(ids['longitude'][i]+r, ids['latitude'][i]-r)
                
                nc = (exy[0]-sxy[0])+1
                nr = (exy[1]-sxy[1])+1
                
                for c in eco:
                    
                    eds = rt.open(idir + 'IUCN_ecoMap_30as-' + c + '_' + ids['year'][i] + '0101_30arcSec.vrt')
                    a = np.max(eds.read(1, window=rt.windows.Window(sxy[1], sxy[0], nc, nr)))
                    val[i,p] = val[i,p] + int(a > 0)
                    eds.close()
                
                val[i,p] = (val[i,p] > 0).astype('uint8')
                
                print(p)
        
        #-----------------------------------------------------------------------------#
        # write valdation per observation
        #-----------------------------------------------------------------------------#
        
        fbn = b(f).split('.')[0]
        df = pd.DataFrame(val)
        df.to_csv(oname1, index=False)
        
        #-----------------------------------------------------------------------------#
        # write overall validation
        #-----------------------------------------------------------------------------#
        
        s = ids['species'][0] # species name
        v = [np.sum(val[:,i]) for i in range(0, val.shape[1])] # sum of correct cases
        e = len(eco) # number of ecosystems
        df = pd.DataFrame({'species':s, 'val':v, 'buffer':bsize, 'nr_eco':e, 'nr_pts':n})
        df.to_csv(oname2, index=False)
        
        rds.close()

print(b(f) + ' done!')
