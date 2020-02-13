#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   www:     christopherbull.com.au
#   Date created: Wed, 22 Jan 2020 17:28:27
#   Machine created on: SB2Vbox
#

"""
Quick script to plot locations of casts
"""
from cb2logger import *
import shareme as sm

import os
import xarray as xr
import collections
import matplotlib.pyplot as plt


if __name__ == "__main__": 
    LogStart('',fout=False)

    # infile='/home/chris/VBoxSHARED/Chris_LADCP/plot_casts/gebco_1min.grd'
    # assert(os.path.exists(infile)),"netCDF file does not exist!"
    # ifile=xr.open_dataset(infile)
    # __import__('pdb').set_trace()


    import shapefile
    sf = shapefile.Reader("/home/chris/VBoxSHARED/Chris_LADCP/trackPoint_JR19002.shp")
    __import__('pdb').set_trace()

    plt.close('all')
    fig=plt.figure()
    ax=fig.add_subplot(1, 1,1)
    ax.contourf(ifile['z'])
    plt.show()
    __import__('pdb').set_trace()


    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
