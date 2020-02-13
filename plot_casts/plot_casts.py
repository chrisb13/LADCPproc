#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   www:     christopherbull.com.au
#   Date created: Tue, 14 Jan 2020 20:58:32
#   Machine created on: SB2Vbox
#
"""
Quick script to plot the processed LADCP casts for JR19002
"""
from cb2logger import *
import scipy.io
import os
import matplotlib.pyplot as plt
import numpy as np

import glob
import collections
import shareme as sm

if __name__ == "__main__": 
    LogStart('',fout=False)

    pfol='./plots/'
    sm.mkdir(pfol)

    ifiles=sorted(glob.glob('/home/chris/VBoxSHARED/Chris_LADCP/processing/processed/*/*pyexport.mat' ))

    lg.info('Found ' + str(len(ifiles)) + ' files to plot.')
    assert(ifiles!=[]),"glob didn't find anything!"

    for infile in ifiles:
        mat = scipy.io.loadmat(infile)

        if 'ctdt' in mat.keys():
            ctd=True
        else:
            ctd=False

        plt.close('all')
        fig=plt.figure()
        if not ctd:
            ax=fig.add_subplot(1, 2,1)
        else:
            ax=fig.add_subplot(2, 2,1)
        ax.plot(mat['u'],mat['press'])
        plt.gca().invert_yaxis()
        ax.grid(True)
        ax.set_ylabel('Pressure')
        ax.set_xlabel('u-velocity')
        ax.set_title(os.path.basename(infile)[0:3])

        if not ctd:
            ax=fig.add_subplot(1, 2,2)
        else:
            ax=fig.add_subplot(2, 2,2)
        ax.plot(mat['v'],mat['press'])
        plt.gca().invert_yaxis()
        ax.grid(True)
        ax.set_xlabel('v-velocity')
        ax.set_title(str(np.round(mat['lat'][0].tolist()[0],1))+' ,'+str(np.round(mat['lon'][0].tolist()[0],1)))
        # plt.show()

        if ctd:
            ax=fig.add_subplot(2, 2,3)
            ax.plot(mat['ctdt'],mat['press'])
            plt.gca().invert_yaxis()
            ax.grid(True)
            ax.set_xlabel('temp')

            ax=fig.add_subplot(2, 2,4)
            ax.plot(mat['ctds'],mat['press'])
            plt.gca().invert_yaxis()
            ax.grid(True)
            ax.set_xlabel('sal')

        efile=pfol+os.path.basename(infile)[0:3]+'.png'

        fig.savefig(efile,dpi=300,bbox_inches='tight')
        lg.info('plot created: '+efile)
        
    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
