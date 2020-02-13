#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   Date created: Mon, 27 Jan 2020 15:20:43
#   Machine created on: SB2Vbox
#
"""
Quick script to try and unravel what's going wrong with CTDs:
        bads=[13, 15, 16, 23, 41, 42, 50, 60, 61, 63]

Solution: needed to manually crop with WinADCP, b/c i_ladcp_near_bottom was not working (see plots and script output).
"""
from cb2logger import *
import glob
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
import datetime
import numpy as np
import shareme as sm
import glob

if __name__ == "__main__": 
    LogStart('',fout=False)

    pfol='./plots/bads/'
    sm.mkdir(pfol)

    f='/home/chris/VBoxSHARED/Chris_LADCP/plot_casts/013_pyexport.mat'
    f='/home/chris/VBoxSHARED/Chris_LADCP/plot_casts/012_pyexport.mat'

    ifiles=sorted(glob.glob('/home/chris/VBoxSHARED/Chris_LADCP/plot_casts/inspectbads/0*_pyexport*'))
    assert(ifiles!=[]),"glob didn't find anything!"
    for f in ifiles:
        mat = scipy.io.loadmat(f)

        bads=[13, 15, 16, 23, 41, 42, 50, 60, 61, 63]
        for bad in bads:
            if str(bad) in os.path.basename(f):
                co='red'
                break
            else:
                co=None

        plt.close('all')
        fig=plt.figure()
        fig=plt.figure(figsize=(14.0,4.0))
        ax=fig.add_subplot(1, 4,1)
        ax.plot(mat['z'].flatten(),color=co)

        try:
            idx=mat['i_ladcp_near_bottom'][0][0]
            ax.scatter(idx,mat['z'][0][idx])
        except Exception as e:
            # raise e
            pass

        ax.set_title('z.  '+ '0.9*zmax = '+str(np.round(mat['zmax'][0][0]*0.9,1)))
        
        ax=fig.add_subplot(1, 4,2)
        ax.plot(mat['timctd'].flatten(),color=co)
        ax.set_title('timctd')
        ax=fig.add_subplot(1, 4,3)
        ax.plot(mat['dt'].flatten(),color=co)
        ax.set_title('dt')
        ax=fig.add_subplot(1, 4,4)
        ax.plot(mat['w'].flatten(),color=co)
        ax.set_title('w')
        # plt.show()

        efile=pfol+os.path.basename(f)+'.png'
        fig.savefig(efile,dpi=300,bbox_inches='tight')
        # lg.info('plot created: '+efile)


        try:
            lg.info(os.path.basename(f)+': '+'Delta_t: ' + str(np.round(mat['delta_t'][0][0],1))+ '. ctd/ladcp (resp) near bottom: '+ str(mat['i_ctd_near_bottom'][0][0]) +'/' + str(mat['i_ladcp_near_bottom'][0][0])+'. zmax and pmax:'+str(np.round(mat['zmax'][0][0],1))+'/'+str(np.round(mat['pmax'][0][0],1)))

        except IndexError as e:
            # raise e
            lg.error(os.path.basename(f)+': i_ladcp_near_bottom: '+str(mat['i_ladcp_near_bottom'][0])+'. zmax and pmax:'+str(np.round(mat['zmax'][0][0],1))+'/'+str(np.round(mat['pmax'][0][0],1)))

    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
