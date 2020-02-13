#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   Date created: Tue, 28 Jan 2020 23:20:25
#   Machine created on: SB2Vbox
#
"""
Quick script to create the summary table for ICEBERGS (JR19002) cruise report
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

    report_tab=collections.OrderedDict()
    for key in sm.cast_type.keys():
        # print key
        castnum=key

        #skip racetrax
        if sm.cast_type[key][0]=='racetrax':
            continue

        try:
            mat = scipy.io.loadmat('/home/chris/VBoxSHARED/Chris_LADCP/processing/processed/JR19002_'+castnum+'/'+castnum+'_pyexport.mat')
            # castnum=os.path.basename(infile)[0:3]
            # print(os.path.basename(infile)[0:3]),sm.cast_type[os.path.basename(infile)[0:3]]
            report_tab[castnum]=sm.cast_type[castnum]
            report_tab[castnum].append(np.round(mat['lat'][0][0],3))
            report_tab[castnum].append(np.round(mat['lon'][0][0],3))
            try:
                report_tab[castnum].append('-'.join([str(m) for m in mat['date'][0][0:3]]) + ' ~ '+':'.join([str(m) for m in mat['date'][0][3:]]))
            except:
                report_tab[castnum].append('')
            # print(report_tab[castnum])
            # print(report_tab[castnum][0],report_tab[castnum][1],report_tab[castnum][2],report_tab[castnum][3])
            # print report_tab[castnum][0]
        except IOError:
            report_tab[castnum]=sm.cast_type[castnum]

        # print castnum,report_tab[castnum][2],report_tab[castnum][3],report_tab[castnum][1],report_tab[castnum][0],report_tab[castnum][4]
        try:
            print castnum,report_tab[castnum][2],report_tab[castnum][3],report_tab[castnum][1],report_tab[castnum][0],report_tab[castnum][4]
        except:
            print report_tab.values()[-1][0],report_tab.values()[-1][1]
        # __import__('pdb').set_trace()
        
        
    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
