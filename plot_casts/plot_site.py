#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   www:     christopherbull.com.au
#   Date created: Thu, 30 Jan 2020 14:15:20
#   Machine created on: SB2Vbox
#

"""
Quick script to plot CTD site locations using coastline data from Alice and Kate
"""
from cb2logger import *
import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np
import shareme as sm
import scipy.io

def plotsites(sites,pltaxis):
    pcnt=0
    for castnum in sites:
        try:
            mat = scipy.io.loadmat('/home/chris/VBoxSHARED/Chris_LADCP/processing/processed/JR19002_'+castnum+'/'+castnum+'_pyexport.mat')
            pltaxis.scatter(mat['lon'][0][0],mat['lat'][0][0],label=castnum+'_'+sm.cast_type[castnum][1])
            # pltaxis.quiver(np.mean(mat['v']),np.mean(mat['v']))
            pcnt+=1
        except:
            lg.warning('could not find a file for: '+castnum)
            pass
    print 'number of plotted sites',pcnt
    return pltaxis

if __name__ == "__main__": 
    LogStart('',fout=False)

    pfol='./plots/'
    sm.mkdir(pfol)

    fol='/home/chris/VBoxSHARED/Chris_LADCP/plot_casts/Coastline_csv/'
    ifiles=sorted(glob.glob(fol+'*.csv'))
    assert(ifiles!=[]),"glob didn't find anything!"

    # # Sheldon_Cove
    sheldon=['002','005','008','009','012','013','014','015','016','017','020','023']

    # # Borgen bay
    borgen=['024','027','028','031','034','035','036','037','038','041','044','045','046']

    # # marian cove
    marian=['047','049','050','051','053','055','057','058','060','061','063']

    for f in ifiles:
        df=pd.read_csv(f)

        if os.path.basename(f)=='Borgen_wgs84.csv':
            plt.close('all')
            fig=plt.figure()
            ax=fig.add_subplot(1, 1,1)
            ax.set_title('Borgen Bay')
            ax.grid(True)
            borgen=['024','027','028','031','034','035','036','037','038','041','044','045','046']
            plotsites(borgen,ax)
            ax.legend()
            ax.scatter(df['Longitude'],df['Latitude'],color='k')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')

        elif os.path.basename(f)=='Marian_wgs84.csv':
            plt.close('all')
            fig=plt.figure()
            ax=fig.add_subplot(1, 1,1)
            ax.set_title('Marian Cove')
            ax.grid(True)
            plotsites(marian,ax)
            ax.legend()
            ax.scatter(df['Longitude'],df['Latitude'],color='k')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
        elif os.path.basename(f)=='Sheldon_wgs84.csv':
            plt.close('all')
            fig=plt.figure()
            ax=fig.add_subplot(1, 1,1)
            # ax.plot(df['Longitude'],df['Latitude'])
            ax.set_title('Sheldon Cove')
            ax.grid(True)
            plotsites(sheldon,ax)
            ax.legend(loc=2)
            ax.scatter(df['Longitude'],df['Latitude'],color='k')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')

        efile=pfol+os.path.basename(f)+'.png'
        fig.savefig(efile,dpi=300,bbox_inches='tight')
        lg.info('plot created: '+efile)

    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
