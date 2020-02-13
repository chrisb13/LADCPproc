#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   Date created: Thu, 30 Jan 2020 18:07:52
#   Machine created on: SB2Vbox
#


"""
Quick script to plot LADCP u-v 

To do: make it call a function to reduce repeated code.
"""
from cb2logger import *
import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np
import shareme as sm
import scipy.io
from matplotlib import gridspec
import matplotlib.pyplot as plt

import collections

def mkaxins(nplots,rows,cols):
    cnt=0
    axins=[]
    for row in range(rows):
        for col in range(cols):
            axins.append(gs[row,col])
            cnt+=1
            if cnt==nplots:
                return axins
    return axins
                


def getsites(sites):
    mats=collections.OrderedDict()
    for castnum in sites:
        try:
            mats[castnum]=scipy.io.loadmat('/home/chris/VBoxSHARED/Chris_LADCP/processing/processed/JR19002_'+castnum+'/'+castnum+'_pyexport.mat')
        except:
            lg.warning('could not find a file for: '+castnum)
            pass
    # print 'number of plotted sites',pcnt
    return mats

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

    #12 borgen, 7 marian, 11 sheldon

    ############
    #  borgen  #
    ############

    mats=getsites(borgen)
    plt.close('all')
    fig=plt.figure(figsize=(9.0,15.0))
    
    gs = gridspec.GridSpec(4, 3,hspace=.555,wspace=0.165)
    axins=mkaxins(12,4,3)

    maguv=[]
    for idx,axee in enumerate(axins):
        # print idx
        ax = plt.subplot(axee)

        ax.plot(mats[mats.keys()[idx]]['u'],mats[mats.keys()[idx]]['press'],color='b')
        ax.tick_params('x', colors='b')
        ax.set_xlabel('u-velocity (m/s)',color='b')
        ax2 = ax.twiny()
        ax2.plot(mats[mats.keys()[idx]]['v'],mats[mats.keys()[idx]]['press'],color='g')
        ax2.tick_params('x', colors='g')
        ax2.set_xlabel('v-velocity (m/s)',color='g')

        plt.gca().invert_yaxis()

        sm.inset_title_box(ax,mats.keys()[idx],bwidth="12%",location=2)

        # maguv[mats.keys()[idx]]=np.sqrt(np.square(mats[mats.keys()[idx]]['u'])+np.square(mats[mats.keys()[idx]]['v']))[:,0]
        maguv.append(np.sqrt(np.square(mats[mats.keys()[idx]]['u'])+np.square(mats[mats.keys()[idx]]['v']))[:,0].tolist())


    maguv_borgen=[]
    for m in maguv:
        maguv_borgen=m+maguv_borgen

    efile=pfol+'borgen_uv.png'
    fig.savefig(efile,dpi=300,bbox_inches='tight')
    lg.info('plot created: '+efile)


    ############
    #  marian  #
    ############

    mats=getsites(marian)
    plt.close('all')
    fig=plt.figure(figsize=(9.0,10.0))
    
    gs = gridspec.GridSpec(3, 3,hspace=.555,wspace=0.165)
    axins=mkaxins(7,3,3)

    maguv=[]
    for idx,axee in enumerate(axins):
        ax = plt.subplot(axee)

        ax.plot(mats[mats.keys()[idx]]['u'],mats[mats.keys()[idx]]['press'],color='b')
        ax.tick_params('x', colors='b')
        ax.set_xlabel('u-velocity (m/s)',color='b')
        ax2 = ax.twiny()
        ax2.plot(mats[mats.keys()[idx]]['v'],mats[mats.keys()[idx]]['press'],color='g')
        ax2.tick_params('x', colors='g')
        ax2.set_xlabel('v-velocity (m/s)',color='g')

        plt.gca().invert_yaxis()

        sm.inset_title_box(ax,mats.keys()[idx],bwidth="12%",location=2)
        maguv.append(np.sqrt(np.square(mats[mats.keys()[idx]]['u'])+np.square(mats[mats.keys()[idx]]['v']))[:,0].tolist())

    maguv_marian=[]
    for m in maguv:
        maguv_marian=m+maguv_marian

    efile=pfol+'marian_uv.png'
    fig.savefig(efile,dpi=300,bbox_inches='tight')
    lg.info('plot created: '+efile)

    #############
    #  sheldon  #
    #############
    
    mats=getsites(sheldon)
    plt.close('all')
    fig=plt.figure(figsize=(9.0,15.0))

    gs = gridspec.GridSpec(4, 3,hspace=.555,wspace=0.165)
    axins=mkaxins(11,4,3)

    maguv=[]
    for idx,axee in enumerate(axins):
        # print idx
        ax = plt.subplot(axee)

        ax.plot(mats[mats.keys()[idx]]['u'],mats[mats.keys()[idx]]['press'],color='b')
        ax.tick_params('x', colors='b')
        ax.set_xlabel('u-velocity (m/s)',color='b')
        ax2 = ax.twiny()
        ax2.plot(mats[mats.keys()[idx]]['v'],mats[mats.keys()[idx]]['press'],color='g')
        ax2.tick_params('x', colors='g')
        ax2.set_xlabel('v-velocity (m/s)',color='g')

        plt.gca().invert_yaxis()

        sm.inset_title_box(ax,mats.keys()[idx],bwidth="12%",location=2)
        maguv.append(np.sqrt(np.square(mats[mats.keys()[idx]]['u'])+np.square(mats[mats.keys()[idx]]['v']))[:,0].tolist())

    maguv_sheldon=[]
    for m in maguv:
        maguv_sheldon=m+maguv_sheldon

    efile=pfol+'sheldon_uv.png'
    fig.savefig(efile,dpi=300,bbox_inches='tight')
    lg.info('plot created: '+efile)

    plt.close('all')
    fig=plt.figure(figsize=(10.5,4.5))
    ax=fig.add_subplot(1, 3,1)
    pd.Series(maguv_sheldon).plot(kind='hist',bins=10,ax=ax)
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Mod of velocity (m/s)')
    ax.grid(True)
    sm.inset_title_box(ax,'a) sheldon',bwidth="35%",location=1)

    ax=fig.add_subplot(1, 3,2)
    pd.Series(maguv_borgen).plot(kind='hist',bins=10,ax=ax)
    ax.set_xlabel('Mod of velocity (m/s)')
    sm.inset_title_box(ax,'b) borgen',bwidth="35%",location=1)
    ax.grid(True)

    ax=fig.add_subplot(1, 3,3)
    pd.Series(maguv_marian).plot(kind='hist',bins=10,ax=ax)
    ax.set_xlabel('Mod of velocity (m/s)')
    sm.inset_title_box(ax,'b) marian',bwidth="35%",location=1)
    ax.grid(True)

    efile=pfol+'modulus_uv.png'
    fig.savefig(efile,dpi=300,bbox_inches='tight')
    lg.info('plot created: '+efile)

    print('median speeds for sheldon, borgen and marian')
    print(np.median(maguv_sheldon),np.median(maguv_borgen),np.median(maguv_marian))

    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
