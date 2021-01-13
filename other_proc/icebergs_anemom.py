#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  Department of Geography and Environmental Sciences, 
#                 Northumbria University, Newcastle upon Tyne, UK
#   Contact: christopher.bull@northumbria.ac.uk
#   Date created: Thu, 23 Apr 2020 16:14:29
#   Machine created on: SB2Vbox
#

"""
Mike M's request to look at ICEBERGS wind data in terms of wind induced mixing

from:	Meredith, Michael P. <mmm@bas.ac.uk>
to:	"Bull, Christopher Y.S." <chbull@bas.ac.uk>
date:	Apr 20, 2020, 4:38 PM
subject:	Re: Quick request for update

Data are in columns, being:-
Year
Days elapsed since start of year
Seconds elapsed since start of year
Wind direction (relative to the ship)
Wind speed (relative to the ship, in knots)
Wind speed (relative to the ship, in m/s)
Wind velocity eastward (relative to the ship)
Wind velocity northward (relative to the ship)
True wind speed (i.e. relative to the solid earth) in m/s
True wind direction TO (direction wind blowing to)
True wind direction FROM (direction wind blowing from)
"""
import sys,os
from cb2logger import *
import glob
import scipy.io
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#conda install -c conda-forge windrose
from windrose import WindroseAxes

def mklne(pltax,tseries,drawmedian=False,drawpercent=False):
    #for rtime in rece_idx:
    pltax.vlines(x=tseries.index[350000], ymin=np.min(tseries), ymax=np.max(tseries),linewidth=6, color='green', zorder=1,alpha=0.5)

    ax.plot(tseries.rolling(window=500,center=True,win_type='hamming').mean(),label='rolling 500 hamming tsteps',color='k')

    #if drawpercent:
        #pme=tseries[tseries>np.percentile(tseries,90)]
        #pltax.plot(pme.index, pme.values,label='> 90 percentile',color='blue',lw=1.2)

    pltax.hlines(y=np.percentile(tseries,90), xmin=tseries.index[0], xmax=tseries.index[-1],linewidth=3, color='yellow', zorder=1,alpha=0.9)
    pltax.hlines(y=np.percentile(tseries,80), xmin=tseries.index[0], xmax=tseries.index[-1],linewidth=3, color='yellow', zorder=1,alpha=0.9)
    pltax.hlines(y=np.percentile(tseries,70), xmin=tseries.index[0], xmax=tseries.index[-1],linewidth=3, color='yellow', zorder=1,alpha=0.9)
    pltax.hlines(y=np.percentile(tseries,60), xmin=tseries.index[0], xmax=tseries.index[-1],linewidth=3, color='yellow', zorder=1,alpha=0.9)

    if drawmedian:
        pltax.hlines(y=np.median(tseries), xmin=tseries.index[0], xmax=tseries.index[-1],linewidth=3, color='yellow', zorder=1,alpha=0.9)
    return

def grab_data(ifol):
    """@todo: Docstring for grab_data
    
    :arg1: @todo
    :returns: @todo
    """
    ifiles=sorted(glob.glob(ifol+'anemom*.txt'))
    assert(ifiles!=[]),"glob didn't find anything!"
    #infile='./anemom/anemom017true.txt'
    for infile in ifiles:
        efile = ifol+os.path.basename(infile)[:-4] +'_table'+ '.h5'

        if os.path.exists(efile):
            lg.info(infile+' already exists, so skipping creation')
            continue

        lg.info('doing :'+infile)

        idx=['Year','Days_elapsed','Seconds_elapsed','Wind_dir_ship','Wind_sp_kn','Wind_sp_ms','Wind_vel_east','Wind_vel_north','True_wind_speed_ms','True_wind_dir_TO','True_wind_dir_FROM']

        df=pd.read_csv(infile,header=None,delim_whitespace=True,names=idx)
        #df=df.ix[0:10]
        # from Mike M:
        #I think yours is offset by 1 day... note that in my plot (attached) the big lump of winds sits between 19.5 and 20, and in yours it sits between 20.5 and 21. This will be duie to the timing convention used...
        df.index=[pd.to_datetime(str(int(yy)),format="%Y")+pd.Timedelta(ss,unit='s')-pd.Timedelta(1,unit='D') for yy,ss in df[['Year','Seconds_elapsed']].values]

        try:
            os.remove(efile)
            lg.info("HDFStore already existed, clobbering!")
        except OSError:
            pass 
        
        store = pd.HDFStore(efile,complevel=9, complib='blosc')
        store.put('df',df)
        store.close()


    ifiles=sorted(glob.glob(ifol+'anemom*.h5'))
    assert(ifiles!=[]),"glob didn't find anything!"
    df=pd.concat([pd.HDFStore(infile).select('df') for infile in ifiles])
    return df


if __name__ == "__main__": 
    LogStart('',fout=False)

    df=grab_data('/home/chris/VBoxSHARED/anemom/')

    plt.close('all')
    fig=plt.figure(figsize=(15.0,5.0))
    ax=fig.add_subplot(1, 1,1)
    #df['True_wind_dir_TO'].plot(color='b')
    #df['True_wind_dir_FROM'].plot(color='r')
    #df['True_wind_speed_ms'].plot(color='r',ax=ax)
    ax.plot(df['True_wind_speed_ms'],label='True_wind_speed_ms',color='r')
    mklne(ax,df['True_wind_speed_ms'],drawmedian=True,drawpercent=True)
    ax.set_ylabel('True_wind_speed (m/s)')
    ax.set_xlabel('Date')
    ax.grid(True)
    ax.legend()
    fig.savefig('./0.png',dpi=300,bbox_inches='tight')
    lg.info('Plot created: '+'./0.png')

    plt.close('all')
    fig=plt.figure(figsize=(15.0,5.0))
    ax=fig.add_subplot(1, 1,1)
    ax.plot(df['True_wind_dir_TO'],label='True_wind_dir_TO',color='r')
    ax.set_ylabel('True_wind_dir_TO ')
    tseries=df['True_wind_dir_TO']
    ax.vlines(x=tseries.index[350000], ymin=np.min(tseries), ymax=np.max(tseries),linewidth=6, color='green', zorder=1,alpha=0.5)
    ax.plot(tseries.rolling(window=500,center=True,win_type='hamming').mean(),label='rolling 500 hamming tsteps',color='k')
    ax.set_xlabel('Date')
    ax.grid(True)
    ax.legend()
    fig.savefig('./0b.png',dpi=300,bbox_inches='tight')
    lg.info('Plot created: '+'./0b.png')

    plt.close('all')
    fig=plt.figure()

    ws = np.random.random(500) * 6
    wd = np.random.random(500) * 360
    ax = WindroseAxes.from_ax(fig=fig)

    # Create wind speed and direction variables
    ax.bar(df['True_wind_dir_TO'], df['True_wind_speed_ms'], normed=True, opening=0.8, edgecolor='white')
    ax.set_legend(loc=2)
    fig.savefig('./1.png',dpi=300,bbox_inches='tight')
    lg.info('Plot created: '+'./1.png')



    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
