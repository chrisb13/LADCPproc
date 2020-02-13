#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   www:     christopherbull.com.au
#   Date created: Sun, 26 Jan 2020 21:47:29
#

"""
Quick script to plot up the 1 Hz CTD files, just go get a sense of what each cast looked like.
"""
from cb2logger import *
import glob
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
import datetime
import numpy as np
import shareme as sm

def gtime(t,delta):
    tinit=datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5])
    tseries=[]
    for delt in delta:
        #seconds offset - could be a bug - ask Mike about time zero
        tseries.append(pd.to_datetime(tinit+datetime.timedelta(0,delt[0])))

    return tseries


if __name__ == "__main__": 
    LogStart('',fout=False)
    ifol='/home/chris/VBoxSHARED/Chris_LADCP/processing/ctd/*1Hz_ascii'
    ifol='/home/chris/VBoxSHARED/Chris_LADCP/processing/ctd/*1hz'
    pfol='./plots/'
    sm.mkdir(pfol)

    ifiles=sorted(glob.glob(ifol))
    assert(ifiles!=[]),"glob didn't find anything!"
    # df=pd.read_csv(ifiles[0])
    for f in ifiles:
        mat = scipy.io.loadmat(f)

        #the problematic ones..
        bads=[13, 15, 16, 23, 41, 42, 50, 60, 61, 63]
        for bad in bads:
            if str(bad) in os.path.basename(f):
                co='red'
                # print 'yay'
                break
            else:
                co=None
        t=gtime(mat['gtime'][0],mat['time_elapsed'])

        df=pd.DataFrame({'press':mat['press'].flatten(),'salin':mat['salin'].flatten(),'temp':mat['temp'].flatten()},index=t)

        plt.close('all')
        fig=plt.figure(figsize=(12.0,4.0))
        ax=fig.add_subplot(1, 3,1)
        ax.plot(df['temp'],df['press'],label='temp')
        plt.gca().invert_yaxis()
        ax.set_title('temp')
        ax.set_ylabel('Pressure')
        ax.set_xlabel('Temperature')
        ax.grid(True)

        ax=fig.add_subplot(1, 3,2)
        ax.plot(df['salin'],df['press'],label='temp')
        plt.gca().invert_yaxis()
        ax.set_title('salinity')
        ax.set_xlabel('Salinity')
        plt.setp(ax.get_yticklabels(),visible=False)
        ax.grid(True)

        ax=fig.add_subplot(1, 3,3)
        ax.plot(df.index,df['press'],label='temp',color=co)
        ax.plot(df.index,[np.max(df['press'])*.9]*len(df['press']),color='k',alpha=0.5,label='90% max press')
        plt.gca().invert_yaxis()
        ax.set_xlabel('Time')
        plt.setp(ax.get_yticklabels(),visible=False)
        ax.grid(True)

        efile=pfol+os.path.basename(f)+'.png'
        fig.savefig(efile,dpi=300,bbox_inches='tight')
        lg.info('plot created: '+efile)

    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
