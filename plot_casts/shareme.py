#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   www:     christopherbull.com.au
#   Date created: Wed, 22 Jan 2020 17:19:00
#   Machine created on: SB2Vbox
#

"""
Python module to contain useful shared parameters
"""

from cb2logger import *
import collections
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.pyplot as plt

cast_type=collections.OrderedDict()


cast_type['001']=['test','']
# sheldon cove
cast_type['002']=['icebergs','SC0'] # no ladcp data
cast_type['003']=['racetrax','SC0']
cast_type['004']=['racetrax','SC0']
cast_type['005']=['icebergs','SC1']
cast_type['006']=['racetrax','SC1']
cast_type['007']=['racetrax','SC1']
cast_type['008']=['icebergs','SC2']
cast_type['009']=['icebergs','SCA']
cast_type['010']=['racetrax','SCA']
cast_type['011']=['racetrax','SCA']
cast_type['012']=['icebergs','SCB']
cast_type['013']=['icebergs','SC6']
cast_type['014']=['icebergs','IF1']
cast_type['015']=['icebergs','IF2']
cast_type['016']=['icebergs','IF3']
cast_type['017']=['icebergs','SCC']
cast_type['018']=['racetrax','SCC']
cast_type['019']=['racetrax','SCC']
cast_type['020']=['icebergs','SCD']
cast_type['021']=['racetrax','SCE']
cast_type['022']=['racetrax','SCE']
cast_type['023']=['icebergs','SCE']

# borgen bay
cast_type['024']=['icebergs','BB0']
cast_type['025']=['racetrax','BB0']
cast_type['026']=['racetrax','BB0']
cast_type['027']=['icebergs','BBX']
cast_type['028']=['icebergs','BBA']
cast_type['029']=['racetrax','BBA']
cast_type['030']=['racetrax','BBA']
cast_type['031']=['icebergs','BBC']
cast_type['032']=['racetrax','BBC']
cast_type['033']=['racetrax','BBC']
cast_type['034']=['icebergs','BB-ICE3']
cast_type['035']=['icebergs','BB-ICE4']
cast_type['036']=['icebergs','BB-ICE2']
cast_type['037']=['icebergs','BB-ICE1']
cast_type['038']=['icebergs','BBD']
cast_type['039']=['racetrax','BBD']
cast_type['040']=['racetrax','BBD']
cast_type['041']=['icebergs','BBE']
cast_type['042']=['racetrax','BBE']
cast_type['043']=['racetrax','BBE']
cast_type['044']=['icebergs','BB1']
cast_type['045']=['EK80cal','']
cast_type['046']=['icebergs','BBB']

# marian cove
cast_type['047']=['icebergs','MC1']
cast_type['048']=['racetrax','MC1']
cast_type['049']=['icebergs','MC0']
cast_type['050']=['icebergs','MC2']
cast_type['051']=['icebergs','MCA']
cast_type['052']=['racetrax','MCA']
cast_type['053']=['icebergs','MCB']
cast_type['054']=['racetrax','MCB']
cast_type['055']=['icebergs','MCC']
cast_type['056']=['racetrax','MCC']
cast_type['057']=['icebergs','MC4'] #4 was possibly a typo?
cast_type['058']=['icebergs','MCD']
cast_type['059']=['racetrax','MCD']
cast_type['060']=['icebergs','MC-ICE1']
cast_type['061']=['icebergs','MC-ICE2']
cast_type['062']=['racetrax','MCE']
cast_type['063']=['icebergs','MCE']

def inset_title_box(ax,title,bwidth="20%",location=1):
    """
    Function that puts title of subplot in a box
    
    :ax:    Name of matplotlib axis to add inset title text box too
    :title: 'string to put inside text box'
    :returns: @todo
    """

    axins = inset_axes(ax,
                       width=bwidth, # width = 30% of parent_bbox
                       height=.30, # height : 1 inch
                       loc=location)

    plt.setp(axins.get_xticklabels(), visible=False)
    plt.setp(axins.get_yticklabels(), visible=False)
    axins.set_xticks([])
    axins.set_yticks([])

    axins.text(0.5,0.3,title,
            horizontalalignment='center',
            transform=axins.transAxes,size=10)

def mkdir(p):
    """make directory of path that is passed"""
    try:
       os.makedirs(p)
       lg.info("output folder: "+p+ " does not exist, we will make one.")
    except OSError as exc: # Python >2.5
       import errno
       if exc.errno == errno.EEXIST and os.path.isdir(p):
          pass
       else: raise

if __name__ == "__main__": 
    print('not supposed to be run like this!')
    # LogStart('',fout=False)
    # #put useful code here!

    # lg.info('')
    # localtime = time.asctime( time.localtime(time.time()) )
    # lg.info("Local current time : "+ str(localtime))
    # lg.info('SCRIPT ended')
