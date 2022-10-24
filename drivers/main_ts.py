#!/opt/anaconda3/envs/bke/bin/python
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import numpy as np
import os
import sys
sys.path.append(r'/Users/uranus/workarea/bke/seabottom2')
from bkenergy.mkbathy.ibathy import Gridinfo, Sbathy
from bkenergy.mktands.itands import Stands
    
if __name__ == '__main__':

#-----------------------------------------------------------------
#/root ----|----------------------------|
#       /data                      /seabottom2
#          |                            |
#     /cmems /emodnet /outs          /bkenergy           /drivers
#        |                                 |
#       /ibi  /med   /backsea ...         mkbathy -->ibathy
#                                        currents --> icke
#                                          tands  --> its
#                                          waves  --> iwke
#-----------------------------------------------------------------
    root = '/Users/uranus/workarea/bke'
    bdata = root + '/data' + '/emodnet'
    out = root + '/data' + '/outs'

    # emodnet bathymetry
    bds = ["A2", "A3", "A4", "B2", "B3", "B4", "C2", "C3", "C4", "D3", "D4"]
    nb = [bdata + "/" + bd + ".mnt" for bd in bds]

    # cmems bottom temperature and salinity
    tsdata = root + '/data' + '/cmems' + '/ibi_ts/'
    fileLists = os.listdir(tsdata)
    nts = []
    for file in fileLists:
        nts.append(file)
    # read target area from the first file
    filename = tsdata + nts[0]  # temperature and salinity
    #generate bathymetry
    g = Sbathy(filename, 0, *nb)
    #Also possible
    #f=Sbathy.from_file("user_bat.nc") #instance from user's bathymetry
    #output bathymetry
    g.to_file(out+'/user_bat_ibi.nc')
    #evaluate bke from currents using interpolated bathymetry
    h=Stands(g.gmed,tsdata,*nts)
    # general instance h=Stands(g.gmed,cdata,*nc,density=1035,percent=90)
    # default values are density=1036.,percent=50.
        #output bathymetry and bke
    h.bts_to_file(out+'/bts_ibi_out.nc')

