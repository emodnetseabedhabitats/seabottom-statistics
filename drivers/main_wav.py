#!/opt/anaconda3/envs/bke/bin/python
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import numpy as np
import os
import sys
sys.path.append(r'/Users/uranus/workarea/bke/seabottom2')
from bkenergy.mkbathy.ibathy import Gridinfo, Sbathy
from bkenergy.waves.iwke import Wbenergy
    
if __name__ == '__main__':

#-----------------------------------------------------------------
#/root ----|-----------------------------|
#       /data                      /seabottom2
#          |                            |
#     /cmems /emodnet /outs          /bkenergy -------->  /drivers
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

    #cmems waves
    cdata=root+'/data'+'/cmems'+'/spa_waves/'
    fileList = os.listdir(cdata)
    nc=[]
    for file in fileList:
        nc.append(file)

    # read target area from the first file
    filename = cdata + nc[0]  # waves
    
    g = Sbathy(filename, 0, *nb)
    #Also possible
    #f=Sbathy.from_file("user_bat.nc") #instance from user's bathymetry
    #output bathymetry
    g.to_file(out+'/user_bat_ibi.nc')
    #evaluate bke from currents using interpolated bathymetry
    h=Wbenergy(g.gmed,cdata,*nc)
    # general instance h=Sbenergy(g.gmed,cdata,*nc,density=1035,percent=90)
    # default values are density=1036.,percent=90.
        #output bathymetry and bke
    h.wbke_to_file(out+'/wbke_spain_out.nc')

