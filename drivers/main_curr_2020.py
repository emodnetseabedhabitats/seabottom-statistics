#!/opt/anaconda3/envs/bke/bin/python
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import numpy as np
import os
import sys
sys.path.append(r'/media/alessandro/sda/script/seabottom3')
from bkenergy.mkbathy.ibathy2020 import Gridinfo, Sbathy
from bkenergy.currents.icke import Sbenergy
    
if __name__ == '__main__':

#-----------------------------------------------------------------
#/root ----|---------------------------|
#       /data                      /sources
#          |                            |
#     /cmems /emodnet /outs          /bkenergy           /drivers
#        |                                 |
#       /ibi  /med   /backsea ...         mkbathy -->ibathy
#                                        currents --> icke
#                                          tands  --> its
#                                          waves  --> iwke
#-----------------------------------------------------------------
	root = '/media/alessandro/sda/script/seabottom3'
	bdata = root + '/data' + '/2020'
	out = root + '/data' + '/outs'
	
	# emodnet bathymetry
	bds = next(os.walk(bdata), (None, None, []))[2]  # [] if no file
	nb = [bdata + "/" + bd for bd in sorted(bds)]
	
	# cmems currents ibi
	cdata = root + '/data' + '/cmems' + '/ibi/curr/'
	fileList = os.listdir(cdata)
	nc = []
	for file in fileList:
		nc.append(file)
	'''
	# cmems bottom temperature and salinity
	tsdata = root + '/data' + '/cmems' + '/ts/'
	fileLists = os.listdir(tsdata)
	nts = []
	for file in fileLists:
		nts.append(file)
	'''
	# read target area from the first file
	filename = cdata + nc[0]  # currents
	# filename1=nb[2]
	#filename_ts = tsdata + nts[0]  # temperature and salinity
	#generate bathymetry
	g = Sbathy(filename, 0, *nb)
	#Also possible
	#f=Sbathy.from_file("user_bat.nc") #instance from user's bathymetry
	#output bathymetry
	g.to_file(out+'/bat_ibi_curr_2020.nc')
	#evaluate bke from currents using interpolated bathymetry
	h=Sbenergy(g.gmed,cdata,*nc)
	# general instance h=Sbenergy(g.gmed,cdata,*nc,density=1035,percent=90)
	# default values are density=1036.,percent=90.
		#output bathymetry and bke
	h.bke_to_file(out+'/bke_ibi_20202_curr.nc')

