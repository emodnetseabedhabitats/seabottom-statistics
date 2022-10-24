#    Copyright 2022 Roberto Inghilesi, Alessandro Mercatini

#    This file is part of seabottom 4.1

#    bchem.py is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    icke.py is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.

#    You should have received a copy of the GNU Lesser General Public License
#    along with bchem.py.  If not, see <http://www.gnu.org/licenses/>.

from netCDF4 import Dataset
import numpy as np


class Sbchem(object):
    """The Class Sbchem is written to provide an evaluation of bottom temperature, salinity and o2 at the seabottom by means of self-contained objects.
    it needs a bathymetry interpolated on the grid and montly/daily data
    from oceanographic models i.e. for IBI files like for example
    dataset-ibi-analysis-forecast-phys-005-001-monthly_1544612310558.nc
    all inputs must be netcdf files. A class (Sbathy) is provided in order to interpolate Emodnet bathymetry on the specific model grid selected.
    usage:
    1st step
    ebaluate the bathymetry by using
    gts=Sbathy(filename_ts,0,*nb
    or gts=Sbathy.fromfile('namefile.nc')
    (see specific documentation)
    2nd step
    evaluate nth percentile of bottom biochemical or temperature
    h=Sbchem(g.gmed,cdata,*nc, parameter='salinity',percent=90)
    (gts.gmed must be generated as an instance gts of class Sbathy).
    percent=xx indicates the xx nth percentile to be evalated, and parameter=xxxx let the user to insert the name of the desired
    parameter to evaluate. Allowed parameters are 'salinity' , 'btemperature', and'oxygen'.
    
    3nd step
    output nth percentile of seabottom energy
    hts.bts_to_file('bchem_ibi_sal_out.nc')
    where 'bchem_ibi_sal_out' is just the name of the output netcdf file.
    Output variables are latitude, longitude, bathymetry, xx percentile field of parameter xxxx, average and std devuation of xxxx.
    Technical reference can be found in the Annex 1 "Compiling oceanographic layers" to the EUSeaMap 2019 - Technical Report
    """

    @property
    def gmed(self):
        return self._gmed

    @property
    def parameter(self):
        return self._parameter

    @property
    def bchm(self):
        return self._bchm
    @property
    def avg_bchm(self):
        return self._avg_bchm
            
    @property
    def std_bchm(self):
        return self._std_bchm

    @property
    def percent(self):
        return self._percent
        # @property

    def __init__(self, gmed, cdata, *args, parameter='salinity', percent=90):
        #parameter= 'salinity' (so),'oxygen'(o2),'Chlorophyll'(chl),btemperature (bottomT)
        #'ammonium'(nh4),'phosphate'(po4),'nitrate'(no3)
        self._gmed = gmed
        self._parameter = parameter
        self._percent = percent
        ftimes = gmed.ntimes
        ndix=len(args)*ftimes+5*ftimes
        #su_t = np.zeros([gmed.nlats, gmed.nlons, len(args) * ftimes], dtype=float)
        su_s = np.zeros([gmed.nlats, gmed.nlons, ndix], dtype=float)
        avgbchm = np.zeros([gmed.nlats, gmed.nlons], dtype=float)
        stdbchm = np.zeros([gmed.nlats, gmed.nlons], dtype=float)
        #tsbot = np.zeros([gmed.nlats, gmed.nlons, gmed.ntimes], dtype=float)

        n1 = 0
        n2=-1
        ltimes=0
        for l in args:  # [:2]:
            print(cdata + l)
            tsbot = self.btands(cdata + l)
            #print('range init:',np.min(tsbot),np.max(tsbot))
            ldims=np.shape(tsbot)
            n1 = n2+1
            ltimes=ldims[2]
            n2 = n1+ltimes-1
            print(ltimes,n1,n2)
            su_s[:, :, n1:n2+1] = tsbot[:, :, :]
            print('range bchm:',np.min(tsbot),np.max(tsbot))
            
        
        #self._bchm = self.process(gmed.nlats, gmed.nlons, su_s, percent)
        
        msu_s=np.where(su_s>0.,su_s,np.nan)
        avgbchm=np.nanmean(msu_s,axis=2)
        stdbchm=np.nanstd(msu_s,axis=2)
        prcdbchm=np.nanpercentile(msu_s,int(percent),axis=2)
        self._std_bchm=stdbchm
        self._avg_bchm=avgbchm
        self._bchm =prcdbchm
        print('range bchm:',np.min(avgbchm),np.max(avgbchm))
        
    def btands(self, filen):
        """The method btands is used to evaluate the salinity on the model level closest to the
        sea bottom at each point of the grid. No BL model is applied.
        It is called by the constructor of the class Btands.
        """
        mv_lst = ['missing_value', '_FillValue']
        #t_lst = ['bottomT', 'temperature']
        if self.parameter.casefold()=='salinity':
            s_lst=['so','salinity']
        elif self.parameter.casefold()=='oxygen':
            s_lst=['o2','oxygen']
        elif self.parameter.casefold()=='chlorophyll':
            s_lst=['chl','clorophyll']
        elif self.parameter.casefold()=='nitrate':
            s_lst=['no3','nitrate']
        elif self.parameter.casefold()=='ammonium':
            s_lst=['nh4','ammonium']
        elif self.parameter.casefold()=='phosphate':
            s_lst=['po4','phosphate']
        elif self.parameter.casefold()=='btemperature':
            t_lst=['bottomT','temperature']
        else:
            print('No suitable parameter found in data file for ',self.parameter)
        #rho = density
        fh = Dataset(filen, mode='r', format="NETCDF4")
        ltimes = int(fh.dimensions['time'].size)
        res = np.zeros([self.gmed.nlats, \
                        self.gmed.nlons,ltimes], dtype=float)
        if self.parameter.casefold()=='btemperature':
            self._gmed.nlevs=1
            for l in t_lst:
                if l in fh.variables.keys():
                    bte = np.array(fh[l])
                    if 'scale_factor' in dir(fh[l]):
                        tsf = np.float64(fh[l].scale_factor)
                    else:
                        tsf=1.
                    if 'add_offset' in dir(fh[l]):
                        tof = np.float64(fh[l].add_offset)
                    else:
                        tof=0.0
                for t in mv_lst:
                    if t in fh[l].__dict__.keys():
                        mvt = np.float64(fh[l].__dict__[t])
                break
            btd = np.where(bte == mvt, bte, bte * tsf + tof)
            print('range btd:',np.min(btd),np.max(btd))
        else:
            nlevs = int(fh.dimensions['depth'].size)
            self._gmed.nlevs=nlevs
            lev = np.array(fh["depth"])
            for l in s_lst:
                print(l)
                if l in fh.variables.keys():
                    print(l,'*')
                    bsa = np.array(fh[l])
                    if 'scale_factor' in dir(fh[l]):
                        ssf = np.float64(fh[l].scale_factor)
                    else:
                        ssf=1.
                    if 'add_offset' in dir(fh[l]):
                        sof = np.float64(fh[l].add_offset)
                    else:
                        sof=0.0
                    for t in mv_lst:
                        if t in fh[l].__dict__.keys():
                            mvs = np.float64(fh[l].__dict__[t])
                    break
            
            bsd = np.where(bsa == mvs, bsa, bsa * ssf + sof)
            print('mvs=',mvs,ssf,sof)
            print('range bsd:',np.min(bsd),np.max(bsd))
        #inc = 0

        for ilat in range(self.gmed.nlats):
            for ilon in range(self.gmed.nlons):
                for it in range(self.gmed.ntimes):
                    if self.parameter.casefold()=='btemperature':
                         res[ilat, ilon, it] = btd[it, ilat, ilon]
                    else:
                        #inc += 1
                        sfloor = self.gmed.z[ilat, ilon]
                        if (sfloor == 0.):
                             continue
                        a = bsd[it, :, ilat, ilon]
                        il=0
                        for ck in a:
                            if ck == mvs:
                                ckp1 = il
                                break
                            il += 1
                        if ckp1 == 0:
                            continue
                        ckp = ckp1 - 1
                        # IBI doesn't work with sf and offset use bsa
                        res[ilat, ilon, it] = bsd[it, ckp, ilat, ilon]
        print('range res:',np.min(res),np.max(res))
        return res

    def process(self, nlats, nlons, su_ke, perc):
        """The method process is used to evaluate the required
        statistics (nth percentile) at each point of the grid.
        It is called by the constructor of the class Stands"""
        ibke = np.zeros([nlats, nlons], dtype=float)
        nt = np.shape(su_ke)
        ntp = nt[2]
        # perc=self.percent/100.
        lep =[]# np.zeros(ntp, dtype=float)
        print('ntp= ', ntp)
        
        #pos = int(ntp * perc / 100.) + 1
        #print('range process:',np.min(su_ke),np.max(su_ke))
        for ilat in range(nlats):
            for ilon in range(nlons):
                lep = np.copy(su_ke[ilat, ilon, :])
                lep[::-1].sort()
                #lep[:] = su_ke[ilat, ilon, :]
                # if ilon==100 and ilat==100:
                #     print(lep)
                ck=np.argmin(lep)
                prc=int(ck*(1.-perc/100.))+1
                #print(ntp,ck,prc)
                if ilon==100 and ilat==100:
                    print(ntp,ck,prc)
                    print(lep)
                ibke[ilat, ilon] = lep[prc]
        return ibke

    def bts_to_file(self, fl):
        """The method bts_to_file is used to output the ncdf4 file
        of statistics of bottom salinity"""
        rot = Dataset(fl, "w", format="NETCDF4")
        lon = rot.createDimension("lon", self.gmed.nlons)
        lat = rot.createDimension("lat", self.gmed.nlats)
        latitudes = rot.createVariable("lat", "f4", ("lat",))
        longitudes = rot.createVariable("lon", "f4", ("lon",))
        latitudes[:] = self.gmed.lats
        longitudes[:] = self.gmed.lons
        bat = rot.createVariable("bat", "f4", ("lon", "lat",))
        bat[:] = self.gmed.z.transpose()
        bchm_perc = rot.createVariable("bchm_perc", "f4", ("lon", "lat",))
        bchm_perc[:] = self.bchm.transpose()
        bchm_avg = rot.createVariable("bchm_avg", "f4", ("lon", "lat",))
        bchm_avg[:] = self.avg_bchm.transpose()
        bchm_std = rot.createVariable("bchm_std", "f4", ("lon", "lat",))
        bchm_std[:] = self.std_bchm.transpose()
        #btemp = rot.createVariable("btemp", "f4", ("lon", "lat",))
        #btemp[:] = self.btemp.transpose()
        rot.close()
