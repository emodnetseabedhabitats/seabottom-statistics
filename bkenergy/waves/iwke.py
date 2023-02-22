#    Copyright 2022 Roberto Inghilesi, Alessandro Mercatini

#    This file is part of seabottom 4.1

#    iwke.py is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    icke.py is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.

#    You should have received a copy of the GNU Lesser General Public License
#    along with iwke.py.  If not, see <http://www.gnu.org/licenses/>.

from netCDF4 import Dataset
import numpy as np

class Wbenergy(object):
    """
    The Class is written to provide an evaluation of current
    kinetic energy at the seabottom by means of
    self-contained objects. it needs a bathymetry
    interpolated on the grid and montly data from
    oceanographic models i.e. for IBI files like
    dataset-ibi-analysis-forecast-phys-005-001-monthly_xxx.nc
    all inputs must be netcdf files. A class (Sbathy) is
    provided in order to interpolate Emodnet bathymetry on
    the specific model grid selected. usage->
    1st step evaluate the bathymetry by using
    g=Sbathy(filename,0,*nb)
    or g=Sbathy.fromfile('namefile.nc')
    (see specific documentation)
    2nd step
    evaluate nth percentile of bottom kinetic energy due to
    waves as h=Wbenergy(g.gmed,wdata,*nw)
    (g.gmed must be generated as an instance of class
    Sbathy), or
    h=Wbenergy(g.gmed,cdata,*nw,density=1035.,percent=90)
    where percent=xx indicates the xx nth percentile to be
    evalated, and density=1035. lets the user to insert a
    value of mean density of the area to be used in the
    evaluation of the bottom enegy.
    3nd step
    output nth percentile of seabottom energy
    h.wbke_to_file('bke_ibi_out.nc')
    Output variables are latitude, longitude, bathymetry,
    xx percentile field of parameter xxxx, average and std devuation of xxxx.
    Technical reference can be found in the Annex 1 "Compiling
    oceanographic layers" to the EUSeaMap 2019 - Technical Report
    """

    @property
    def gmed(self):
        return self._gmed

    @property
    def wbke(self):
        return self._wbke

    @property
    def avg_wbke(self):
        return self._avg_wbke
            
    @property
    def std_wbke(self):
        return self._std_wbke
        
    @property
    def percent(self):
        return self._percent

    @property
    def density(self):
        return self._density

    def __init__(self, gmed, cdata, *args, density=1036., percent=90):
        self._gmed = gmed
        self._density = density
        self._percent = percent
        ftimes = gmed.ntimes
        ndix=len(args)*ftimes+5*ftimes
        su_ke=np.zeros([gmed.nlats,gmed.nlons,ndix],dtype=float)
        avgbke = np.zeros([gmed.nlats, gmed.nlons], dtype=float)
        stdbke = np.zeros([gmed.nlats, gmed.nlons], dtype=float)
        
        n1=0
        n2=-1
        ltimes=0
        for l in args:  # [:2]:
            print(cdata + l)
            tbke = self.bken(cdata + l, density)
            ldims=np.shape(tbke)
            n1 = n2+1
            ltimes=ldims[2]
            n2 = n1+ltimes-1
            print(ltimes,n1,n2)
            su_ke[:,:,n1:n2+1]=tbke[:,:,:]
            # print('res', n,n1,n2,tbke[100,100,:])
            # print('shp', np.shape(tbke),np.shape(su_ke[:,:,n1:n2]))
            #su_ke[:, :, n1:n2] = tbke[:, :, :]
            #n += 1
        #self._bke = self.process(gmed.nlats, gmed.nlons, su_ke, percent)
        #avgbke=np.average(su_ke,axis=2)
        #stdbke=np.std(su_ke,axis=2)
        msu_ke=np.where(su_ke>0.,su_ke,np.nan)
        avgbke=np.nanmean(msu_ke,axis=2)
        stdbke=np.nanstd(msu_ke,axis=2)
        prcdbke=np.nanpercentile(msu_ke,int(percent),axis=2)
        self._std_wbke=stdbke
        self._avg_wbke=avgbke
        self._wbke=prcdbke
        
    def __str__(self):
        return "([{} \n {}])".format(self.density, self.gmed)

    def process(self, nlats, nlons, su_ke, perc):
        """The method process is used to evaluate the
        required statistics (nth percentile) at each point
        of the grid"""
        ibke = np.zeros([nlats, nlons], dtype=float)
        nt = np.shape(su_ke)
        ntp = nt[2]
        # perc=self.percent/100.
        lep = []#np.zeros(ntp, dtype=float)
        print('ntp= ', ntp)
        #prc = int(ntp * perc / 100.) + 1
        for ilat in range(nlats):
            for ilon in range(nlons):
                #lep[:] = su_ke[ilat, ilon, :]
                lep = np.copy(su_ke[ilat, ilon, :])
                lep[::-1].sort()
                ck=np.argmin(lep)
                #prc=int(ck*perc/100.)+1
                prc=int(ck*(1.-perc/100.))+1
                if ilon==100 and ilat==100:
                    print(ntp,ck,prc)
                    print(lep)
                ibke[ilat,ilon]=lep[prc]
                
        return ibke

    def bken(self, filen, density):
        """The method bken is used to evaluate the kinetic
        energy close (zfix=1m) to the sea bottom at each
        point of the grid. It is called by the constructor
        of the class Benergy. The simple boundary layer
        model for waves is based on a simplified method
        proposed by Soulby (1987) and described in the Annex
        1- Compiling oceanographic layers to the EUSeaMap
        2019 - Technical Report
        """
        tp_lst=['VTPK',]
        hs_lst=['VHM0',]
        mv_lst=['missing_value','_FillValue']
        #res = np.zeros([self.gmed.nlats, \
        #self.gmed.nlons, self.gmed.ntimes], dtype=float)
        rho = density
        fh = Dataset(filen, mode='r', format="NETCDF4")
        # nlevs=int(fh.dimensions['depth'].size)
        ltimes = int(fh.dimensions['time'].size)
        res=np.zeros([self.gmed.nlats,self.gmed.nlons,ltimes],dtype=float)
        # lev=np.array(fh["depth"])
        for l in tp_lst:
            if l in fh.variables.keys():
                tp = np.array(fh[l])
                if 'scale_factor' in dir(fh[l]):
                    tpsf = np.float64(fh[l].scale_factor)
                else:
                    tpsf=1.
                if 'add_offset' in dir(fh[l]):
                    tpof = np.float64(fh[l].add_offset)
                else:
                    tpof=0.0
                for t in mv_lst:
                    if t in fh[l].__dict__.keys():
                        mvtp = np.float64(fh[l].__dict__[t])
                break
        for l in hs_lst:
            if l in fh.variables.keys():
                hs = np.array(fh[l])
                if 'scale_factor' in dir(fh[l]):
                    hssf = np.float64(fh[l].scale_factor)
                else:
                    hssf=1.
                if 'add_offset' in dir(fh[l]):
                    hsof = np.float64(fh[l].add_offset)
                else:
                    hsof=0.0
                #hssf = np.float64(fh[l].scale_factor)
                #hsof = np.float64(fh[l].add_offset)
                for t in mv_lst:
                    if t in fh[l].__dict__.keys():
                        mvhs = np.float64(fh[l].__dict__[t])
                break
        # hs=np.where(hs==mvhs,hs,hs*hssf+hsof)
        # tp=np.where(tp==mvtp,tp,tp*tpsf+tpof)

        # inc=0
        for ilat in range(self.gmed.nlats):
            for ilon in range(self.gmed.nlons):
                # for it in range(self.gmed.ntimes):
                #    inc+=1
                dd = self.gmed.z[ilat, ilon]
                if (dd == 0.):
                    continue
                
                hm0 = hs[:, ilat, ilon]
                tpz = tp[:, ilat, ilon]
                hm0 = np.where(hm0 == mvhs, 0, hm0)
                tpz = np.where(tpz < 0.01, 0.01, tpz)
                k = np.shape(hm0)
                dep = np.zeros(k)
                dep = dep + abs(dd)
                dep = np.where(hm0 > 0.8 * dep, 10., dep)
                dep = np.where(dep < 4., 4., dep)
                # if hm0==mvhs or tpz == mvtp:
                #    continue

                # if hm0 > 0.8*dd:
                #    dd=10.

                kub = self.burms(hm0, tpz, dep)
                kk = 0.5 * rho * kub ** 2
                # if kk>2000.:
                #    print('check kk', kk, ilon,ilat,it,hm0,tpz)
                res[ilat, ilon, :] = kk[:]  # limiter: every time block must conform to the first
        return res

    def burms(self, hs, tz, hb):
        """ implementation of Soulsby (1987) simplified algorithm"""
        grav = 9.81
        tn = np.sqrt(hb / grav)
        t = tn / tz
        a = (6500. + (0.56 + 15.54 * t) ** 6.) ** (1. / 6.)
        b = (1. + a * t ** 2.) ** 3.
        # if b<0.01:
        #    print('check b',b)
        # if hb<0.01:
        #    print('check hb',hb)
        urms = hs / (4. * tn * b)
        return urms

    def wbke_to_file(self, fl):
        """The method bke_to_file is used to output the
        ncdf4 file of statistics of bottom kinetic energy"""
        rot = Dataset(fl, "w", format="NETCDF4")
        lon = rot.createDimension("lon", self.gmed.nlons)
        lat = rot.createDimension("lat", self.gmed.nlats)
        latitudes = rot.createVariable("lat", "f4", ("lat",))
        longitudes = rot.createVariable("lon", "f4", ("lon",))
        latitudes[:] = self.gmed.lats
        longitudes[:] = self.gmed.lons
        bat = rot.createVariable("bat", "f4", ("lon", "lat",))
        bat[:] = self.gmed.z.transpose()
        wbke_perc = rot.createVariable("wbke_perc", "f4", ("lon", "lat",))
        wbke_perc[:] = self.wbke.transpose()
        wbke_avg = rot.createVariable("wbke_avg", "f4", ("lon", "lat",))
        wbke_avg[:] = self.avg_wbke.transpose()
        wbke_std = rot.createVariable("wbke_std", "f4", ("lon", "lat",))
        wbke_std[:] = self.std_wbke.transpose()
        
        rot.close()
