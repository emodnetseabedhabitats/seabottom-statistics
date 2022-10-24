#    Copyright 2022 Roberto Inghilesi, Alessandro Mercatini

#    This file is part of seabottom 4.1

#    ibathy.py is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    icke.py is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.

#    You should have received a copy of the GNU Lesser General Public License
#    along with icke.py.  If not, see <http://www.gnu.org/licenses/>.

from netCDF4 import Dataset
import numpy as np

class Sbenergy(object):
    """
    The Class is written to provide an evaluation of current kinetic energy at the
    seabottom by means of self-contained objects.
    it needs a bathymetry interpolated on the grid and montly data
    from oceanographic models i.e. for IBI files like
    dataset-ibi-analysis-forecast-phys-005-001-monthly_xxx.nc
    all inputs must be netcdf files. A class (Sbathy) is provided in order
    to interpolate Emodnet bathymetry on the specific model grid selected.
    usage:
    1st step
    evaluate the bathymetry by using
    g=Sbathy(filename,0,*nb)
    or g=Sbathy.fromfile('namefile.nc')
    (see specific documentation)
    2nd step
    evaluate nth percentile of bottom energy
    h=Sbenergy(g.gmed,cdata,*nc,) #instance from emodnet 2016 bathymetry
    (g.gmed must be generated as an instance of class Sbathy),
    or
    h=Sbenergy(g.gmed,cdata,*nc,density=1035.,percent=90)
    where percent=xx indicates the xx nth percentile to be evalated, and
    density=1035. let the user to insert a value of mean density of the area to be used in the evaluation of the bottom enegy.
    3nd step
    output nth percentile of seabottom energy
    h.bke_to_file('bke_ibi_out.nc')
    Output variables are latitude, longitude, bathymetry, xx percentile field of bke,average and std deviation of bke.
    Technical reference can be found in the Annex 1 "Compiling oceanographic layers" to the EUSeaMap 2019 - Technical Report
    """

    @property
    def gmed(self):
        return self._gmed

    @property
    def bke(self):
        return self._bke
        
    @property
    def avg_bke(self):
        return self._avg_bke
            
    @property
    def std_bke(self):
        return self._std_bke
        
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
        #su_ke = np.zeros([gmed.nlats, gmed.nlons, len(args) * ftimes], dtype=float)
        #tbke = np.zeros([gmed.nlats, gmed.nlons, gmed.ntimes], dtype=float)
        #n = 0
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
        #self._bke = self.process(gmed.nlats, gmed.nlons, su_ke, percent)
        #avgbke=np.average(su_ke,axis=2)
        #stdbke=np.std(su_ke,axis=2)
        #self._std_bke=stdbke
        #self._avg_bke=avgbke
        msu_ke=np.where(su_ke>0.,su_ke,np.nan)
        avgbke=np.nanmean(msu_ke,axis=2)
        stdbke=np.nanstd(msu_ke,axis=2)
        prcdbke=np.nanpercentile(msu_ke,int(percent),axis=2)
        self._std_bke=stdbke
        self._avg_bke=avgbke
        self._bke =prcdbke

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
                lep = np.copy(su_ke[ilat, ilon, :])
                #lep[:] = su_ke[ilat, ilon, :]
                lep[::-1].sort()
                #if ilat==100 and ilon==100:
                #    print('su_ker',lep)
                ck=np.argmin(lep)
                prc=int(ck*(1.-perc/100.))+1
                #print(ntp,ck,prc)
                if ilon==100 and ilat==100:
                    print(ntp,ck,prc)
                    print(lep)
                #ll = np.sort(lep)
                # if ilon==100 and ilat==100:
                #     print(ll)
                ibke[ilat, ilon] = lep[prc]
        return ibke

    def bken(self, filen, density):
        """The method bken is used to evaluate the kinetic
        energy close (zfix=1m) to the sea bottom at each
        point of the grid. It is called by the constructor
        of the class Benergy. The simple boundary layer
        model is described in the Annex 1 "Compiling
        oceanographic layers" to the EUSeaMap 2019 -
        Technical Report"
        """
        uc_lst = ['uo', 'vozocrtx']
        vc_lst = ['uo', 'vomecrty']
        mv_lst = ['missing_value', '_FillValue']
        #res = np.zeros([self.gmed.nlats, \
        #    self.gmed.nlons, self.gmed.ntimes], dtype=float)
        om = 2 * np.pi / 24. / 3600.
        cdmin = 0.0025
        zob = 0.0035
        zfix = 1.0
        # rho=1036.
         
        rho = density
        fh = Dataset(filen, mode='r', format="NETCDF4")
        nlevs = int(fh.dimensions['depth'].size)
        self._gmed.nlevs=nlevs
        ltimes = int(fh.dimensions['time'].size)
        lev = np.array(fh["depth"])
        res=np.zeros([self.gmed.nlats,\
                      self.gmed.nlons,ltimes],dtype=float)
        
        for l in uc_lst:
            if l in fh.variables.keys():
                cu=np.array(fh[l])
                if 'scale_factor' in dir(fh[l]):
                    usf=np.float64(fh[l].scale_factor)
                else:
                    usf=1.
                if 'add_offset' in dir(fh[l]):
                    uof=np.float64(fh[l].add_offset)
                else:
                    uof=0.0
                for t in mv_lst:
                    if t in fh[l].__dict__.keys():
                        mvu=np.float64(fh[l].__dict__[t])
                break
        for l in vc_lst:
            if l in fh.variables.keys():
                cv=np.array(fh[l])
                if 'scale_factor' in dir(fh[l]):
                    vsf=np.float64(fh[l].scale_factor)
                else:
                    vsf=1.
                if 'add_offset' in dir(fh[l]):
                    vof=np.float64(fh[l].add_offset)
                else:
                    vof=0.0
                for t in mv_lst:
                    if t in fh[l].__dict__.keys():
                        mvv=np.float64(fh[l].__dict__[t])
                break
        cu=np.where(cu==mvu,cu,cu*usf+uof)
        cv=np.where(cv==mvv,cv,cv*vsf+vof)
        ##print('maxcu=',np.max(cu))
        # print('maxcv=',np.max(cv))
        ang = [np.pi / 180. * ilt for ilt in self.gmed.lats]
        inc = 0
        fco = 2. * om * np.sin(ang)
        for ilat in range(self.gmed.nlats):
            for ilon in range(self.gmed.nlons):
                for it in range(ltimes):
                    inc += 1
                    sfloor = self.gmed.z[ilat, ilon]
                    if (sfloor == 0.):
                        continue
                    a = cu[it, :, ilat, ilon]
                    b = cv[it, :, ilat, ilon]
                    # to be substituted/check immediately
                    il = 0
                    ckp1 = 0
                    for ck in a:
                        if ck == mvu:
                            ckp1 = il
                            break
                        il += 1

                    # ckp1=np.argmin(a)
                    il = 0
                    ckp2 = 0
                    for ck in b:
                        if ck == mvv:
                            ckp2 = il
                            break
                        il += 1

                    # ckp2=np.argmin(a)
                    if ckp1 == 0 or ckp2 == 0:
                        continue

                    if ckp1 == ckp2:
                        ckp = ckp1 - 1
                    else:
                        print('warning ckp', ckp1, ckp2)

                    u = [cu[it, ckp, ilat, ilon] ** 2 + \
                         cv[it, ckp, ilat, ilon] ** 2]
                    g = np.sqrt(u)

                    layer = sfloor + lev[ckp]

                    dh = lev[ckp + 1] - lev[ckp]
                    ustar = self.bleach(ckp, sfloor, g, dh, zob)
                    delta = ustar / fco[ilat]
                    zbt = 0.4 * delta
                    # if layer>200.:
                    #   print(ckp,lev[ckp],layer,ilat,ilon)

                    if abs(zbt) > abs(layer):
                        for j in range(ckp):
                            ckp = ckp - 1
                            u = [cu[it, ckp, ilat, ilon] ** 2 + \
                                 cv[it, ckp, ilat, ilon] ** 2]
                            g = np.sqrt(u)
                            layer = sfloor + lev[ckp]
                            dh = lev[ckp + 1] - lev[ckp]
                            ustar = self.bleach(ckp, sfloor, g, dh, zob)
                            delta = ustar / fco[ilat]
                            zbt = 0.4 * delta
                            if (abs(zbt) <= abs(layer)):
                                break
                    ubt = 2.5 * ustar * np.log(zfix / zob)
                    # ubt1=g+2.5*ustar*np.log(zbt/delta)
                    res[ilat, ilon, it] = 0.5 * ubt ** 2. * rho
        return res

    def bke_to_file(self, fl):
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
        bke_perc = rot.createVariable("bke_perc", "f4", ("lon", "lat",))
        bke_perc[:] = self.bke.transpose()
        bke_avg = rot.createVariable("bke_avg", "f4", ("lon", "lat",))
        bke_avg[:] = self.avg_bke.transpose()
        bke_std = rot.createVariable("bke_std", "f4", ("lon", "lat",))
        bke_std[:] = self.std_bke.transpose()
        
        rot.close()


    def bleach(self,ckp, sfloor, g, dh, zob):
        """The method bleach  is used to apply a simple BL
        model. References in Annex 1 Compiling
        oceanographic layers to the EUSeaMap 2019.
        Technical Report and Maraldi et al., 2013, """
        ndz = dh / (2 * zob)
        cd = 1. / (2.5 * np.log(ndz)) ** 2.
        ustar = np.sqrt(cd * g ** 2.)
        return ustar
