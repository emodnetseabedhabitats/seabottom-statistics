#    Copyright 2022 Roberto Inghilesi, Alessandro Mercatini

#    This file is part of seabottom 4.1

#    ibathy2020.py is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    icke.py is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.

#    You should have received a copy of the GNU Lesser General Public License
#    along with ibathy2020.py.  If not, see <http://www.gnu.org/licenses/>.

from netCDF4 import Dataset
import numpy as np

class Node:
    """The Class Node is written to provide an object with attributes latitude and longitude of a grid node. It is used in the Classes Box, Gridinfo and Sbathy."""
    @property
    def plon(self):
        return self._plon

    @property
    def plat(self):
        return self._plat

    def __init__(self, plon, plat):
        self._plon = plon
        self._plat = plat

    def __str__(self):
        return "([lon {} - lat {}])".format(self.plon, self.plat)


class Box:
    """The Class Box is written to provide an object with attributes two grid nodes p0,p1 defining a grid box.
    It is used in the Classes Gridinfo and Sbathy."""
    @property
    def p0(self):
        return self._p0

    @property
    def p1(self):
        return self._p1

    def __init__(self, p0, p1):
        self._p0 = p0
        self._p1 = p1

    def __str__(self):
        return "([p0{} \n p1{}])".format(self.p0, self.p1)


class Gridinfo(object):
    """The Class Gridinfo is written to provide an object with attributes defining the properties of the fields data read either from a model output or from the Emodnet bathymetry files. It is used in the Classes Sbathy, Sbenergy, Stands, Wbenergy."""
    @property
    def nlons(self):
        return self._nlons

    @property
    def ntimes(self):
        return self._ntimes

    @property
    def nlats(self):
        return self._nlats

    @property
    def nlevs(self):
        return self._nlevs

    @property
    def lons(self):
        return self._lons

    @property
    def lats(self):
        return self._lats

    @property
    def nresx(self):
        return self._nresx
    
    @property
    def nresy(self):
        return self._nresy
        
    @property
    def z(self):
        return self._z

    @property
    def box(self):
        return self._box

    @property
    def missval(self):
        return self._missval
    
        
    @nlevs.setter
    def nlevs(self, value):
        self._nlevs = value
    
    @z.setter
    def z(self, value):
        self.__z[i][j] = value

    def __init__(self, filename, ig):
        lon_lst = ['lon', 'longitude', 'rlon']
        lat_lst = ['lat', 'latitude', 'rlat']
        mv_lst = ['missing_value', '_FillValue']
        if ig == 0:
            fh = Dataset(filename, mode='r', format="NETCDF4")
            ntimes = fh.dimensions['time'].size
            # dimensions
            for l in lon_lst:
                if l in fh.dimensions.keys():
                    nlons = fh.dimensions[l].size
                    break
            for l in lat_lst:
                if l in fh.dimensions.keys():
                    nlats = fh.dimensions[l].size
                    break
            # variables
            for l in lon_lst:
                if l in fh.variables.keys():
                    lons = np.array(fh[l])
                    break
            for l in lat_lst:
                if l in fh.variables.keys():
                    lats = np.array(fh[l])
                    break
            self._nlons = int(nlons)
            self._nlats = int(nlats)
            #self._nlevs = int(dlev)
            self._ntimes = int(ntimes)
            z = np.zeros([nlats, nlons], dtype=float)
            self._missval = 0
        elif ig > 0:
            gh = Dataset(filename, mode='r', format="NETCDF4")
            nlons = gh.dimensions['lon'].size
            nlats = gh.dimensions['lat'].size
            lons = np.array(gh['lon'])
            lats = np.array(gh['lat'])
            dlev = 0
            self._nlons = int(nlons)
            self._nlats = int(nlats)
            self._nlevs = int(dlev)
            mv=gh['elevation']._FillValue
            dept = gh['elevation']
            dep=np.where(np.isnan(dept), 0., dept)
            z = np.array(dep)  # *np.float64(sf)+np.float64(off)
            z = np.where(z > -2., 0., z) #check limiter
            self._missval = 0
        self._lons = lons
        self._lats = lats
        self._z = z
        ndel1=nlons/(lons[-1]-lons[0])
        ndel2=nlats/(lats[-1]-lats[0])
        self._nresx=int(ndel1+.1)
        self._nresy=int(ndel2+.1)
        p0 = Node(lons[0], lats[0])
        p1 = Node(lons[nlons - 1], lats[nlats - 1])
        bx = Box(p0, p1)
        self._box = bx


class Sbathy(object):
    """The Class Sbathy is written to provide an evaluation of bathymetry interpolated on the
    oceanographic model grid based on EMODNET bathymetry (2020). It is meant as an aid to the
    evaluation of seabottom kinetic energy or temperature and salinity obtained from CMEMS
    models with fixed z levels.usage- 1st step evaluate the bathymetry by using an instance of
    Sbathy taking inpufile of the oceanographic model describing the grid and a list of names
    of EMODNETbathymetry files. Not all files belonging to the EMODNET distribution need
    necessarily to be referenced but the referenced files must completely cover the area. See
    EMODNET documentation or reference all files, the program will consider only those that
    are necessary. g=Sbathy(filename_ts,0,*nb or g=Sbathy.fromfile('namefile.nc') if a
    previous instance has produced the interpolated bathymetry using for example
    g.to_file('ibi_bat-XX.nc') g = Sbathy(filename, 0, *nb, printcheck=1) prints the boundary
    index check for each eumodnet bathymetry block. Technical reference can be found in the Anne
    "Compiling oceanographic layers" to the EUSeaMap 2019 - Technical Report"""
    
    @property
    def gmed(self):
        return self._gmed
    @classmethod
    def from_file(cls, nfile):
        """The Class method from_file is useed to build an istance Sbathy
        using an available file of bathymetry instead of re-interpolating from Emodnet data."""
        fh = Dataset(out + '/' + nfile, mode='r', format="NETCDF4")
        nlats = fh.dimensions['lat'].size
        nlons = fh.dimensions['lon'].size
        lons = np.array(fh["lon"])
        lats = np.array(fh["lat"])
        instance = cls(filename, 0)
        bat = np.array(fh["bat"]).transpose()
        for ip in range(nlats):
            for jp in range(nlons):
                instance.gmed.z[ip, jp] = bat[ip][jp]
        return instance

    def __init__(self, filename, inc, *args,checkb=0):
        gmed = Gridinfo(filename, 0)
        self._gmed = gmed
        for g in args:
            print(g)
            gemo = Gridinfo(g, 1)
            pn = self.get_ibox(gmed, gemo)
            ibox = Box([pn[0], pn[2]], [pn[1], pn[3]])
            if pn[0] == -1:
                continue
            for ip in range(pn[2], pn[3]):
                for jp in range(pn[0], pn[1]):
                    p = Node(gmed.lons[jp], gmed.lats[ip])
                    gmed.z[ip, jp] = self.h480(p, gemo)

    def to_file(self, fl):
        """The method to_file is used to output the ncdf4 file of bathymetry interpolated from Emodnet data on the model grid."""
        rot = Dataset(fl, "w", format="NETCDF4") #renoved out path in rot A.
        lon = rot.createDimension("lon", self.gmed.nlons)
        lat = rot.createDimension("lat", self.gmed.nlats)
        latitudes = rot.createVariable("lat", "f4", ("lat",))
        longitudes = rot.createVariable("lon", "f4", ("lon",))
        latitudes[:] = self.gmed.lats
        longitudes[:] = self.gmed.lons
        bat = rot.createVariable("bat", "f4", ("lon", "lat",))
        bat[:] = self.gmed.z.transpose()
        rot.close()

    def h480(self, p, gemo):
        """The method h480 is used to interpolate from Emodnet data on the model grid.
        it averages those available of the 9 points of the Emodnet file closer to the
        gridpoint of the model. Reference system is in the gemo (Emodnet) framework.The method
        is used in the constructor of Sbathy"""
        # reference system is in the gemo framework
        j = int(gemo.nresx * (p.plon - gemo.box.p0.plon))
        i = int(gemo.nresy * (p.plat - gemo.box.p0.plat))
        h480 = 0.
        dumb = 0.
        n = 0
        if i == 0 or j == 0:
            if (gemo.z[i, j] != gemo.missval):
                n = 1
                dumb = gemo.z[i, j]
        elif j == gemo.nlons:
            if (gemo.z[i, j] != gemo.missval):
                n = 1
                dumb = gemo.z[i, j]
        elif i == gemo.nlats:
            if (gemo.z[i, j] != gemo.missval):
                n = 1
                dumb = gemo.z[i, j]
        else:
            for ic in range(-1, 1, 1):
                for jc in range(-1, 1, 1):
                    ic = int(ic)
                    jc = int(jc)
                    if (gemo.z[i + ic, j + ic] != gemo.missval):
                        dumb += gemo.z[i + ic, j + jc]
                        n += 1
        if n == 0:
            dumb = 0.
            return 0.
        else:
            return dumb / n
    def get_ibox(self, gmed, gemo):
        """The method get_ibox is used to interpolate from Emodnet data on the model grid.it defines a box in the Emodnet framework containing the points of the model grid to be interpolated. If no point is present the flag -1 is returned. The method is used in the constructor of Sbathy"""
        fmed = gmed.box.p1.plon  # end of med bathy long.
        xmed = gmed.box.p0.plon  # start of med bathy long.
        xmod = gemo.box.p0.plon  # start of emodnet block long
        fmod = gemo.box.p1.plon  # end of emodnet block long
        res1 = 0
        res2 = 0
        if fmed < xmod or xmed > fmod:  # no intersection #check:
            return -1, -1, -1, -1
        else:
            if xmed < fmod:
                if xmed <= xmod:
                    temp = (xmod - xmed) * gmed.nresx
                    i = int(temp+1)
                    j = (fmod - xmed) * gmed.nresx + 1
                    res1 = max(i, 0)
                    res2 = min(int(j), gmed.nlons)
                elif xmed > xmod:
                    i = 0
                    j = (fmod - xmed) * gmed.nresx +1
                    res1 = i
                    res2 = min(int(j), gmed.nlons)
        fmed = gmed.box.p1.plat  # send of med bathy lat.
        xmed = gmed.box.p0.plat  # start of med bathy lat
        xmod = gemo.box.p0.plat  # start of emodnet block lat
        fmod = gemo.box.p1.plat  # end of emodnet block lat
        #print(fmed, xmed, xmod, fmod, '#####################')
        res3 = 0
        res4 = 0

        if fmed < xmod or xmed > fmod:  # no intersection
            return -1, -1, -1, -1
        else:
            if xmed < fmod:
                if xmed >= xmod:
                    i = 0
                    j = (fmod - xmed) * gmed.nresy +1
                elif xmed < xmod:
                    temp = (xmod - xmed) * gmed.nresy
                    i = int(temp)
                    j = (fmod - xmed) * gmed.nresy +1
                res3 = max(i, 0)
                res4 = min(int(j), gmed.nlats)
        return res1, res2, res3, res4

