# Seabottom Statistical Libraries v4.1
This software has been developed in the course of several different phases of "The
European Marine Observation and Data Network (EMODnet) broad-scale seabed habi-
tat map for Europe", known as EUSeaMap since the start of phase 2 in 2013. (see https:
//www.emodnet-seabedhabitats.eu). It started as a series of dfferent FORmula TRANs-
lator (FORTRAN) codes aimed to provide an estimate of the kinetic energy close to the
sea bottom due to the effects of waves and currents in the Mediterranean Sea.
The process of evaluating the value of a parameter on a curvilinear
surface close to the bottom from a more or less regular three-dimensional grid is not a
simple problem of interpolation, since near the sea bottom important physical processes
occur that make the boundary layer a place where the vertical variations of `dynamical'
variables are generally non-linear. Also the depth of the boundary layer is not constant,
so that the only way to infer a reasonable value for the velocity or the temperature near
the bottom (here the distance of 1 m from the bottom surface is taken as a reference
value for energy) is to try to simulate the physical processes and build a model of the
boundary layer. (or the boundary layers, since waves and currents have their own way
of dealing with the bottom). Said that, there are still parameters which are
'less dynamic' and have a sufficiently smooth variation in space that the interpolation is
sometimes still a viable option (see the Bchem library section and the discussion on the
validity of results for biochemical parameters). In time, codes were modified to cover
also Black Sea, Iberia-Biscay-Ireland, and Macaronesia. The code for the evaluation of
wave kinetic energy depends on time series of significant wave and peak period fields at
the surface as obtained from high resolution statistical wave models of fourth generation
like Wave Model (WAM) or Wave Watch III (WWIII) and a high resolution
bathymetry in the area of interest. The estimate is obtained by a simplified method,
and described in the Annex "Compiling oceanographic layers". 
In order to obtain an estimate of the kinetic energy due to current at a small distance
from the bottom, a simple adaptation of the boundary layer algorithm used in the oceanographic 
model used at Copernicus Marine Service (CMEMS) for the Mediterranean and
Black Sea has been devised, with the idea that the application of the same scheme as a
postprocessing of the 3-dimensional current field would provide a more consistent evaluation 
near the bottom ([8, 1]). Also in this case, the method is highly dependent on the
available bathymetry of the local area considered. All methods are clearly gross 
simplifications of the real physical processes, acceptable only in providing long term statistics.
Since 2016 CMEMS has provided sufficiently high resolution simulations of physical
characteristics of the European seas to be successfully used to evaluate the bottom ki-
netic energy. In 2022 it was decided to implement all software in the form of a library
which would provide a way to easily extract the desired estimates using only CMEMS
oceanographic outputs and EMODnet bathymetryes.
The code is written in Python for several (good) reasons. The first is that it is a widely 
used language which has been adopted by almost all scientific communities. Second is that 
is particularly apt to work in object-oriented environments, has effective methods to deal 
with Network Common Data Form (NetCDF) and leads to surprisingly simple installation 
procedures for the libraries. Actually only the path of the library needs to be set in 
the script/notebook. 
On the other hand, even though some efforts have been spent to improve vectorialization, 
the execution of a long term statistic evaluation could take some time. Not too
long though (something between a coffee and a compact lunch), so in order to keeps the
installation of the library and the code as simple as possible, no parallelisation has been
considered at the present stage. Another practical need is the evaluation of biochemical
parameters (like sea water salinity (salinity), mole concentration of dissolved molecular
oxygen in sea water (O2) (oxygen) or mole concentration of chlorophyll(a) in sea water
(CHL) (chlorophyll)) near the bottom. In this case no effort has been made to provide a
boundary layer model, the user must be aware that the evaluation is based on the simple
extraction of the parameter data at the level which is closest to the bottom. In the case
of the temperature, CMEMS provides for the temperature at the bottom, so the library
provides only the framework for the evaluation of the desired statistics. In order to keep
the management of data sufficiently simple, the adopted bathymetry is the 2016 version,
based on 12 data blocks. The library for interpolating the European Marine Observation
and Data Network (EMODnet) 2020 product works on many blocks of high resolution
data. Users may decide to use this version in order to gain some accuracy at the price
of a longer pre-processing. The libraries have been devised to work on several CMEMS
oceanographic products at regional scale (i.e. all products based on a latitude/longitude
WGS84 reference coordinate system - The Artic Ocean product is an example of product
which cannot be processed as being Polar Stereographic. See Test section) the simplicity
of usage in terms of auto-consistency, means that the user is not required to specify the
area of application or the time periods, the library simply works on all the files and reads
the parameters from the Network Common Data Form (NetCDF) input files. This of
course requires some extra care in the preparation of input data obtained from CMEMS.
