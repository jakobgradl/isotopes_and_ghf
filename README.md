# Isotopes and GHF
Stable isotopes as a proxy for assessing geothermal heat flow in Greenland
Gradl et al., manuscript submitted to Frontiers in Earth Science

--------------

Code files for the subglacial hydrological routing and isotope transport

The implemented routing algorithm is the MD-inf algorithm of Seibert & McGlynn, 2007, https://doi.org/10.1029/2006WR005128.

We adapt their algorithm to the irregular triangular ISSM-mesh used in our study. 

--------------

This is the order in which to run the different functions:

av, rv, vw = flowpath(md, phi) # IrrDMinf.py, need to provide an estimate of subglacial hydropotential phi

bua, bvol, bflux, wshed, rv_new, vw_new, melt_new = basalrouting(md, av, rv, vw, va) # BasalRouting.py

age = nyeplusmelt(md, melt_new) # NyePlusMelt.py

isocomp = AgeToIso(md,age) # AgeToIso.py

basaliso, vertexiso = basalisotopes(md, av, rv_new, vw_new, va, melt_new, bvol, isocomp) # BasalIsotopeRouitng.py

--------------

md: loaded ISSM model providing the mesh geometry ery

av: adjacent vertices (mesh vertices adjacent to each vertex)

rv: receiving vertices (vertices receiving flux from a given vertex)

vw: vertex weights (fraction of flux received by every rv from a given vertex)

bua: basal upstream area (area upstream of a given vertex)

bvol: basal volume (accumulated water volume at each vertex, not normalized!)

bflux: basal flux (accumulated volume normalized by the associated vertex area)

rv_new, vw_new, melt_new: correct spurious sources/sinks arising from the irregular mesh geometry

age: age of basal ice

isocomp: isotope composition of the basal ice

basaliso: isotope composition of the basal meltwater