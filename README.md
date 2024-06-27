# GHF-from-isotopes
Refining Greenland Geothermal heat flow through stable isotope analysis

This is the order in which to run the different funcitons:

av, rv, vw = flowpath(md, phi) # IrrDMinf.py, you will need to provide an estimate of subglacial hydropotential

bua, bvol, bflux, wshed, rv_new, vw_new, melt_new = basalrouting(md, av, rv, vw, va) # BasalRouting.py

age = nyeplusmelt(md, melt_new) # NyePlusMelt.py

isocomp = AgeToIso(md,age) # AgeToIso.py

basaliso, vertexiso = basalisotopes(md, av, rv_new, vw_new, va, melt_new, bvol, isocomp) # BasalIsotopeRouitng.py
