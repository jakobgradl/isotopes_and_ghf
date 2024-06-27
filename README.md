# GHF-from-isotopes
Refining Greenland Geothermal heat flow through stable isotope analysis

This is the order in which to run the different funcitons:

av, rv, vw, dw = flowpath(md, phi)
bua, bvol, bflux, wshed, rv_new, vw_new, melt_new = basalrouting(md, av, rv, vw, va)
age = nyeplusmelt(md, melt_new)
isocomp = AgeToIso(md,age)
basaliso, vertexiso = basalisotopes(md, av, rv_new, vw_new, va, melt_new, bvol, isocomp)
