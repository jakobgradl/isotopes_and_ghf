
# BasalIsotopeRouting.py

# Compute the isotopic composition of the basal meltwater
# d18O values are treated as a conservative tracer to compute the transport and mixing of isotope signals

# The structure of the routing algorithm is the same as BasalRouting.py
# You will need the output of BasalRouting.py for this computation

#----------------




# python packages

import numpy as np
import csv
import math as m
from importlib import reload
from copy import deepcopy
from netCDF4 import Dataset

# ISSM packages
from model import *
from InterpFromGridToMesh import InterpFromGridToMesh







        # Main function to compute the isotope routing
        # The isotope routing is set up in the same way as the basalrouting() and dflux() in BasalRouting.py

def basalisotopes(LoadedModel, AdjacentVertices, ReceivingVertices, VertexWeights, VertexArea, Melt, BasalRouting, IsoComp, Test='noTest'): # {{{
    # LoadedModel is your ISSM model
    # AdjacentVertices, ReceivingVertices, VertexWeights is the output of IrrMDinf.py
    # VertexArea is the output of VertexArea.py


        # Variables
    
    md = deepcopy(LoadedModel)
    av = AdjacentVertices # List of adjacent vertices at each vertex
    rv = ReceivingVertices # List of all adjacent vertices at vertex that receive water from vertex
    vw = VertexWeights # Corresponding weights of how to partition the routing between all receiving vertices
    melt = deepcopy(Melt) # Basal melt rate
    # ReceivingVertices, VertexWeights, and Melt are part of the output of BasalRouting.py, they have received treatment for edge-cases
    va = VertexArea
    totalflux = deepcopy(BasalRouting[:md.mesh.numberofvertices2d]) # that's bvol in the basalrouting()
    isotopes = deepcopy(IsoComp[:md.mesh.numberofvertices2d]) # This is the d18O composition of the basal ice, the output of AgeToIso.py


        # Load data:

    # Fill holes in the basal ice d18O data by linear interpolation
    for j in np.arange(md.mesh.numberofvertices2d):
        if m.isnan(isotopes[j]):
            isotopes[j] = np.mean(isotopes[av[j][1:]])

            
        # data arrays:
        
    meltmask = np.where(melt > 0, 1, 0) # mask of the basal meltarea

    # arrays for the computation
    vertexiso = np.nan * np.ones(md.mesh.numberofvertices2d) # total isotope amount at each vertex (needed for computation, no physical meaning)
    isoflag = np.ones(md.mesh.numberofvertices2d,dtype=int) # a marker used in the recursion
    meltflux = melt[:md.mesh.numberofvertices2d] * va # areal basal meltwater flux associated with each mesh vertex

    # array for the output
    basalisotopes =  np.zeros(md.mesh.numberofvertices) # d18O values of the water in the basal drainage system


        # Main computation:

    for vertex in np.arange(md.mesh.numberofvertices2d):
        if isoflag[vertex] and meltmask[vertex] and not m.isnan(isotopes[vertex]):
            vertexiso[vertex], isoflag[vertex], meltflux[vertex] = disotope(md, vertex, vertexiso, isoflag, meltflux, va, av, rv, vw, isotopes, meltmask)

    basalisotopes[:md.mesh.numberofvertices2d] = vertexiso / totalflux # composite d18O value of basal water per unit volume of water


    print('    Calculation Basal Composite Isotope Value: Done')
    return basalisotopes, vertexiso

# }}}






def disotope(LoadedModel, Vertex, VertexIso, IsoFlag, MeltFlux, VertexArea, AdjacentVertices, ReceivingVertices, VertexWeight, IsoComp,  MeltMask): # {{{
    # adjacent to dflux() in BasalRouting.py (Tarboton, 1997, https://doi.org/10.1029/96WR03137)
    # isntead of upstream area, here we accumulate the d18O values of all basal meltwater produced upstream
    
    md = LoadedModel
    vertex = Vertex
    vertexiso = VertexIso
    isoflag = IsoFlag
    va = VertexArea
    meltflux = MeltFlux
    av = AdjacentVertices
    rv = ReceivingVertices
    vw = VertexWeight
    isocomp = IsoComp
    meltmask = MeltMask

    if m.isnan(isocomp[vertex]):
        count = 0
        isocomp[vertex] = 0
        for j in av[vertex][1:]:
            if m.isnan(isocomp[j]):
                pass
            else:
                isocomp[vertex] = isocomp[vertex] + isocomp[j]
                count = count + 1
        isocomp[vertex] = isocomp[vertex] / count
        
    if m.isnan(meltflux[vertex]):
        meltflux[vertex] = 0
        
    vertexiso[vertex] = isocomp[vertex] * meltflux[vertex]
    isoflag[vertex] = 0
    
    for i in av[vertex][1:]:
        if vertex in rv[i]:
            p = vw[i][np.where(rv[i]==vertex)]
            if isoflag[i] and meltmask[i] and not m.isnan(isocomp[i]):
                vertexiso[i], isoflag[i], meltflux[i] = disotope(md, i, vertexiso, isoflag, meltflux, va, av, rv, vw, isocomp, meltmask)

            if not m.isnan(vertexiso[i]):
                try:
                    vertexiso[vertex] = vertexiso[vertex] + (p * vertexiso[i])
                except:
                    print('error with', vertex, i, p)

    return vertexiso[vertex], isoflag[vertex], meltflux[vertex]

# }}}


