
# BasalRouting.py
# Calculate the basal routing and water flux

# You need to have run VertexArea.py and IrrMDinf.py beforehand as you need their outputs

# --------------------



# python packages
import numpy as np
import math as m
from importlib import reload
from copy import deepcopy

# ISSM packages
from model import *
from ll2xy import ll2xy # this one's not strictly necessary





        # Main function to compute the water routing

def basalrouting(LoadedModel, AdjacentVertices, ReceivingVertices, VertexWeights, VertexArea): #{{{
    # LoadedModel is your ISSM model
    # AdjacentVertices, ReceivingVertices, VertexWeights is the output of IrrMDinf.py
    # VertexArea is the output of VertexArea.py


        # Variables
    
    md = deepcopy(LoadedModel)
    av = AdjacentVertices # List of adjacent vertices at each vertex
    rv = deepcopy(ReceivingVertices) # List of all adjacent vertices at vertex that receive water from vertex
    vw = deepcopy(VertexWeights) # Corresponding weights of how to partition the routing between all receiving vertices
    va = VertexArea

    melt = deepcopy(md.basalforcings.groundedice_melting_rate) # basal meltrate
    # melt = deepcopy(md.results.SteadystateSolution.BasalforcingsGroundediceMeltingRate)

    
        # data arrays

    # arrays for the computation
    upstreamarea = np.zeros(md.mesh.numberofvertices2d) # accumulated upstream area at each vertex (the primary output of the routing algorithm)
    vertexflux = np.zeros(md.mesh.numberofvertices2d) # water flux at each vertex
    fluxflag = np.ones(md.mesh.numberofvertices2d,dtype=int) # a little marker used in the computation
    flowstop = np.zeros(md.mesh.numberofvertices2d) # counter for watershed

    # arrays for the output
    basalvolume =  np.zeros(md.mesh.numberofvertices) # routed water volume
    basalupstreamarea = np.zeros(md.mesh.numberofvertices) # routed basal area
    basalflux = np.zeros(md.mesh.numberofvertices) # vertex-area normalised basal water flux
    watershed = np.zeros(md.mesh.numberofvertices) # highlights vertices that don't receive any water

    posmelt = np.where(melt > 0, 1, 0) # mask of positive non-zero basal melt
    
    
        # Handling of Edge-Cases

    # 1)
    # if a vertex shows a local high in hydropotential, i.e., doesn't get any are because it gives to al av's
    # smooth that out by interpolation
    # this is already done in flowpath in IrrMDinf.py, just as a reminder here :)
#    for j in np.arange(md.mesh.numberofvertices2d):
#        for i in av[j][1:]:
#            if phi[j] > phi[i]:
#                localhigh = 1
#            else:
#                localhigh = 0
#                break
#        if localhigh:
#            phi[j] = np.mean(phi[av[j][1:]])


    # 2)
    # if melt is negative on a single vertex but positive all around (probably an artifact of the simulation)
    # fill in the gap by interpolation
    count = 0
    for j in np.arange(md.mesh.numberofvertices2d):
        if not posmelt[j]:
            for i in av[j][1:]:
                if posmelt[i]:
                    ismelt = 1
                else:
                    ismelt = 0
                    break
            if ismelt:
                melt[j] = np.mean(melt[av[j][1:]])
                posmelt[j] = 1
                count += 1
    print("    Negative Melt filled in on ", count, "nodes")


    # 3)
    # currently, there are vertices with 1) individual vw's greater than 1 and 2) sum(vw) greater than 1
    # 1) these vertices give more than 100% of their area to one av and an equal negative amount to another av
    # ...so they balance out to 100% overall
    # 2) these vertices give more than 100% of their area to their av's
    # 3) there are no vertices that give less than 100% of their area to their av's, so that's good (except for 0%)
    # relay those to the following edge-case by removing their rv's and vw's
    count = 0
    for j in np.arange(md.mesh.numberofvertices2d):
        if np.sum(vw[j]) > 1.01 or np.any(vw[j] > 1):
#            rv[j] = np.empty(0)
#            vw[j] = np.empty(0)
            count += 1
    print("    ", count, "nodes that route MORE than 100 % of area")
    # that is certainly not the proper way to do it, it's a bug in the flowpath() that needs fixing
    # but on account of being short on time I'll do it like this for now
    # that way, at least the area is conserved during routing...

    count = 0
    for j in np.arange(md.mesh.numberofvertices2d):
        if (np.sum(vw[j]) < 0.98 and not rv[j].size == 0) or np.any(vw[j] < 0):
            count += 1
    print("    ", count, "nodes that route LESS than 100 % of area")



    # 4)
    # if a vertex doesn't route anything to another vertex (not really sure how that can happen...)
    # give area to all adjacent vertices in equal parts that don't themselves route to this vertex
    count = 0
    for j in np.arange(md.mesh.numberofvertices2d):
        if rv[j].size == 0:# and not md.mesh.vertexonboundary[j]:
            
            giveareathere = np.empty(0)
            count += 1
            for i in av[j][1:]:
                if j in rv[i]:
                    pass
                else:
                    giveareathere = np.append(giveareathere, i)

            try:
                rv[j] = giveareathere
                vw[j] = (1 / giveareathere.size) * np.ones(giveareathere.size)
            except:
                pass
#                print(j, '\t No receiving vertices', '\t', rv[j])
    print("    ", count, "nodes that don't route anything to other nodes")

            
    # Note: 4) needs to be before 5) so that vw[p].size is never 0


    # 5)
    # if a vertex doesn't get any area from another vertex (that's an artifact of the triangular grid, refer to supplements of Gradl et al., 2024)
    # get 0.1 of area of each adjacent vertex, that aren't themselves receiving vertices (this is totally arbitrary)
    count = 0
    for j in np.arange(md.mesh.numberofvertices2d):
        for i in av[j][1:]:
            if j in rv[i]:
                getsarea = 1
                break
            else:
                getsarea = 0
        if not getsarea:
            count += 1
            getareahere = np.isin(av[j][1:], rv[j], invert=True) # i.e., is adjacent vertex not in receiving vertex
            if np.squeeze(np.where(getareahere)).size > 1:
                for p in av[j][np.squeeze(np.where(getareahere))+1]:
                    #print(j, p, getareahere)
                    rv[p] = np.append(rv[p], j)
                    vw[p] = vw[p] - (0.1 / vw[p].size)
                    vw[p] = np.append(vw[p], 0.1)
            else:
                pass # still have to figure out how that's possible
    print("    ", count, "nodes that don't get any area")





    # 6)
    # Test whether any nodes route water to each other (which is not possible)
    count = 0
    for vertex in np.arange(md.mesh.numberofvertices2d):
        for i in av[vertex][1:]:
            if i in rv[vertex] and vertex in rv[i]:
                count += 1
    print("    ", int(count/2), "pairs of nodes that route area to each other")

    
        # End Edge Cases




        # If you want, you can only compute a single drainage basin inside the domain
        # Drainage basin is computed recursively, so you define the end-point/outlet
        # Here are two examples:

    # Area draining to EastGRIP
    # Position of EastGRIP as 5x5 km cell, from Lat: 75°38'N, Lon: 35°60'W
    [egripx, egripy] = ll2xy(75.5, -35.5, 1)
    egrip = np.where(np.logical_and(md.mesh.vertexonbase[:md.mesh.numberofvertices2d], np.logical_and(np.logical_and(md.mesh.x2d>(egripx-5000),md.mesh.x2d<(egripx+5000)), np.logical_and(md.mesh.y2d>(egripy-5000),md.mesh.y2d<(egripy+5000)))))
    egrip = np.squeeze(egrip)

    # Area draining to NEGIS outlet
    # NEGIS outlet 1x4 ° area at 79NGlac and ZachIsstr, from Lat: 78.5-79.5 °N, Lon: 20-24 °W (that is square in areal image)
    [negisxhigh, negisyhigh] = ll2xy(79.5, -20.0, 1)
    [negisxlow, negisylow] = ll2xy(78.5, -30.0, 1) #-24.0, 1)
    negis = np.where(np.logical_and(md.mesh.vertexonbase[:md.mesh.numberofvertices2d], np.logical_and(np.logical_and(md.mesh.x2d>negisxlow,md.mesh.x2d<negisxhigh), np.logical_and(md.mesh.y2d>negisylow,md.mesh.y2d<negisyhigh))))
    negis = np.squeeze(negis)




    

        # Main computation:

        
    # This is for routing through No-Melt areas:
    melt = np.where(posmelt, melt, 0)
    # To enable routing through No-Melt areas, comment-out the posmelt condition below and in dflux()
    # this doesn't work for the isotope routing, though


    # un-comment desired drainage basin
    
#    for vertex in negis:    # NEGIS drainage basin
#    for vertex in egrip:    # EGRIP drainage basin
    for vertex in np.arange(md.mesh.numberofvertices2d):    # drainage on full domain
        # cycle through all vertices to compute their upstream area with dflux()
        
        if fluxflag[vertex] and posmelt[vertex]:
            # fluxflag marks if a vertex has already been computed to avoid double computation
            # posmelt ensures no routing through no-melt areas
            
            recursioncount = 0
            upstreamarea[vertex], vertexflux[vertex], fluxflag[vertex], flowstop[vertex] = dflux(md, vertex, upstreamarea, vertexflux, fluxflag, flowstop, va, av, rv, vw, recursioncount, posmelt, melt)




        # put 2d arrays into 3d arrays for the output
    
    basalupstreamarea[:md.mesh.numberofvertices2d] = upstreamarea
    basalvolume[:md.mesh.numberofvertices2d] = vertexflux
    basalflux[:md.mesh.numberofvertices2d] = basalvolume[:md.mesh.numberofvertices2d] / va
    watershed[:md.mesh.numberofvertices2d] = flowstop

    print('    Calculation Basal Water Routing: Done')
    return basalupstreamarea, basalvolume, basalflux, watershed, rv, vw, melt

# }}}












def dflux(LoadedModel, Vertex, UpstreamArea, VertexFlux, FluxFlag, FlowStop, VertexArea, AdjacentVertices, ReceivingVertices, VertexWeights, RecursionCount, PosMelt, Melt): # {{{
    # This is the actual algorithm for recursively computing the accumulated upstream area at a given vertex
    # The structure of the algorithm is described in Tarboton (1997), https://doi.org/10.1029/96WR03137
    
    md = LoadedModel
    currentvertex = deepcopy(Vertex) # either the vertex passed on by basalrouting() or the vertex of the preceeding recursion step passed on by dflux()
    upstreamarea = UpstreamArea
    vertexflux = VertexFlux
    fluxflag = FluxFlag
    flowstop = FlowStop
    vertexarea = VertexArea
    av = AdjacentVertices
    rv = ReceivingVertices
    vw = VertexWeights
    rc = RecursionCount
    posmelt = PosMelt
    melt = Melt
    # !!! It's important to NOT deepcopy vertexflux and fluxflag, so that the recursion accesses the same objects in each step

    upstreamarea[currentvertex] = vertexarea[currentvertex] # set the current vertex's area as its upstream area
    vertexflux[currentvertex] = melt[currentvertex] * vertexarea[currentvertex] # same with the water flux

    if m.isnan(upstreamarea[currentvertex]):
        upstreamarea[currentvertex] = 0
    if m.isnan(vertexflux[currentvertex]):
        vertexflux[currentvertex] = 0
        
    fluxflag[currentvertex] = 0 # mark that we have already done this current vertex
    adjacentflow = 0

    # cycle through all adjacent vertices, check if they route water to current vertex, and if the do call dflux() on them to compute their area
    # this is the recursive part
    # add all upstream area to the upstreamarea of current vertex
    for i in av[currentvertex][1:]:
        if currentvertex in rv[i]:
            adjacentflow = adjacentflow + 1
            p = vw[i][np.where(rv[i]==currentvertex)]
            if fluxflag[i] and posmelt[i]:
                rc = rc + 1
                # recursion:
                upstreamarea[i], vertexflux[i], fluxflag[i], flowstop[i] = dflux(md, i, upstreamarea, vertexflux, fluxflag, flowstop, vertexarea, av, rv, vw, rc, test, posmelt, melt)

            try:
                upstreamarea[currentvertex] = upstreamarea[currentvertex] + (p * upstreamarea[i])
                if not m.isnan(vertexflux[i]):
                    vertexflux[currentvertex] = vertexflux[currentvertex] + (p * vertexflux[i])
            except:
                print('error with', currentvertex, i, p)

    if adjacentflow == 0:
        flowstop[currentvertex] = 1 # i.e., this vertex doesn't get any area from another vertex
        
    return upstreamarea[currentvertex], vertexflux[currentvertex], fluxflag[currentvertex], flowstop[currentvertex]

# }}}
