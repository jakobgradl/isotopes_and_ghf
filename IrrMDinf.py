

# IrrMDinf.py

# Hydrological routing algorithm MDinf (Seibert and McGlynn, 2007, doi:10.1029/2006WR005128)
# adapted for irregular triangular meshes.

# This algorithm determines which vertices receive water from any given vertex and the partitioning of flow for multi-directional routing
# it produces the input for BasalRouting.py

# Specifically, the algorithm calculates the upstream area at each vertex.
# Meaning, the algorithm doesn't route waterflow but mesharea

# You will need to provide a map of the subglacial hydropotential

# Original use-case was built on ISSM with python interface
# So the data structure of the model is the one used in ISSM
# That you may need to adjust to your specific case

# Most variables that are constructed here have the 2d-structure var=(x,y)
# x being the mesh vertices, y being a list of data associated with each vertex

# --------------------



# Python packages
import numpy as np
import math as m
from importlib import reload
from copy import deepcopy

# ISSM packages
from model import *
from ll2xy import ll2xy






        # Main function

def flowpath(LoadedModel, HydroPot):
    # LoadedModel is your ISSM model, HydroPot is the basal hydropotential map
    # refer to the subroutines for explanations of variables
    
    md = deepcopy(LoadedModel)
    phi = deepcopy(HydroPot)

    vertcon, e2d, adjvert = adjacentvertices(md)
    nphi = new_phi(md, phi, adjvert)
    norm = normals(md, nphi, e2d)
    elor = elementorientation(md, norm)
    vertvect = vertexvectors(md, adjvert)
    eflowdir, dirflags = elementflowdirections(md, vertcon, e2d, adjvert, elor, vertvect, nphi)
    sdfdirs, flowtoedge = steepestslopes(md, eflowdir, dirflags, adjvert, vertvect, nphi)
    dirslope, vertID = directionslopes(md, sdfdirs, adjvert, vertvect, nphi, vertcon, eflowdir, flowtoedge, norm)
    dirweight = directionweight(md, sdfdirs, dirslope)
    receivvert, vertweights = receivingvertices(md, dirweight, vertID, flowtoedge, e2d,  vertvect, adjvert, sdfdirs, nphi)

    print('    Calculation Flowpath: Done')

    return adjvert, receivvert, vertweights
    # adjvert > vertices adjacent to vertex
    # receivvert > subset of adjacent vertices that receive flow from vertex
    # vertweights > partitioning of flow from vertex between the receiving vertices








        # Sub-routines of the main function

def adjacentvertices(LoadedModel):
    # determine all vertices that are adjacent to vertex
    
    md = deepcopy(LoadedModel)

    basalvertices = np.squeeze(np.where(md.mesh.vertexonbase)) # collapse 3d model to 2d, we are only interested in the basal layer
    vertexconnectivity = NodeConnectivity(md.mesh.elements2d, md.mesh.numberofvertices2d)
    # vertexconnectivity is a two-dimensional array.
    # First dimension is the total number of mesh vertices.
    # Second dimention is an array of length 100. The first positions list the mesh elements adjacent to vertex, the last position lists how many elements are adjacent to vertex. 

    # convert indexing from matlab to python (ISSM is built on matlab)
    vertexconnectivity[:,:100] = vertexconnectivity[:,:100] - 1
    elements2d = deepcopy(md.mesh.elements2d)
    elements2d = elements2d - 1

    # convert vertexconnectivity to list adjacent vertices instead of elements
    adjacentvertices = np.zeros(md.mesh.numberofvertices2d,dtype=object)
    for vertex in np.arange(md.mesh.numberofvertices2d): # loop through all mesh vertices
        vertexlist = np.zeros((vertexconnectivity[vertex,100]*3),dtype=int)
        count = 0
        for i in np.arange(vertexconnectivity[vertex,100]): # loop through all adjacent elements (vertexconnectivity[vertex,i])
            for j in np.arange(3): # loop through all three vertices of each element and store their indices in vertexlist
                vertexlist[count] = elements2d[vertexconnectivity[vertex,i],j]
                count = count + 1
        vertexlist = np.where(vertexlist == vertex, -1, vertexlist) # replace current vertex in vertexlist with -1
        vertexlist = np.unique(vertexlist) # crop vertexlist, right now there are lots of double listings, first position will be listing -1
        vertexlist[0] = vertexlist.size - 1 # replace -1 in first position with the number of adjacent vertices
        adjacentvertices[vertex] = vertexlist
        # !!! vertexlist is a list of the python indices of the vertices, not the original matlab vertex-numbers

    print('    Calculating Adjacent Vertices: Done')

    return vertexconnectivity, elements2d, adjacentvertices
    # Note: all variables that are created in the following subroutines have the same structure as adjacentvertices, but they drop the first position, i.e., don't list the total number of vertices

# }}}





def new_phi(LoadedModel, HydroPot, AdjacentVertices): #{{{
    # HydroPot is the basal hydropotential prescribed by you
    
    # some edge-case handling for BasalRouting.py
    # if a vertex shows a local high in hydropotential it won't receive any area because it routes to all adjacent vertices
    # conversely for local lows, it will receive area from everywhere and not route anything
    # both cases lead to discontinuities in the flow which is a probem for the routing of isotope values in IsotopeRouitng.py
    # smooth out by interpolation

    md = deepcopy(LoadedModel)
    phi = deepcopy(HydroPot)
    adjacentvertices = deepcopy(AdjacentVertices)
    
    for j in np.arange(md.mesh.numberofvertices2d):
        for i in adjacentvertices[j][1:]:
            if phi[j] > phi[i]:
                localhigh = 1
            else:
                localhigh = 0
                break
        if not localhigh:
            for i in adjacentvertices[j][1:]:
                if phi[j] < phi[i]:
                    locallow = 1
                else:
                    locallow = 0
                    break
        if localhigh or locallow:
            phi[j] = np.mean(phi[adjacentvertices[j][1:]]) # remember that the first entry in adjacentvertices is the total number of adjacent vertices which must be skipped here

    return phi

#}}}





def normals(LoadedModel, HydroPot, Elements2D): #{{{
    # Calculation of element normal vectors: norm = (norm_x, norm_y, norm_z)
    # important for next step
    
    md = deepcopy(LoadedModel)
    phi = deepcopy(HydroPot)
    elements2d = deepcopy(Elements2D)
    
    norm = np.zeros((md.mesh.numberofelements2d, 3))

    for i in np.arange(md.mesh.numberofelements2d):
        z1 = phi[elements2d[i,1]] - phi[elements2d[i,0]]
        z2 = phi[elements2d[i,2]] - phi[elements2d[i,0]]

        x1 = md.mesh.x2d[elements2d[i,1]] - md.mesh.x2d[elements2d[i,0]]
        x2 = md.mesh.x2d[elements2d[i,2]] - md.mesh.x2d[elements2d[i,0]]

        y1 = md.mesh.y2d[elements2d[i,1]] - md.mesh.y2d[elements2d[i,0]]
        y2 = md.mesh.y2d[elements2d[i,2]] - md.mesh.y2d[elements2d[i,0]]
        
        norm[i,0] = (z1 * y2) - (z2 * y1)
        norm[i,1] = (z1 * x2) - (z2 * x1)
        norm[i,2] = (y1 * x2) - (y2 * x1)

        if norm[i,2] < 0:
            # if norm_z < 0,i.e., if the normal is pointing down -> flip that around
            norm[i] = -1 * norm[i]

    print('    Calculation Element Normals: Done')

    return norm

# }}}





def elementorientation(LoadedModel, Norm): #{{{
    # Calculation of normal-derived direction d of each element
    # i.e., determine the horizontal component of each element normal
    # Seibert & McGlynn 2007, doi:10.1029/2006WR005128

    md = deepcopy(LoadedModel)
    norm = deepcopy(Norm)
        
    elementorientation = np.zeros((md.mesh.numberofelements2d))

    for i in np.arange(md.mesh.numberofelements2d):
        if norm[i,0] > 0 and norm[i,1] >= 0:
            elementorientation[i] = np.arctan(norm[i,1] / norm[i,0])
        elif norm[i,0] > 0 and norm[i,1] < 0:
            elementorientation[i] = np.arctan(norm[i,1] / norm[i,0]) + 2 * np.pi
        elif norm[i,0] < 0:
            elementorientation[i] = np.arctan(norm[i,1] / norm[i,0]) + np.pi
        elif norm[i,0] == 0 and norm[i,1] > 0:
            elementorientation[i] = np.pi/2
        elif norm[i,0] == 0 and norm[i,1] < 0:
            elementorientation[i] = 3 * np.pi/2

    print('    Calculation Steepest Slope Direction: Done')

    return elementorientation

#}}}






def vertexvectors(LoadedModel, AdjacentVertices): #{{{
    # Calculate the horizontal component of the mesh-edge vectors from vertex to each neighbouring vertex
    # this will be needed for the partitioning of flow

    md = deepcopy(LoadedModel)
    adjacentvertices = deepcopy(AdjacentVertices)

    vertexvectors = np.zeros(md.mesh.numberofvertices2d,dtype=object)
    
    for vertex in np.arange(md.mesh.numberofvertices2d):
        vectorlist = np.zeros((adjacentvertices[vertex][0].astype(int)))
        for i in np.arange(adjacentvertices[vertex][0].astype(int)):
            # Position 0 in adjacentvertices[:] gives amount of adjacent vertices
            # but the vertex indices are at Position [1:], i.e., np.arange(amount+1)
            # Therefore, iterate through i+1
            x = md.mesh.x2d[adjacentvertices[vertex][i+1].astype(int)] - md.mesh.x2d[vertex]
            y =	md.mesh.y2d[adjacentvertices[vertex][i+1].astype(int)] - md.mesh.y2d[vertex]

            if x > 0 and y >= 0:
                direction = np.arctan(y/x)
            elif x > 0 and y < 0:
                direction = np.arctan(y/x) + 2 * np.pi
            elif x < 0:
                direction = np.arctan(y/x) + np.pi
            elif x == 0 and y > 0:
                direction = np.pi/2
            elif x == 0 and y < 0:
                direction = 3 * np.pi/2

            direction = 2 * np.pi - direction # for some reason it's necessary to flip these upside-down...
            vectorlist[i] = direction

        vertexvectors[vertex] = vectorlist
        # vertexvectors has the same structure as adjacentvertices, but the secondary arrays are one index shorter, as they don't list the amount of adjacent vertices
        # instead of the indices of the adjacent vertices it gives the horiz. direction from the central vertex to these vertices
        
    print('    Calculation Edge-Vertex Directions: Done')

    return vertexvectors

#}}}






def elementflowdirections(LoadedModel, VertexConnectivity, Elements2D, AdjacentVertices, ElementOrientation, VertexVectors, HydroPot): #{{{
    # for each mesh element surrounding the current vertex, find out whether the element routes water
    #     - away or towards the centre vertex,
    #     - to one of its two neighboring elements, or
    #     - away from this region of the mesh.

    # Consider a triangular mesh element with vertices M, P1, and P2; M being the current mesh vertex, P1 and P2 the far vertices
    # Find out if elementorientation of this element, i.e., the direction of its downslope gradient with origin at M, is directed in between P1 and P2
    # If elementorientation is in between P1 and P2, elementorientation will be stored in elementflowdirections as the direction of the steepest slope of this element
    # If elementorientation is not in between P1 and P2, whichever edge ([M,P1] or [M,P2]) is closer to elementorientation will be set as the steepest slope of this element

    md = deepcopy(LoadedModel)
    vertexconnectivity = deepcopy(VertexConnectivity)
    elements2d = deepcopy(Elements2D)
    adjacentvertices = deepcopy(AdjacentVertices)
    elementorientation = deepcopy(ElementOrientation)
    vertexvectors = deepcopy(VertexVectors)
    phi = deepcopy(HydroPot)

    elementflowdirections = np.zeros(md.mesh.numberofvertices2d,dtype=object)
    directionflags = np.zeros(md.mesh.numberofvertices2d,dtype=object)

    for vertex in np.arange(md.mesh.numberofvertices2d):
        slopes = np.zeros(vertexconnectivity[vertex,100])
        flags = np.zeros(vertexconnectivity[vertex,100],dtype=int)
        for i in np.arange(vertexconnectivity[vertex,100]):
            element = vertexconnectivity[vertex,i]
            elementvertices = elements2d[element]
            P = np.squeeze(elementvertices[np.squeeze(np.where(elementvertices != vertex))])
            M = vertex
            P1 = P[0]
            P2 = P[1]


            if vertexvectors[M][np.squeeze(np.where(adjacentvertices[M][1:]==P[0]))] == elementorientation[element]:
                # Is elementorientation equal to edge [M,P1]?
                slopes[i] = elementorientation[element]

            elif vertexvectors[M][np.squeeze(np.where(adjacentvertices[M][1:]==P[1]))] == elementorientation[element]:
                # Is elementorientation equal to edge [M,P2]?
                slopes[i] = elementorientation[element]

            else:

                # The following is a geometric investigation of the directions based on the radian circel interval [0,2pi] with M being in the centre of the circle
                # We need to figure out if elementorientation is inside the small angle in between (M,P1) and (M,P2)
                # Let's check every possible configuration of flow directions individually!
                
                dP1 = vertexvectors[M][np.squeeze(np.where(adjacentvertices[M][1:]==P[0]))] # == d(M,P1) == A, i.e., the edge [M,P1]
                dP2 = vertexvectors[M][np.squeeze(np.where(adjacentvertices[M][1:]==P[1]))] # == d(M,P2) == C, i.e., the edge [M,P2]
                dEl = elementorientation[element] # == d == B

                # Find out if B is in between A and C
                # If it's not, then find out if it's closer to A or C
                
                # from the prior if clauses we already know that B is not equal to A or C
                # and since A and C span a triangle, so the angle between them is never equal Pi, as that would imply a 180° angle in one corner
                
                if dP2 > dP1:
                    # C is left of A
                    # Note: left meaning counterclockwise on the radian circle, right meaning clockwise on the radian circle
                    # Note: We don't know in which geometric order P1 and P2 are.
                    # The element triangle could be either [M,P1,P2] or [M,P2,P1], so we need to check both cases.
                    if (dP2 - dP1) < np.pi:
                        #small angle to the left of A
                        if dEl < dP1 or dEl > dP2:
                            # B is outside
                            # large angle contains 0/2pi crossing
                            # Note: The radian circle is a circle, i.e., it loops back on itself such that 0 equal 2pi
                            steepedge = steepestedge(md.mesh.x2d, md.mesh.y2d, M, P1, P2, phi)
                            if steepedge == "dP1":
                                slopes[i] = dP1
                            else:
                                slopes[i] = dP2

                        else:
                            # B is inside
                            slopes[i] = elementorientation[element]
                            flags[i] = 1
                    else:
                         #small angle to the right of A
                        if dEl < dP1 or dEl > dP2:
                            # B is inside
                            slopes[i] = elementorientation[element]
                            flags[i] = 1
                        else:
                            # B is outside
                            # large angle is contiuous (i.e., does not contain 0/2pi crossing)
                            steepedge = steepestedge(md.mesh.x2d, md.mesh.y2d, M, P1, P2, phi)
                            if steepedge == "dP1":
                                slopes[i] = dP1
                            else:
                                slopes[i] = dP2

                else:
                    # A is left of C
                    if (dP1 - dP2) < np.pi:
                        # small angle to the left of C
                        if dEl < dP2 or dEl > dP1:
                            # B is outisde
                            # large angle contains the 0/2pi crossing
                            steepedge = steepestedge(md.mesh.x2d, md.mesh.y2d, M, P1, P2, phi)
                            if steepedge == "dP1":
                                slopes[i] = dP1
                            else:
                                slopes[i] = dP2

                        else:
                            # B is inside
                            slopes[i] = elementorientation[element]
                            flags[i] = 1
                    else:
                        # small angle to the right of C
                        if dEl < dP2 or dEl > dP1:
                            # B is inside
                            slopes[i] = elementorientation[element]
                            flags[i] = 1
                        else:
                            # B is outside
                            # large angle is continuous
                            steepedge = steepestedge(md.mesh.x2d, md.mesh.y2d, M, P1, P2, phi)
                            if steepedge == "dP1":
                                slopes[i] = dP1
                            else:
                                slopes[i] = dP2

                # the case in which B is outside and exactly in the middle between A and C is equal to P1,P2 > M
                

        elementflowdirections[vertex] = slopes # the steepest slopes for each element, either one of the edges or in bewteen the edges
        directionflags[vertex] = flags # a flag that marks when the slope direction is in between the edges 

    print('    Calculation Element Flow Directions Ver. 2: Done')
    
    # route edge-flow to steepest edge, not nearest edge
    return elementflowdirections, directionflags

#}}}






def steepestedge(XCoor, YCoor, M, P1, P2, phi): #{{{
    # which edge is the steepest slope, [M,P1] or [M,P2]?
    # needed above in elementflowdirections()

    dz1 = phi[P1] - phi[M]
    dz2 = phi[P2] - phi[M]
    dy1 = YCoor[P1] - YCoor[M]
    dy2 = YCoor[P2] - YCoor[M]
    dx1 = XCoor[P1] - XCoor[M]
    dx2 = XCoor[P2] - XCoor[M]

    horlength1 = np.sqrt(np.square(dx1)+np.square(dy1))
    horlength2 = np.sqrt(np.square(dx2)+np.square(dy2))

    slope1 = -dz1 / horlength1
    slope2 = -dz2 / horlength2

    if slope1 > slope2: # slopes are going down, but they have positive values so the larger slope is the steeper slope
        steepedge = "dP1"
    else:
        steepedge = "dP2"

    return steepedge

#}}}






def steepestslopes(LoadedModel, ElementFlowDirections, DirectionFlags, AdjacentVertices, VertexVectors, HydroPot): #{{{
    # In elementflowdirections() we determined the downslope direction for each element adjacent to vertex
    # Now, we filter out the overall steepest downslope directions (can be more than one!)
    # Check for each element:
    #    1) if the elementflowdirection is in between P1 and P2 ("in-element direction"), save as a steepest slope for vertex
    #    2) if the elementflowdirection is one of the element edges
    #        - if the adjacent element has the same edge as its elementflowdirection, then save as a steepest slope for vertex
    #        - otherwise discard

    md = deepcopy(LoadedModel)
    elementflowdirections = deepcopy(ElementFlowDirections)
    directionflags = deepcopy(DirectionFlags) # flags in-element directions determined in elementflowdirections()
    adjacentvertices = deepcopy(AdjacentVertices)
    vertexvectors = deepcopy(VertexVectors)
    phi = deepcopy(HydroPot)
    
    steepestdownflowdirections = np.zeros(md.mesh.numberofvertices2d,dtype=object)
    flowtoedge = np.zeros(md.mesh.numberofvertices2d,dtype=object)

    for vertex in np.arange(md.mesh.numberofvertices2d):
        alldirections = elementflowdirections[vertex]
        allflags = directionflags[vertex]

        steepestdirections = np.zeros(1)

        for i in np.arange(alldirections.size):
            if allflags[i] == 1: # i.e., if alldirections[i] is an in-element direction
                steepestdirections = np.append(steepestdirections, alldirections[i])
            elif allflags[i] == 0 and alldirections[i] in np.delete(alldirections, i): # if alldirections[i] exists more than once then two elements point to that edge and it will be considered as an edge-flow
                vertexnumber = adjacentvertices[vertex][np.squeeze(np.where(vertexvectors[vertex] == alldirections[i]))+1]
                # +1 because adjacentvertices has the number of vertices at position [0]
                if phi[vertex] > phi[vertexnumber]: # double-check that the edge-flow is not pointing uphill
                    steepestdirections = np.append(steepestdirections, alldirections[i])

        steepestdirections = np.unique(steepestdirections)
        steepestdirections[0] = steepestdirections.size - 1
        steepestdownflowdirections[vertex] = steepestdirections
        
        edgeflags = np.zeros(steepestdirections.size)
        for i in np.arange(1,steepestdirections.size):
            edgeflags[i] = np.any(vertexvectors[vertex] == steepestdirections[i])
        flowtoedge[vertex] = edgeflags

    print('    Filter Steepest Directions: Done')

    return steepestdownflowdirections, flowtoedge # steepestdownflowdirections has size of array at position [0]

#}}}






def directionslopes(LoadedModel, SteepestDownflowDirections, AdjacentVertices, VertexVectors, HydroPot, VertexConnectivity, ElementFlowDirections, FlowToEdge, Norm): #{{{
    # Now that we know the steepest downflow directions at each vertex, we need to partition the flow between the different directions
    # Assign weight to the downslope vertices based on steepnes of the slope
    # 1) Calculate Slope for the downslope directions (slope based on hydropotential, not topography!)

    md = deepcopy(LoadedModel)
    steepestdownflowdirections = deepcopy(SteepestDownflowDirections)
    adjacentvertices = deepcopy(AdjacentVertices)
    vertexvectors = deepcopy(VertexVectors)
    phi = deepcopy(HydroPot)
    vertexconnectivity = deepcopy(VertexConnectivity)
    elementflowdirections = deepcopy(ElementFlowDirections)
    flowtoedge = deepcopy(FlowToEdge)
    norm = deepcopy(Norm)
    
    directionslope = np.zeros(md.mesh.numberofvertices2d,dtype=object)
    vertexidentifier = np.zeros(md.mesh.numberofvertices2d,dtype=object)
    for vertex in np.arange(md.mesh.numberofvertices2d):
        slopes = np.zeros(steepestdownflowdirections[vertex].size)
        numberofelement = np.zeros(steepestdownflowdirections[vertex].size)

        if steepestdownflowdirections[vertex][0] == 0: # Is there no steepst downslope direction at vertex?
            if md.mesh.vertexonboundary[vertex]: # Maybe vertex is at the domain boundary?
                directionslope[vertex] = 'boundary'
                vertexidentifier[vertex] = 'boundary'
            else:
                directionslope[vertex] = 'sink' # Maybe vertex is a local sink? This shouldn't actually occur as it's the edge-case dealt with in new_phi()
                vertexidentifier[vertex] = 'sink'
                
        else:
            for i in np.arange(1,steepestdownflowdirections[vertex].size):

                # if the downflow direction is an "edge-flow":
                if flowtoedge[vertex][i]:
                    # calculate slope from vertices of that edge
                    vertexnumber = adjacentvertices[vertex][np.squeeze(np.where(vertexvectors[vertex] == steepestdownflowdirections[vertex][i]))+1]
                    # +1 because adjvert has the number of vertices at position [0]
                    dz = phi[vertexnumber] - phi[vertex]
                    dy = md.mesh.y2d[vertexnumber] - md.mesh.y2d[vertex]
                    dx = md.mesh.x2d[vertexnumber] - md.mesh.x2d[vertex]
                    horlength = np.sqrt(np.square(dx)+np.square(dy))
                    slopes[i] = -dz / horlength # the slope is defined as the tan() of the slope angle
                    numberofelement[i] = vertexnumber

                # if the downflow direction is an "in-element flow":
                else:
                    # calculate slope from element normal
                    element = vertexconnectivity[vertex,np.where(elementflowdirections[vertex] == steepestdownflowdirections[vertex][i])]
                    # the normal is perpendicular to the element surface: dz -> dhor, horlength -> vertlength
                    dhor = np.sqrt(np.square(norm[element,0])+np.square(norm[element,1]))
                    vertlength = norm[element,2]
                    slopes[i] = dhor / vertlength
                    numberofelement[i] = element

            directionslope[vertex] = slopes
            vertexidentifier[vertex] = numberofelement

    print('    Calculation Slope: Done')

    return directionslope, vertexidentifier

#}}}






def directionweight(LoadedModel, SteepestDownflowDirections, DirectionSlope): #{{{
    # 2) Calculate the weight of each downflow direction (Quinn et al, 1991, https://doi.org/10.1002/hyp.3360050106)

    md = deepcopy(LoadedModel)
    steepestdownflowdirections = deepcopy(SteepestDownflowDirections)
    directionslope = deepcopy(DirectionSlope)
    
    directionweight = np.zeros(md.mesh.numberofvertices2d,dtype=object)
    for vertex in np.arange(md.mesh.numberofvertices2d):
        weights = np.zeros(steepestdownflowdirections[vertex].size)
        for i in np.arange(1,steepestdownflowdirections[vertex].size):
            if type(directionslope[vertex]) != str and steepestdownflowdirections[vertex].size > 1:
                # Option 1:
#                weights[i] = directionslope[vertex][i] / np.sum(directionslope[vertex][1:])
                # Option 2, following Freeman, 1991, https://doi.org/10.1016/0098-3004(91)90048-I:
                p = 1.1
                weights[i] = np.power(directionslope[vertex][i], p) / np.sum(np.power(directionslope[vertex][1:], p))

        directionweight[vertex] = weights # weights has a 0 at position [0], where steepestdownfl. has the array size

    return directionweight

#}}}
    




    
def receivingvertices(LoadedModel, DirectionWeight, VertexIdentifier, FlowToEdge, Elements2D, VertexVectors, AdjacentVertices, SteepestDownflowDirections, HydroPot): #{{{
    # 3) Identify the adjacent vertices that vertex routes to and determine the weighting for each flow direction 

    # vertexidentifier contains the indices of the vertices (for edge-flow) and elements (for in-element flow) that the downslope directions point to
    # flowtoedge flags 1 if vertexidentifier refers to a vertex, 0 if to an element

    md = deepcopy(LoadedModel)
    directionweight = deepcopy(DirectionWeight)
    vertexidentifier = deepcopy(VertexIdentifier)
    flowtoedge = deepcopy(FlowToEdge)
    elements2d = deepcopy(Elements2D)
    vertexvectors = deepcopy(VertexVectors)
    adjacentvertices = deepcopy(AdjacentVertices)
    steepestdownflowdirections = deepcopy(SteepestDownflowDirections)
    phi = deepcopy(HydroPot)
    
    receivingvertices = np.zeros(md.mesh.numberofvertices2d,dtype=object)
    vertexweights = np.zeros(md.mesh.numberofvertices2d,dtype=object)
    
    for vertex in np.arange(md.mesh.numberofvertices2d):
        weights = np.zeros(directionweight[vertex].size - 1,dtype=object)
        vertices = np.zeros(directionweight[vertex].size - 1,dtype=object)
        
        for i in np.arange(directionweight[vertex].size - 1):

            # edge-flow: 
            if flowtoedge[vertex][i+1]:
                vertices[i] = vertexidentifier[vertex][i+1]
                weights[i] = directionweight[vertex][i+1]

            # in-element flow:
            # the routing needs to be split between P1 and P2 (Tarboton, 1997, https://doi.org/10.1029/96WR03137)
            else:
                elementvertices = elements2d[vertexidentifier[vertex][i+1].astype(int)]
                P1P2 = -1 * np.ones(2)
                k = 0
                for j in np.arange(3):
                    if elementvertices[j] != vertex:
                        P1P2[k] = elementvertices[j]
                        k = k + 1
                vertices[i] = P1P2
                weights[i] = np.ones(2)
                vP1 = vertexvectors[vertex][np.squeeze(np.where(adjacentvertices[vertex][1:]==P1P2[0]))]
                vP2 = vertexvectors[vertex][np.squeeze(np.where(adjacentvertices[vertex][1:]==P1P2[1]))]
                angleP1P2 = vP1 - vP2

                if angleP1P2 < 0:
                    angleP1P2 = vP2 - vP1
                if angleP1P2 > np.pi:
                    angleP1P2 = 2 * np.pi - angleP1P2
                    
                vflow = steepestdownflowdirections[vertex][i+1]
                angleP1flow = vP1 - vflow
                
                if angleP1flow < 0:
                    angleP1flow = vflow - vP1
                if angleP1flow > np.pi:
                    angleP1flow = 2 * np.pi - angleP1flow
                    
                angleflowP2 = vflow - vP2
                
                if angleflowP2 < 0:
                    angleflowP2 = vP2 - vflow
                if angleflowP2 > np.pi:
                    angleflowP2 = 2 * np.pi - angleflowP2
                    
                weights[i][0] = (angleflowP2 / angleP1P2) * directionweight[vertex][i+1]
                weights[i][1] = (angleP1flow / angleP1P2) * directionweight[vertex][i+1]

                # Consider this edge-case:
                # direction of [M,P1] is 2.76, direction of [M,P2] is 4.5, elementorientation is 2.77
                # i.e., elementorientation is in in-element direction, and the flow will be split between P1 and P2
                # However, in this case P2 is actually upstream of M, so if we split the flow like above, some flow will go upstream
                # Therefore:

                if phi[int(P1P2[0])] > phi[vertex]:
                    weights[i][0] = 0
                    weights[i][1] = 1 * directionweight[vertex][i+1]
                elif phi[int(P1P2[1])] > phi[vertex]:
                    weights[i][0] = 1 * directionweight[vertex][i+1]
                    weights[i][1] = 0

        # flatten out vertices and weights, they are a mess of nested arrays and have all sorts of dimensions...
        verticessorted = np.empty(0)
        weightssorted = np.empty(0)

        if vertices.shape[0] == 0:
            verticessorted = verticessorted
            weightssorted = weightssorted
        elif vertices.shape[0] == 1:
            if vertices[0].shape == ():
                verticessorted = np.append(verticessorted, vertices[0])
                weightssorted = np.append(weightssorted, weights[0])
            else:
                for f in np.arange(vertices[0].shape[0]):
                    verticessorted = np.append(verticessorted, vertices[0][f])
                    weightssorted = np.append(weightssorted, weights[0][f])
        elif vertices.shape[0] >= 2:
            for g in np.arange(vertices.shape[0]):
                if vertices[g].shape == ():
                    verticessorted = np.append(verticessorted, vertices[g])
                    weightssorted = np.append(weightssorted, weights[g])
                else:
                    for h in np.arange(vertices[g].shape[0]):
                        verticessorted = np.append(verticessorted, vertices[g][h])
                        weightssorted = np.append(weightssorted, weights[g][h])

        # remove double entries, recalculate weights for twice-appointed vertices
        finalvertexlist = np.empty(0)
        finalweights = np.empty(0)
        indexofdouble = np.empty(0)
        if verticessorted.size == 0:
            pass
        else:
            for e in np.arange(verticessorted.size):
                if np.any(indexofdouble == e):
                    pass
                elif np.any(np.delete(verticessorted,e) == verticessorted[e]):
                    indexofdouble = np.append(indexofdouble, e)
                    double = np.where(np.delete(verticessorted,e) == verticessorted[e])[0][0] + 1
                    indexofdouble = np.append(indexofdouble, double)

                    finalvertexlist = np.append(finalvertexlist, verticessorted[e])
                    finalweights = np.append(finalweights, (weightssorted[e] + weightssorted[double]))
                else:
                    finalvertexlist = np.append(finalvertexlist, verticessorted[e])
                    finalweights = np.append(finalweights, weightssorted[e])

        # remove the 0-weights from the case above
        if np.any(finalweights==0):
            thisone = np.where(finalweights==0)
            finalvertexlist = np.delete(finalvertexlist, thisone)
            finalweights = np.delete(finalweights, thisone)
            
        receivingvertices[vertex] = finalvertexlist
        vertexweights[vertex] = finalweights

    return receivingvertices, vertexweights
    # receivingvertices are the adjacent vertices that receive flow from vertex
    # vertexweights are the corresponding weights of how much flow they receive from vertex
# }}}




        # That's it :)
