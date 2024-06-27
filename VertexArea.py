
# VertexArea.py

# Calculate Vertex Area
# Necessary for BasalRouting.py (flux is basal melt times area)

# The area of each vertex is defined to be the area of the polygon described by
    # 1) the centroids of the connected elements and
    # 2) the midpoints of the edges of the elements
    
# The intersection of each connected element with this polygon creates an irregular quadrilateral inside each element
# The polygon-area is the sum of these quadrilaterals

# The quadrilaterals themselves can be broken down into two triangles with known side-lengths (mesh-coordinates are known)
# The triangle area can be calculated with Heron's formula sqrt(p*(p-a)*(p-b)+(p-c))
# with a,b,c being the lengths of the triangle-sides, p is the perimeter 0.5*(a+b+c)
# The coordinates of the centroid can be calculated as xc = 1/3 * (x1+x2+x3) etc

# There is probably an easier method using vectors or whatever, so feel free to do that instead
# If your model is not ISSM, it's probably easier to write your own script for this step

# --------------


# python packages
import numpy as np
import math as m
from importlib import reload
from copy import deepcopy

# ISSM packages
from model import * 




def vertexarea(LoadedModel): # LoadedModel is your ISSM model

        # variables
    
    md = deepcopy(LoadedModel)
    vertexconnectivity = NodeConnectivity(md.mesh.elements2d, md.mesh.numberofvertices2d)
    # vertexconnectivity is a 2d array of shape (numberofvertices2d, 100).
    # It lists the indices of all adjacent elements at each vertex. Position (:,100) is the number of adjacent elements at each vertex.
    elements2d = deepcopy(md.mesh.elements2d)
    # elements2d are the indices of the 2d mesh-elements

    # transform Matlab indexing to python indexing
    vertexconnectivity[:,:100] = vertexconnectivity[:,:100] - 1
    elements2d = elements2d - 1


        # data arrays
    
    vertexarea = np.zeros(md.mesh.numberofvertices2d)
    centroid = np.zeros(5)
    halfv2 = np.zeros(5)
    halfv3 = np.zeros(5)


        # compute vertexarea
        # at each vertex, we cycle through the adjacent elements constructing the two triangles inside each element and adding all triangles together
    
    for vertex in np.arange(md.mesh.numberofvertices2d):
        for element in vertexconnectivity[vertex,:(vertexconnectivity[vertex,100])]:
            v1 = vertex
            v23 = np.squeeze(elements2d[element,np.where(elements2d[element]!=vertex)])
            v2 = v23[0]
            v3 = v23[1]

            # centroid coordinates:
            centroid[0] = 1/3 * (md.mesh.x2d[v1] + md.mesh.x2d[v2] + md.mesh.x2d[v3])
            centroid[1] = 1/3 * (md.mesh.y2d[v1] + md.mesh.y2d[v2] + md.mesh.y2d[v3])
            centroid[2] = 1/3 * (md.mesh.z[v1] + md.mesh.z[v2] + md.mesh.z[v3])
            # distance vertex-centroid:
            centroid[3] = np.sqrt(np.square(centroid[0] - md.mesh.x2d[v1]) + np.square(centroid[1] - md.mesh.y2d[v1]) + np.square(centroid[2] - md.mesh.z[v1]))

            # element - side midpoint coordinates:
            halfv2[0] = 1/2 * (md.mesh.x2d[v1] + md.mesh.x2d[v2])
            halfv2[1] = 1/2 * (md.mesh.y2d[v1] + md.mesh.y2d[v2])
            halfv2[2] = 1/2 * (md.mesh.z[v1] + md.mesh.z[v2])
            halfv3[0] = 1/2 * (md.mesh.x2d[v1] + md.mesh.x2d[v3])
            halfv3[1] = 1/2 * (md.mesh.y2d[v1] + md.mesh.y2d[v3])
            halfv3[2] = 1/2 * (md.mesh.z[v1] + md.mesh.z[v3])
            
            # distance vertex - side midpoints:
            halfv2[3] = np.sqrt(np.square(halfv2[0] - md.mesh.x2d[v1]) + np.square(halfv2[1] - md.mesh.y2d[v1]) + np.square(halfv2[2] - md.mesh.z[v1]))
            halfv3[3] = np.sqrt(np.square(halfv3[0] - md.mesh.x2d[v1]) + np.square(halfv3[1] - md.mesh.y2d[v1]) + np.square(halfv3[2] - md.mesh.z[v1]))

            # distance side midpoints - centroid
            halfv2[4] = np.sqrt(np.square(halfv2[0] - centroid[0]) + np.square(halfv2[1] - centroid[1]) + np.square(halfv2[2] - centroid[2]))
            halfv3[4] = np.sqrt(np.square(halfv3[0] - centroid[0]) + np.square(halfv3[1] - centroid[1]) + np.square(halfv3[2] - centroid[2]))

            # triangle perimeters:
            pv2 = 0.5 * (centroid[3] + halfv2[3] + halfv2[4])
            pv3 = 0.5 *	(centroid[3] + halfv3[3] + halfv3[4])

            # triangle areas, Heron's method:
            areav2 = np.sqrt(pv2 * (pv2 - centroid[3]) * (pv2 - halfv2[3]) * (pv2 - halfv2[4]))
            areav3 = np.sqrt(pv3 * (pv3 - centroid[3]) * (pv3 - halfv3[3]) * (pv3 - halfv3[4]))

            vertexarea[vertex] = vertexarea[vertex] + areav2 + areav3

    print('    Calculation Vertexarea: Done')
    return vertexarea

