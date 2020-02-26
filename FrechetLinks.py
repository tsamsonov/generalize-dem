import os
import sys
import arcpy
import traceback
import math
import numpy
from scipy.spatial.distance import cdist

# Frechet distance implementation is borrowed from
# https://gist.github.com/MaxBareiss/ba2f9441d9455b56fbc9
def euc_dist(p1, p2):
    return math.sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]))

def _c(ca, i, j, P, Q):
    if ca[i,j] > -1:
        return ca[i,j]
    elif i == 0 and j == 0:
        ca[i,j] = euc_dist(P[0],Q[0])
    elif i > 0 and j == 0:
        ca[i,j] = max(_c(ca,i-1,0,P,Q),euc_dist(P[i],Q[0]))
    elif i == 0 and j > 0:
        ca[i,j] = max(_c(ca,0,j-1,P,Q),euc_dist(P[0],Q[j]))
    elif i > 0 and j > 0:
        ca[i,j] = max(min(_c(ca,i-1,j,P,Q),_c(ca,i-1,j-1,P,Q),_c(ca,i,j-1,P,Q)),euc_dist(P[i],Q[j]))
    else:
        ca[i,j] = float("inf")
    return ca[i,j]

def frechet_dist(P,Q):
    ca = numpy.full((len(P),len(Q)), -1)

    return _c(ca, len(P)-1, len(Q)-1, P, Q)

def frechet_matrix(P,Q):
    n = len(P)
    m = len(Q)
    ca = numpy.full((n, m), -1)

    ca[0, 0] = euc_dist(P[0], Q[0])

    for i in range(1, n):
        ca[i, 0] = max(ca[i - 1, 0], euc_dist(P[i], Q[0]))
    for j in range(1, m):
        ca[0, j] = max(ca[0, j - 1], euc_dist(P[0], Q[j]))
    for i in range(1, n):
        for j in range(1, m):
            ca[i, j] = max(min(ca[i - 1, j],
                               ca[i, j - 1],
                               ca[i - 1, j - 1]),
                           euc_dist(P[i], Q[j]))
    return ca

def euc_matrix(P, Q):
    mdist = cdist(P, Q, 'euclidean')
    return mdist

def execute(in_hydrolines, hydro_field, in_counterparts, count_field):

    hcursor = arcpy.da.SearchCursor(in_hydrolines, ['SHAPE@', hydro_field])

    for row in hcursor:
        hydro_coords = []
        for pnt in row[0].getPart().next():
            hydro_coords.append([pnt.X, pnt.Y])

        id = row[1]
        ccursor = arcpy.da.SearchCursor(in_counterparts, ['SHAPE@', count_field], where_clause = count_field + ' = ' + str(id))
        counterpart = ccursor.next()

        count_coords = []
        for pnt in counterpart[0].getPart().next():
            count_coords.append([pnt.X, pnt.Y])

        # fmatrix = frechet_matrix(hydro_coords, count_coords)

        eucs = euc_matrix(hydro_coords, count_coords)

        nb_forward = numpy.argmin(eucs, 0)
        nb_back = numpy.argmin(eucs, 1)

        arcpy.AddMessage('FORWARD: ' + str(nb_forward))
        arcpy.AddMessage('BACK: ' + str(nb_back))

        # arcpy.AddMessage('Last euc distance: ' + str(eucs[-1, -1]))

    return

if __name__ == 'main':
    in_hydrolines = arcpy.GetParameterAsText(0)
    hydro_field = arcpy.GetParameterAsText(1)
    in_counterparts = int(arcpy.GetParameterAsText(2))
    count_field = int(arcpy.GetParameterAsText(3))

    try:
        execute(in_hydrolines, hydro_field, in_counterparts, count_field)

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        print("Processing failed")