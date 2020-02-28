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

def frechet_dist(P,Q):
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
    return ca[n-1, m-1]

def euc_matrix(P, Q):
    mdist = cdist(P, Q, 'euclidean')
    return mdist

def execute(in_hydrolines, hydro_field, in_counterparts, count_field, out_links, out_area):

    arcpy.CreateFeatureclass_management(os.path.dirname(out_links), os.path.basename(out_links),
                                        geometry_type='POLYLINE', spatial_reference=in_hydrolines)

    arcpy.AddField_management(out_links, 'ID', 'LONG')
    arcpy.AddField_management(out_links, 'DIR', 'TEXT', field_length=8)

    insertcursor = arcpy.da.InsertCursor(out_links, ['SHAPE@', 'ID', 'DIR'])

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

        eucs = euc_matrix(count_coords, hydro_coords)

        ni = len(count_coords)
        nj = len(hydro_coords)

        # find basic min j for each i
        minjays = []
        for i in range(ni):
            minj = numpy.argmin(eucs[i, :])
            for k in range (i, ni):
                curj = numpy.argmin(eucs[k, :])
                if  minj > curj:
                    minj = curj
            minjays.append(minj)

        # fill empty min j by connecting to nearest i
        jbacks = []
        ibacks = []
        curj = 0
        for i in range(1, ni-1):
            nextj = minjays[i]
            if nextj - curj > 1:
                for j in range(curj+1, nextj):
                    iback = numpy.argmin(eucs[(i-1):(i+1), j]) + i - 1
                    jbacks.append(j)
                    ibacks.append(iback)
            curj = nextj

        # check if the last points are connected
        if nj-1 not in minjays:
            for j in range(curj + 1, nj):
                jbacks.append(j)
                ibacks.append(ni-1)

        # check if the first points are connected
        if 0 not in minjays:
            ibacks.insert(0, 0)
            jbacks.insert(0, 0)

        pairs = zip(range(ni), minjays)
        backpairs = zip(ibacks, jbacks)

        # pairs = [[0, 0]]
        #
        # i = 0
        # j = 0
        #
        # np = 1
        #
        # while i < ni - 1 or j < nj - 1:
        #     if j == nj - 1:
        #         pairs.append([i + 1, j])
        #         i += 1
        #     elif i == ni - 1:
        #         pairs.append([i, j + 1])
        #         j += 1
        #     elif eucs[i + 1, j] < eucs[i, j + 1] and eucs[i + 1, j] < eucs[i + 1, j + 1]:
        #         if np > 1:
        #             if pairs[-2][0] == pairs[-1][0] and pairs[-2][1] != pairs[-1][1]:
        #                 pairs = pairs[:-1]
        #                 np -= 1
        #                 pairs.append([i + 1, j + 1])
        #                 i += 1
        #                 j += 1
        #             else:
        #                 pairs.append([i + 1, j])
        #                 i += 1
        #         else:
        #             pairs.append([i + 1, j])
        #             i += 1
        #     elif eucs[i, j + 1] < eucs[i + 1, j] and eucs[i, j + 1] < eucs[i + 1, j + 1]:
        #         if np > 1:
        #             if pairs[-2][0] != pairs[-1][0] and pairs[-2][1] == pairs[-1][1]:
        #                 pairs = pairs[:-1]
        #                 np -= 1
        #                 pairs.append([i + 1, j + 1])
        #                 i += 1
        #                 j += 1
        #             else:
        #                 pairs.append([i, j + 1])
        #                 j += 1
        #         else:
        #             pairs.append([i, j + 1])
        #             j += 1
        #     else:
        #         pairs.append([i + 1, j + 1])
        #         i += 1
        #         j += 1
        #     np += 1

        fdist = frechet_dist(count_coords, hydro_coords)
        arcpy.AddMessage('Frechet distance: ' + str(fdist) + ' (' + hydro_field + ' = ' + str(id) + ')')

        features = []
        for pair in pairs:
            line = [arcpy.Point(*count_coords[pair[0]]), arcpy.Point(*hydro_coords[pair[1]])]
            features.append(arcpy.Polyline(arcpy.Array(line)))

        backfeatures = []
        for pair in backpairs:
            line = [arcpy.Point(*count_coords[pair[0]]), arcpy.Point(*hydro_coords[pair[1]])]
            backfeatures.append(arcpy.Polyline(arcpy.Array(line)))

        for feature in features:
            insertcursor.insertRow([feature, id, 'Forward'])

        for backfeature in backfeatures:
            insertcursor.insertRow([backfeature, id, 'Backward'])

    if out_area is not None:
        polys = 'in_memory/polys'
        arcpy.FeatureToPolygon_management([in_hydrolines, in_counterparts, out_links], polys)
        arcpy.Dissolve_management(polys, out_area)

    return

if __name__ == 'main':
    in_hydrolines = arcpy.GetParameterAsText(0)
    hydro_field = arcpy.GetParameterAsText(1)
    in_counterparts = int(arcpy.GetParameterAsText(2))
    count_field = int(arcpy.GetParameterAsText(3))
    out_links = arcpy.GetParameterAsText(4)
    out_area = arcpy.GetParameterAsText(5)

    try:
        execute(in_hydrolines, hydro_field, in_counterparts, count_field, out_links, out_area)

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        print("Processing failed")