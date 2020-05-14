# -*- coding: cp1251 -*-
import arcpy
import sys
import traceback
import math
import numpy
from arcpy.sa import *
from datetime import datetime

shift = []
window = []
noData = -9999
# 2019, Timofey Samsonov, Lomonosov Moscow State University


def process_raster(inraster, niter, nfilt):
    global noData

    ni = inraster.shape[0]
    nj = inraster.shape[1]

    arcpy.SetProgressor("step", "Processing rows", 0, ni-1, 1)
    raster = inraster

    outraster = numpy.zeros((ni, nj))

    return outraster

def euc_distance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

def path_length(path):
    L = 0
    n = len(path)
    (ic, jc) = path[0]
    for k in range(1, n):
        (i, j) = path[k]
        L += euc_distance((i, j), (ic, jc))
        ic = i
        jc = j
    return L


def execute(in_raster, in_streams, in_field, out_raster):

    EPS = 1

    inraster = arcpy.sa.Raster(in_raster)
    arcpy.env.extent = inraster.extent  # Very important!
    arcpy.env.snapRaster = inraster  # Very important!
    cell_size = float(arcpy.GetRasterProperties_management(in_raster, "CELLSIZEX").getOutput(0).replace(',', '.'))
    ncols = int(arcpy.GetRasterProperties_management(in_raster, "COLUMNCOUNT").getOutput(0))
    nrows = int(arcpy.GetRasterProperties_management(in_raster, "ROWCOUNT").getOutput(0))

    lowerleft = arcpy.Point(inraster.extent.XMin, inraster.extent.YMin)
    crs = inraster.spatialReference

    stream = 'in_memory/stream'
    spnt = 'in_memory/spnt'

    features = []
    ids = []
    strlyr = 'strlyr'
    arcpy.MakeFeatureLayer_management(in_streams, strlyr)

    npdem = arcpy.RasterToNumPyArray(in_raster, ncols = ncols, nrows = nrows)

    arcpy.AddMessage('RASTERIZING STREAMS' + str(datetime.now()))

    with arcpy.da.SearchCursor(in_streams, in_field) as rows:
        for row in rows:
            id = row[0]
            ids.append(id)
            arcpy.AddMessage('ID = ' + str(id))
            arcpy.SelectLayerByAttribute_management(strlyr, 'NEW_SELECTION', in_field + ' = ' + str(id))
            arcpy.FeatureVerticesToPoints_management(strlyr, spnt, 'START')

            arcpy.PolylineToRaster_conversion(strlyr, in_field, stream, cellsize = cell_size)
            str_mask = Power(stream, 0)

            cost = CostDistance(spnt, str_mask)
            npcost = arcpy.RasterToNumPyArray(cost, nodata_to_value = -1)

            cells = numpy.argwhere(npcost >= 0)
            values = npcost[npcost >= 0]

            idx = numpy.argsort(values)

            feature = cells[idx, :]

            features.append(feature)

    arcpy.AddMessage('CARVING' + str(datetime.now()))

    finished = False

    ifeat = 0

    while not finished:
        for feature in features:
            ni = len(feature)
            cell = feature[0]
            zprev = npdem[cell[0], cell[1]]
            ishill = False

            carved = [cell]

            ncarved = 0

            for i in range(1, ni):
                cell = feature[i]
                z = npdem[cell[0], cell[1]]

                if z >= zprev:
                    ishill = True

                    if i == ni - 1:
                        z = zprev - EPS
                        npdem[cell[0], cell[1]] = z
                    else:
                        carved.append(cell)

                if z < zprev:
                    if ishill:
                        ncarved += 1

                        ishill = False
                        carved.append(cell)
                        carlength = path_length(carved)

                        ncar = len(carved)
                        dz = z - zprev

                        for j in range(1, ncar-1):
                            npdem[carved[j][0], carved[j][1]] = zprev + dz * path_length(carved[0:j]) / carlength

                    carved = [cell]
                    zprev = z

            arcpy.AddMessage('ID = ' + str(ids[ifeat]) + ': carved ' + str(ncarved) + ' sections')
            ifeat += 1
        finished = True

    outraster = arcpy.NumPyArrayToRaster(npdem, lowerleft, cell_size)
    arcpy.DefineProjection_management(outraster, crs)

    outraster.save(out_raster)

    arcpy.AddMessage('END' + str(datetime.now()))

    return

if __name__ == "__main__":
    try:
        in_raster = arcpy.GetParameterAsText(0)
        in_streams = arcpy.GetParameterAsText(1)
        in_field = arcpy.GetParameterAsText(2)
        out_raster = arcpy.GetParameterAsText(3)

        execute(in_raster, in_streams, in_field, out_raster)
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)