# -*- coding: cp1251 -*-
# Raster stream network generalization by Leonowicz-Jenny algorithm
# 2020, Timofey Samsonov, Lomonosov Moscow State University
import sys
import arcpy
import numpy
import traceback
import math

MAXACC = 0

def find_up_cell(accraster, i, j):
    w = [[0.70710678, 1, 0.70710678],[1, 1, 1], [0.70710678, 1, 0.70710678]]  # distance weights
    shift = [-1, 0, 1]
    minmax = 4000000000
    a = 0
    kmin = 1
    lmin = 1

    # finding differences in 3x3 neighbourhood

    for k in shift:
        b = 0
        ik = i+k
        for l in shift:
            jl = j+l
            temp = (accraster[i, j] - accraster[ik, jl]) * w[a][b]
            if 0 < temp < minmax:
                minmax = temp 
                kmin = a
                lmin = b      
            b += 1
        a += 1

    iUp = i + shift[kmin]
    jUp = j + shift[lmin]

    return iUp, jUp


def find_down_cell(accraster, i, j):
    w = [[0.70710678, 1, 0.70710678],[1, 1, 1], [0.70710678, 1, 0.70710678]]  # distance weights
    shift = [-1, 0, 1]
    minmax = 0
    a = 0
    kmin = 1
    lmin = 1

    # finding differences in 3x3 neighbourhood

    for k in shift:
        b = 0
        ik = i+k
        for l in shift:
            jl = j+l
            temp = (accraster[ik, jl] - accraster[i, j]) * w[a][b]
            if temp > minmax:
                minmax = temp
                kmin = a
                lmin = b
            b += 1
        a += 1

    iUp = i + shift[kmin]
    jUp = j + shift[lmin]

    return iUp, jUp

def find_cell(accraster, i, j, Down = True):
    if Down:
        return find_down_cell(accraster, i, j)
    else:
        return find_up_cell(accraster, i, j)

def trace_flow_cells(accraster, i, j, minacc, neigh, down = True):
    acc = accraster[i, j]
    ik = i
    jk = j
    n = 0
    stream = []
    selcells = []
    seln = []
    inside = False

    try:
        if acc >= minacc:
            while True:
                current = (ik, jk)

                if current in neigh:
                    selcells.append(neigh.index(current))
                    seln.append(n)
                    inside = True
                elif inside: # we previously get into the neighborhood
                    break

                stream.append(current)

                inext, jnext = find_cell(accraster, ik, jk, down)

                if inext == ik and jnext == jk:
                    break

                ik = inext
                jk = jnext
                n += 1

            if len(selcells) > 0:
                end = seln[numpy.argsort(selcells)[0]] + 1
                return stream[0:end]
            else: return []
        else: return []

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        raise Exception

def extend_array(array, nx, ny, value):
    
    ni = array.shape[0]
    nj = array.shape[1]

    extarray = numpy.empty((ni+ny, nj+nx))
    extarray.fill(value)
    for i in range(ni):
        for j in range(nj):
            extarray[i, j] = array[i, j]

    return extarray

def get_neighborhood(i, j, radius, cellsize, ni, nj):
    w = int(math.ceil(radius / cellsize))  # calculate kernel radius (rounded)
    l = range(-w, w + 1)  # calculate kernel indices
    idx = numpy.meshgrid(l, l)  # generate coordinate matrices

    flt = (idx[0] ** 2 + idx[1] ** 2 <= w ** 2)  # filter by distance

    x = idx[0][flt] + i
    y = idx[1][flt] + j

    flt_xy = (x >= 0) * (x < ni) * (y >= 0) * (y < nj)  # filter by domain

    x = x[flt_xy]
    y = y[flt_xy]

    order = numpy.argsort(((x - i) ** 2 + (y - j) ** 2) ** 0.5)

    neigh = list(map(lambda a, b: (a, b), x[order], y[order]))

    return neigh

def process_raster(inraster, minacc, radius, startxy, endxy, minx, miny, cellsize):

    try:
        global MAXACC
        ni = inraster.shape[0]
        nj = inraster.shape[1]
        n = len(startxy)

        outraster = numpy.zeros((ni, nj))
        extinraster = extend_array(inraster, 1, 1, 0)

        arcpy.SetProgressor("step", "Processing rivers", 0, n - 1, 1)

        streams = []
        for k in range(n):

            iend = ni - math.trunc((endxy[k][1] - miny) / cellsize) - 1
            istart = ni - math.trunc((startxy[k][1] - miny) / cellsize) - 1

            jend = math.trunc((endxy[k][0] - minx) / cellsize)
            jstart = math.trunc((startxy[k][0] - minx) / cellsize)

            endneigh = get_neighborhood(iend, jend, radius, cellsize, ni, nj)
            startneigh = get_neighborhood(istart, jstart, radius, cellsize, ni, nj)

            arcpy.SetProgressorLabel("Finding closest stream for river " + str(k) + " from " + str(n))

            stream = []
            for (i, j) in startneigh:
                if  minacc < inraster[i, j] > minacc:
                    stream = trace_flow_cells(extinraster, i, j, minacc, endneigh)
                    if len(stream) > 0:
                        break
            streams.append(stream)

            arcpy.SetProgressorPosition(k)

        streams.sort(key = len)

        for k in range(n):
            for l in range(len(streams[k])):
                outraster[streams[k][l][0], streams[k][l][1]] = k

        return outraster
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        raise Exception

def execute(instreams, inraster, outstreams, minacc, radius):
    global MAXACC
    MAXACC = float(str(arcpy.GetRasterProperties_management(inraster, "MAXIMUM")))
    desc = arcpy.Describe(inraster)
    lowerleft = arcpy.Point(desc.extent.XMin, desc.extent.YMin)
    cellsize = desc.meanCellWidth
    crs = desc.spatialReference


    rasternumpy = arcpy.RasterToNumPyArray(inraster, nodata_to_value = MAXACC + 1)

    startpts = 'in_memory/startpts'
    endpts = 'in_memory/endpts'

    arcpy.FeatureVerticesToPoints_management(instreams, startpts, point_location='START')
    startxy = []
    for row in arcpy.da.SearchCursor(startpts, ["SHAPE@XY"]):
        x, y = row[0]
        startxy.append([x, y])


    arcpy.FeatureVerticesToPoints_management(instreams, endpts, point_location = 'END')
    endxy = []
    for row in arcpy.da.SearchCursor(endpts, ["SHAPE@XY"]):
        x, y = row[0]
        endxy.append([x, y])

    # Tracing stream lines
    arcpy.AddMessage("Searching for closest stream lines...")
    newrasternumpy = process_raster(rasternumpy, minacc, radius, startxy, endxy, lowerleft.X, lowerleft.Y, cellsize)

    # Convert python list to ASCII
    arcpy.AddMessage("Writing streams...")
    outinnerraster = arcpy.sa.Int(arcpy.NumPyArrayToRaster(newrasternumpy, lowerleft, cellsize, value_to_nodata=0))
    arcpy.DefineProjection_management(outinnerraster, crs)
    # outinnerraster.save(outraster)

    arcpy.RasterToPolyline_conversion(outinnerraster, outstreams, simplify='NO_SIMPLIFY')

if __name__ == "__main__":
    try:
        inStreams = arcpy.GetParameterAsText(0)
        inRaster = arcpy.GetParameterAsText(1)
        outStreams = arcpy.GetParameterAsText(2)
        minAcc = float(arcpy.GetParameterAsText(3))
        radius = float(arcpy.GetParameterAsText(4))

        execute(inStreams, inRaster, outStreams, minAcc, radius)
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)