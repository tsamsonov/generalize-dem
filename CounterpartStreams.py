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

def extend_array(array, nx, ny, value):

    ni = array.shape[0]
    nj = array.shape[1]

    extarray = numpy.empty((ni + ny, nj + nx))
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

def cell_distance(cell1, cell2):
    return math.sqrt((cell1[0] - cell2[0])**2 + (cell1[1] - cell2[1])**2)

def path_length(path):
    L = 0
    n = len(path)
    (ic, jc) = path[0]
    for k in range(1, n):
        (i, j) = path[k]
        L += cell_distance((i,j), (ic,jc))
        ic = i
        jc = j
    return L


def trace_flow_cells(accraster, euc, i, j, minacc, endneigh, down = True):
    acc = accraster[i, j]
    ik = i
    jk = j
    n = 0
    stream = []
    endcells = []
    endn = []
    e = []
    in_end = False

    try:
        if acc >= minacc:
            while True:
                current = (ik, jk)

                if current in endneigh:
                    endcells.append(endneigh.index(current))
                    endn.append(n)
                    in_end = True
                elif in_end: # we previously get into the neighborhood
                    break

                stream.append(current)
                e.append(euc[ik, jk])

                inext, jnext = find_cell(accraster, ik, jk, down)

                if inext == ik and jnext == jk:
                    break

                ik = inext
                jk = jnext
                n += 1

            if len(endcells) > 0:
                end = endn[numpy.argsort(endcells)[0]] + 1
                return stream[0:end], e[0:end]
            else: return [], []
        else: return [], []

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        raise Exception

def process_raster(inraster, eucs, minacc, radius, deviation, startxy, endxy, ids, minx, miny, cellsize):

    try:
        global MAXACC
        ni = inraster.shape[0]
        nj = inraster.shape[1]
        n = len(startxy)

        extinraster = extend_array(inraster, 1, 1, 0)

        arcpy.SetProgressor("step", "Processing rivers", 0, n - 1, 1)

        streams = []
        failed = []
        succeeded = []
        for k in range(n):

            iend = ni - math.trunc((endxy[k][1] - miny) / cellsize) - 1
            istart = ni - math.trunc((startxy[k][1] - miny) / cellsize) - 1

            jend = math.trunc((endxy[k][0] - minx) / cellsize)
            jstart = math.trunc((startxy[k][0] - minx) / cellsize)

            endneigh = get_neighborhood(iend, jend, radius, cellsize, ni, nj)
            startneigh = get_neighborhood(istart, jstart, radius, cellsize, ni, nj)

            arcpy.AddMessage("Finding closest stream for river " + str(k) + " from " + str(n))

            stream = []
            weight = float('Inf')
            for (i, j) in startneigh:
                if  inraster[i, j] > minacc:
                    s, e = trace_flow_cells(extinraster, eucs[k,:,:], i, j, minacc, endneigh)
                    ncells = len(e)

                    if ncells > 0 and max(e) <= deviation:
                        l = 0
                        cur = s[l]
                        startstream = []
                        while cur in startneigh:
                            startstream.append(cur)
                            l += 1
                            if (l < ncells):
                                cur = s[l]
                            else:
                                break

                        nstart = l-1

                        L = path_length(startstream) + 1

                        E = sum(e[0:nstart]) / (cellsize * L)
                        D = cell_distance((i, j), startneigh[0]) / L
                        w = math.sqrt(E + 1) * (D + 1)
                        if w < weight:
                            stream = s
                            weight = w

            if len(stream) == 0:
                failed.append(ids[k])
            else:
                streams.append(stream)
                succeeded.append(ids[k])

            arcpy.SetProgressorPosition(k)

        outraster = None
        nodatavalue = 0 if (min(ids) > 0) else min(ids) - 1
        N = len(streams)

        if (N > 0):
            outraster = numpy.full((N, ni, nj), nodatavalue)

            for l in range(N):
                for ncells in range(len(streams[l])):
                    outraster[l, streams[l][ncells][0], streams[l][ncells][1]] = succeeded[l]

        return outraster, nodatavalue, failed

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        raise Exception

def execute(in_streams, inIDfield, inraster, outstreams, minacc, radius, deviation):
    global MAXACC
    MAXACC = float(str(arcpy.GetRasterProperties_management(inraster, "MAXIMUM")))
    desc = arcpy.Describe(inraster)
    lowerleft = arcpy.Point(desc.extent.XMin, desc.extent.YMin)
    cellsize = desc.meanCellWidth
    crs = desc.spatialReference

    rasternumpy = arcpy.RasterToNumPyArray(inraster, nodata_to_value = MAXACC + 1)

    domain = 'in_memory/domain'
    arcpy.RasterDomain_3d(inraster, domain, 'POLYGON')

    instreams_crop = 'in_memory/str_cropped'
    arcpy.Clip_analysis(in_streams, domain, instreams_crop)

    startpts = 'in_memory/startpts'
    endpts = 'in_memory/endpts'


    arcpy.FeatureVerticesToPoints_management(instreams_crop, startpts, point_location='START')
    startxy = []
    startID = []
    for row in arcpy.da.SearchCursor(startpts, ["SHAPE@XY", inIDfield]):
        x, y = row[0]
        startxy.append([x, y])
        startID.append(row[1])


    arcpy.FeatureVerticesToPoints_management(instreams_crop, endpts, point_location = 'END')
    endxy = []
    endID = []
    for row in arcpy.da.SearchCursor(endpts, ["SHAPE@XY", inIDfield]):
        x, y = row[0]
        endxy.append([x, y])
        endID.append(row[1])

    eucs = numpy.zeros((len(startxy), rasternumpy.shape[0], rasternumpy.shape[1])).astype(float)

    i = 0
    rrast = arcpy.sa.Raster(inraster)
    arcpy.env.extent = rrast.extent  # Very important!
    arcpy.env.snapRaster = rrast  # Very important!
    instreamslyr = 'strlyr'
    arcpy.MakeFeatureLayer_management(instreams_crop, instreamslyr)
    ids = []
    for row in arcpy.SearchCursor(instreams_crop):
        id = row.getValue(inIDfield)
        ids.append(id)
        arcpy.SelectLayerByAttribute_management(instreamslyr, 'NEW_SELECTION', '"' + inIDfield + '" = ' + str(id))
        euc = arcpy.sa.EucDistance(instreamslyr, cell_size=cellsize)
        eucs[i, :, :] = arcpy.RasterToNumPyArray(euc)
        i += 1

    arcpy.AddMessage(startID)
    arcpy.AddMessage(endID)
    arcpy.AddMessage(ids)

    # return

    # Tracing stream lines
    arcpy.AddMessage("Searching for closest stream lines...")
    newrasternumpy, nodatavalue, failed = process_raster(rasternumpy, eucs, minacc, radius, deviation, startxy,
                                                         endxy, ids, lowerleft.X, lowerleft.Y, cellsize)

    nfailed = len(failed)
    nsucc = i - nfailed

    result = 'in_memory/result'
    if nsucc == 0:
        arcpy.CreateFeatureclass_management('in_memory', 'result', 'POLYLINE', spatial_reference=crs)
        arcpy.AddField_management(result, 'type', 'TEXT', field_length=10)
    else:
        # Convert python list to ASCII
        outinnerraster = arcpy.sa.Int(arcpy.NumPyArrayToRaster(newrasternumpy[0, :, :],
                                                               lowerleft, cellsize, value_to_nodata=nodatavalue))
        arcpy.DefineProjection_management(outinnerraster, crs)
        arcpy.RasterToPolyline_conversion(outinnerraster, result, background_value='NODATA', simplify='NO_SIMPLIFY')
        arcpy.AddField_management(result, 'type', 'TEXT', field_length=10)

        templines = 'in_memory/templine'
        for k in range(1, nsucc):
            outinnerraster = arcpy.sa.Int(arcpy.NumPyArrayToRaster(newrasternumpy[k, :, :],
                                                                   lowerleft, cellsize, value_to_nodata=nodatavalue))
            arcpy.DefineProjection_management(outinnerraster, crs)
            arcpy.RasterToPolyline_conversion(outinnerraster, templines, background_value='NODATA', simplify='NO_SIMPLIFY')
            arcpy.Append_management(templines, result, schema_type = 'NO_TEST')

        arcpy.CalculateField_management(result, 'type', '"Real"')

    # Process failed streams using shortest path strategy
    if len(failed) > 0:

        arcpy.AddMessage("Searching for spurious counterpart streams. IDs: " + str(failed))

        startlyr = 'startlyr'
        arcpy.MakeFeatureLayer_management(startpts, startlyr)
        endlyr = 'endlyr'
        arcpy.MakeFeatureLayer_management(endpts, endlyr)
        i = 0
        for row in arcpy.SearchCursor(instreams_crop):
            id = row.getValue(inIDfield)
            if id in failed:

                arcpy.AddMessage("ID =  " + str(id))

                arcpy.SelectLayerByAttribute_management(instreamslyr, 'NEW_SELECTION',
                                                        '"' + inIDfield + '" = ' + str(id))
                arcpy.SelectLayerByAttribute_management(startlyr, 'NEW_SELECTION',
                                                        '"' + inIDfield + '" = ' + str(id))
                arcpy.SelectLayerByAttribute_management(endlyr, 'NEW_SELECTION',
                                                        '"' + inIDfield + '" = ' + str(id))
                euc = arcpy.NumPyArrayToRaster(eucs[i,:,:], lowerleft, cellsize)
                arcpy.DefineProjection_management(euc, crs)

                euc_mask = arcpy.sa.Reclassify(euc, "value", arcpy.sa.RemapRange([[0,deviation, 100000],[deviation,euc.maximum,'NODATA']]))

                strs = arcpy.sa.Reclassify(inraster, "value", arcpy.sa.RemapRange([[0,minacc,'NODATA'],[minacc,MAXACC,1]]))

                cost = arcpy.sa.ExtractByMask(strs, euc_mask) * euc

                arcpy.Mosaic_management(euc_mask, cost, 'MINIMUM')

                backlink = arcpy.sa.CostBackLink(startlyr, cost)
                costpath = arcpy.sa.CostPath(endlyr, cost, backlink)

                if arcpy.sa.IsNull(costpath).minimum == 0:
                    costlines = 'in_memory/costlines'
                    arcpy.RasterToPolyline_conversion(costpath + 1, costlines,
                                                      background_value='NODATA',
                                                      simplify='NO_SIMPLIFY')
                    arcpy.CalculateField_management(costlines, 'grid_code', id)
                    arcpy.Append_management(costlines, result, schema_type = 'NO_TEST')

            i += 1

        reslyr = 'reslyr'
        arcpy.MakeFeatureLayer_management(result, reslyr)
        arcpy.SelectLayerByAttribute_management(reslyr, 'NEW_SELECTION', 'type is null')
        arcpy.CalculateField_management(reslyr, 'type', '"Imitating"')

    arcpy.UnsplitLine_management(result, outstreams, ['grid_code', 'type'])
    arcpy.Densify_edit(outstreams, 'DISTANCE', cellsize)

if __name__ == "__main__":
    try:
        inStreams = arcpy.GetParameterAsText(0)
        inIDfield = arcpy.GetParameterAsText(1)
        inRaster = arcpy.GetParameterAsText(2)
        outStreams = arcpy.GetParameterAsText(3)
        minAcc = float(arcpy.GetParameterAsText(4))
        radius = float(arcpy.GetParameterAsText(5))
        deviation = float(arcpy.GetParameterAsText(6))

        execute(inStreams, inIDfield, inRaster, outStreams, minAcc, radius, deviation)

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)