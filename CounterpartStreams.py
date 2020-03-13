# -*- coding: cp1251 -*-
# Counterpart streams
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

# borrowed from https://gis.stackexchange.com/questions/150200/reversing-polyline-direction-based-on-raster-value-using-arcpy
def FlipLine(Line):
    rPnts = arcpy.Array()
    for i in range(len(Line)):
        rPnts.append(Line[len(Line) - i - 1]) # flip the points in the array
    OutShape = arcpy.Polyline(rPnts)
    return OutShape

def get_values(features, field):
    return numpy.asarray([row[0] for row in arcpy.da.SearchCursor(features, field)])

def set_values(features, field, values):
    with arcpy.da.UpdateCursor(features, field) as rows:
        i = 0
        for row in rows:
            row[0] = values[i]
            rows.updateRow(row)
            i += 1
    return features

def process_raster(instreams, inIDfield, in_raster, minacc, radius, deviation, demraster, penalty, startpts, endpts,
                   ids, ordids, ordnears, lowerleft, cellsize, crs, outstreams):

    try:
        global MAXACC

        rrast = arcpy.sa.Raster(in_raster)
        arcpy.env.extent = rrast.extent  # Very important!
        arcpy.env.snapRaster = rrast  # Very important!

        inraster = arcpy.RasterToNumPyArray(in_raster, nodata_to_value=MAXACC + 1)

        ni = inraster.shape[0]
        nj = inraster.shape[1]
        n = len(ids)

        minx = lowerleft.X
        miny = lowerleft.Y

        eucs = numpy.zeros((len(ids), inraster.shape[0], inraster.shape[1])).astype(float)

        instreamslyr = 'strlyr'
        arcpy.MakeFeatureLayer_management(instreams, instreamslyr)

        arcpy.AddMessage('CALCULATING EUCLIDEAN DISTANCE RASTERS...')
        for i in range(n):
            arcpy.SelectLayerByAttribute_management(instreamslyr, 'NEW_SELECTION',
                                                    '"' + inIDfield + '" = ' + str(ordids[i]))
            euc = arcpy.sa.EucDistance(instreamslyr, cell_size=cellsize)
            eucs[i, :, :] = arcpy.RasterToNumPyArray(euc)

        idx = numpy.zeros(n).astype(int)

        for i in range(n):
            idx[i] = numpy.where(ids == ordids[i])[0]

        startxy = []
        for row in arcpy.da.SearchCursor(startpts, ["SHAPE@XY"]):
            x, y = row[0]
            startxy.append([x, y])

        endxy = []
        for row in arcpy.da.SearchCursor(endpts, ["SHAPE@XY"]):
            x, y = row[0]
            endxy.append([x, y])

        startxy = [startxy[i] for i in idx]
        endxy = [endxy[i] for i in idx]

        extinraster = extend_array(inraster, 1, 1, 0)

        arcpy.AddMessage("TRACING COUNTERPARTS...")

        arcpy.SetProgressor("step", "Processing rivers", 0, n - 1, 1)

        streams = []
        types = []

        for k in range(n):

            iend = ni - math.trunc((endxy[k][1] - miny) / cellsize) - 1
            istart = ni - math.trunc((startxy[k][1] - miny) / cellsize) - 1

            jend = math.trunc((endxy[k][0] - minx) / cellsize)
            jstart = math.trunc((startxy[k][0] - minx) / cellsize)

            endneigh = get_neighborhood(iend, jend, radius, cellsize, ni, nj)
            startneigh = get_neighborhood(istart, jstart, radius, cellsize, ni, nj)

            depid = ordnears[k]
            isdep = False
            dep = []
            if depid != -1:
                isdep = True
                dep = streams[numpy.where(ordids == depid)[0].tolist()[0]]
                for cell in endneigh:
                    if cell in dep:
                        endneigh = get_neighborhood(cell[0], cell[1], radius, cellsize, ni, nj)
                        break

            arcpy.AddMessage("ID = " + str(ordids[k]) + ' (' + str(k + 1) + " from " + str(n) + ')')

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
                        D = euc_distance((i, j), startneigh[0]) / L
                        w = math.sqrt(E + 1) * (D + 1)
                        if w < weight:
                            stream = s
                            weight = w
            extend = True

            if len(stream) > 0:

                if isdep:
                    nl = len(stream)

                    for i in range(nl):
                        if (stream[i] in dep):
                            nl = i + 1
                            extend = False
                            break

                    if not extend:
                        streams.append(stream[:nl])
                        types.append('Stream')

                else:
                    streams.append(stream)
                    types.append('Stream')
                    extend = False

            else:
                extend = False

            if (len(stream) == 0) or extend:

                npstart = numpy.full((ni, nj), -1).astype(int)

                if extend:
                    npstart[stream[-1]] = ordids[k]
                    arcpy.AddMessage("Extending by shortest path")
                else:
                    npstart[startneigh[0]] = ordids[k]
                    arcpy.AddMessage("Using shortest path")

                startlyr = arcpy.NumPyArrayToRaster(npstart, lowerleft, cellsize, value_to_nodata=-1)
                arcpy.DefineProjection_management(startlyr, crs)

                npend = numpy.full((ni, nj), -1).astype(int)
                npend[endneigh[0]] = ordids[k]
                endlyr = arcpy.NumPyArrayToRaster(npend, lowerleft, cellsize, value_to_nodata=-1)
                arcpy.DefineProjection_management(endlyr, crs)

                euc = arcpy.NumPyArrayToRaster(eucs[k, :, :], lowerleft, cellsize)
                arcpy.DefineProjection_management(euc, crs)

                euc_mask = (euc + 1) * arcpy.sa.Reclassify(euc, "value",
                                                           arcpy.sa.RemapRange([[0, deviation, penalty],
                                                                                [deviation, euc.maximum,
                                                                                 'NODATA']]))

                strs = arcpy.sa.Reclassify(in_raster, "value",
                                           arcpy.sa.RemapRange([[0, minacc, 'NODATA'], [minacc, MAXACC, 1]]))

                cost = euc * arcpy.sa.ExtractByMask(strs, euc_mask)

                arcpy.Mosaic_management(euc_mask * arcpy.sa.Raster(demraster), cost, 'MINIMUM')

                backlink = arcpy.sa.CostBackLink(startlyr, cost)
                costpath = arcpy.sa.CostPath(endlyr, cost, backlink)
                costdist = arcpy.sa.CostDistance(startlyr, cost)

                nppath = arcpy.RasterToNumPyArray(costpath, nodata_to_value=-1)
                npdist = arcpy.RasterToNumPyArray(costdist, nodata_to_value=-1)

                cellidx = numpy.where(nppath >= 0)
                cells = numpy.argwhere(nppath >= 0)

                values = npdist[cellidx]

                idx = numpy.argsort(values)

                path = list(map(tuple, cells[idx, :]))

                if isdep:
                    nl = len(path)

                    for i in range(nl):
                        if (path[i] in dep):
                            nl = i + 1
                            break
                    if extend:
                        streams.append(stream + path[1:nl])
                        types.append('Stream + Path')
                    else:
                        streams.append(path[:nl])
                        types.append('Path')
                else:
                    streams.append(path)
                    types.append('Path')

            arcpy.SetProgressorPosition(k)

        outraster = None
        nodatavalue = 0 if (min(ordids) > 0) else min(ordids) - 1

        N = len(streams)

        if (N > 0):
            outraster = numpy.full((N, ni, nj), nodatavalue)
            for l in range(N):
                for ncells in range(len(streams[l])):
                    outraster[l, streams[l][ncells][0], streams[l][ncells][1]] = ordids[l]

        outinnerraster = arcpy.sa.Int(arcpy.NumPyArrayToRaster(outraster[0, :, :],
                                                               lowerleft, cellsize, value_to_nodata=nodatavalue))
        arcpy.DefineProjection_management(outinnerraster, crs)

        result = 'in_memory/result'
        arcpy.RasterToPolyline_conversion(outinnerraster, result, background_value='NODATA', simplify='NO_SIMPLIFY')

        templines = 'in_memory/templine'
        for k in range(1, N):
            outinnerraster = arcpy.sa.Int(arcpy.NumPyArrayToRaster(outraster[k, :, :],
                                                                   lowerleft, cellsize, value_to_nodata=nodatavalue))
            arcpy.DefineProjection_management(outinnerraster, crs)
            arcpy.RasterToPolyline_conversion(outinnerraster, templines, background_value='NODATA',
                                              simplify='NO_SIMPLIFY')
            arcpy.Append_management(templines, result, schema_type='NO_TEST')

        arcpy.AddField_management(result, 'type', 'TEXT', field_length=16)
        result = set_values(result, 'type', types)

        arcpy.UnsplitLine_management(result, outstreams, ['grid_code', 'type'])

        # ensure right direction
        arcpy.AddMessage('ENSURING RIGHT DIRECTION...')

        with  arcpy.da.UpdateCursor(outstreams, ["SHAPE@", 'grid_code']) as rows:
            for row in rows:
                line = row[0].getPart(0)
                count_start = [line[0].X, line[0].Y]
                id = row[1]

                idx = numpy.where(ordids == id)[0].tolist()[0]

                hydro_start = startxy[idx]
                hydro_end = endxy[idx]

                if euc_distance(count_start, hydro_start) > euc_distance(count_start, hydro_end):
                    row[0] = FlipLine(line)
                    rows.updateRow(row)

        arcpy.Densify_edit(outstreams, 'DISTANCE', cellsize)

        return

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        raise Exception

def execute(in_streams, inIDfield, inraster, demRaster, outstreams, minacc, penalty, radius, deviation):
    global MAXACC
    MAXACC = float(str(arcpy.GetRasterProperties_management(inraster, "MAXIMUM")))
    desc = arcpy.Describe(inraster)
    lowerleft = arcpy.Point(desc.extent.XMin, desc.extent.YMin)
    cellsize = desc.meanCellWidth
    crs = desc.spatialReference

    domain = 'in_memory/domain'
    arcpy.RasterDomain_3d(inraster, domain, 'POLYGON')

    instreams_crop = 'in_memory/str_cropped'
    arcpy.Clip_analysis(in_streams, domain, instreams_crop)


    # Get start and endpoints of rivers
    startpts = 'in_memory/startpts'
    endpts = 'in_memory/endpts'

    arcpy.FeatureVerticesToPoints_management(instreams_crop, startpts, point_location = 'START')
    arcpy.FeatureVerticesToPoints_management(instreams_crop, endpts, point_location = 'END')

    # Calculate distance to determine hierarchy
    dist_lines = 'in_memory/dist_lines'

    arcpy.GenerateNearTable_analysis(endpts, instreams_crop, dist_lines, closest=False, closest_count=2)
    ids = get_values(instreams_crop, inIDfield)
    ins = get_values(dist_lines, 'IN_FID')
    nears = get_values(dist_lines, 'NEAR_FID')
    dist = get_values(dist_lines, 'NEAR_DIST')

    # dependent are streams which endpoints are located
    # less or equal to cellsize from another

    ordids = None
    ordnears = None

    if len(ids) == 1:
        ordids = ids
        ordnears = numpy.asarray([-1])
    else:
        flt = numpy.logical_and(ins != nears, dist <= cellsize)

        depids = ids[ins[flt] - 1]
        depnears = ids[nears[flt] - 1]

        # construct the order
        flt2 = numpy.logical_and(ins != nears, dist > cellsize)
        ordids = ids[ins[flt2] - 1]
        ordnears = numpy.full(len(ordids), -1).astype(int) # TODO: make this value more robust

        arcpy.AddMessage('INDEPENDENT STREAMS: ' + str(ordids))

        while(len(depnears) > 0):
            ord = numpy.logical_not(numpy.in1d(depnears, depids))
            ordids = numpy.append(ordids, depids[ord])
            ordnears = numpy.append(ordnears, depnears[ord])
            depids = depids[numpy.logical_not(ord)]
            depnears = depnears[numpy.logical_not(ord)]

    # Tracing stream lines

    process_raster(instreams_crop, inIDfield, inraster, minacc, radius, deviation, demRaster, penalty,
                   startpts, endpts, ids, ordids, ordnears, lowerleft, cellsize, crs, outstreams)

    return

if __name__ == "__main__":
    try:
        inStreams = arcpy.GetParameterAsText(0)
        inIDfield = arcpy.GetParameterAsText(1)
        inRaster = arcpy.GetParameterAsText(2)
        demRaster = arcpy.GetParameterAsText(3)
        outStreams = arcpy.GetParameterAsText(4)
        minAcc = float(arcpy.GetParameterAsText(5))
        penalty = int(arcpy.GetParameterAsText(6))
        radius = float(arcpy.GetParameterAsText(7))
        deviation = float(arcpy.GetParameterAsText(8))

        execute(inStreams, inIDfield, inRaster, demRaster, outStreams, minAcc, penalty, radius, deviation)

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)