# -*- coding: cp1251 -*-
# Counterpart streams
# 2020, Timofey Samsonov, Lomonosov Moscow State University
import sys
import arcpy
import numpy
import traceback
import math
import Utils

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

def get_window(ij, ni, nj, size=3):
    w = int((size - 1) / 2)
    l = range(-w, w + 1)  # calculate kernel indices
    idx = numpy.meshgrid(l, l)  # generate coordinate matrices

    x = idx[0] + ij[0]
    y = idx[1] + ij[1]

    flt_xy = (x >= 0) * (x < ni) * (y >= 0) * (y < nj) #* numpy.logical_or(x != ij[0], y != ij[1])

    x = x[flt_xy]
    y = y[flt_xy]

    npdist = ((x - ij[0]) ** 2 + (y - ij[1]) ** 2) ** 0.5

    order = numpy.argsort(npdist)

    dist = npdist[order].tolist()

    neigh = list(map(lambda a, b: (a, b), x[order], y[order]))

    return neigh[1:], dist[1:]

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

def get_coordinates(features):
    lines = []
    with arcpy.da.SearchCursor(features, "SHAPE@") as rows:
        for row in rows:
            coords = []
            for pnt in row[0].getPart().next():
                coords.append([pnt.X, pnt.Y])
            lines.append(coords)
    return lines

def set_values(features, field, values):
    with arcpy.da.UpdateCursor(features, field) as rows:
        i = 0
        for row in rows:
            row[0] = values[i]
            rows.updateRow(row)
            i += 1
    return features

def get_neighbor(path, npdist, npcomp, ni, nj, i, j, minimax, prohibited):

    win = get_window(i, j, ni, nj)

    comp = npcomp[i][j]
    curmin = min(comp)

    win = set(win) - set(path)

    d = float("inf")
    selected = None
    for cell in win:
        if cell not in prohibited:
            nextcomp = npcomp[cell[0]][cell[1]]
            nextmin = min(nextcomp)
            if (nextmin > curmin):
                continue

            if (minimax in nextcomp):
                dnext = npdist[cell]
                if (dnext < d):
                    d = dnext
                    selected = cell
    return selected

def invback(a, b):

    offset = (a[0] - b[0], a[1] - b[1])

    back = {
        6: (-1, -1),
        7: (-1, 0),
        8: (-1, 1),
        5: (0, -1),
        0: (0, 0),
        3: (0, 1),
        4: (1, -1),
        1: (1, 0),
        2: (1, 1),
    }

    for k in range(9):
        if back[k] == offset:
            return k

    return None

def cost_distance(source, npcost, npcomp, cellsize, nodatavalue=-1):
    visited = []
    calculated = [source]

    ni = npcost.shape[0]
    nj = npcost.shape[1]

    total = (npcost != nodatavalue).sum()

    arcpy.AddMessage('Need to visit ' + str(total) + ' cells')

    npdist = numpy.full((ni, nj), float('Inf'))
    npback = numpy.full((ni, nj), nodatavalue)

    npdist[source] = 0
    npback[source] = 0

    while len(calculated) > 0:
        newcalc = []
        newdist = []
        for cell in calculated:
            celldist = npdist[cell]
            nb, dist = get_window(cell, ni, nj)
            for (ij, d) in zip(nb, dist):
                if (ij not in visited) and (npcost[ij] != nodatavalue):
                    accum_dist = celldist + d * npcost[ij]
                    if accum_dist < npdist[ij]:
                        npdist[ij] = accum_dist
                        npback[ij] = invback(cell, ij)
                        if (ij not in newcalc):
                            newcalc.append(ij)
                            newdist.append(accum_dist)
                        else:
                            newdist[newcalc.index(ij)] = accum_dist

            visited.append(cell)
        calculated = [ij for dist, ij in sorted(zip(newdist, newcalc))]
        perc = 100 * float(len(visited)) / float(total)
        arcpy.AddMessage('Percentage processed: ' + str(perc)) # TODO: correct estimate!

    return npdist, npback


def cost_path(coords, cost, radius, source, destination):
    back = {
        0: (0,  0),
        1: (1,  0),
        2: (1,  1),
        3: (0,  1),
        4: (1, -1),
        5: (0, -1),
        6: (-1,-1),
        7: (-1, 0),
        8: (-1, 1)
    }

    desc = arcpy.Describe(cost)
    lowerleft = arcpy.Point(desc.extent.XMin, desc.extent.YMin)
    cellsize = desc.meanCellWidth
    minx = lowerleft.X
    miny = lowerleft.Y

    # npdist  = arcpy.RasterToNumPyArray(distance)
    # npback = arcpy.RasterToNumPyArray(backlink)
    npcost = arcpy.RasterToNumPyArray(cost, nodata_to_value = -1)


    ni = npcost.shape[0]
    nj = npcost.shape[1]

    # idest = ni - math.trunc((destination[1] - miny) / cellsize) - 1
    # jdest = math.trunc((destination[0] - minx) / cellsize)

    # Generate compatibility raster

    arcpy.AddMessage('Generating compatibility sets')
    npcomp = [[set() for j in range(nj)] for i in range(ni)] # compatibility raster

    arcpy.AddMessage('Empty sets created')

    k = 0
    for pnt in coords:
        ip = ni - math.trunc((pnt[1] - miny) / cellsize) - 1
        jp = math.trunc((pnt[0] - minx) / cellsize)
        nb = get_neighborhood(ip, jp, radius, cellsize, ni, nj)
        for (i, j) in nb:
            npcomp[i][j] = npcomp[i][j].union([k])
        k +=1

    arcpy.AddMessage('Calculating distance and backlink rasters')
    npdist, npback = cost_distance(source, npcost, npcomp, cellsize)

    ij = destination
    path = [ij]
    type = npback[ij]
    k = 0

    arcpy.AddMessage('Tracing')
    # minimax = [min(npcomp[ij[0]][ij[1]])]
    # kprob = None
    # minprob = None
    # prohibited = []
    # problem = False

    while type != 0:
        shift = back[npback[ij]]
        ij = (ij[0] + shift[0], ij[1] + shift[1])

        # ij = get_neighbor(path, npdist, npcomp, ni, nj, ij[0], ij[1], minimax[k], prohibited)

        # if ijback != ij:
        #     kprob = k - 1
        #
        # if ij == None:
        #     minprob = minimax[k-1]
        #     for l in range(kprob, k):
        #         prohibited.append(path.pop())
        #         minimax.pop()
        #     problem = True
        #     k = kprob
        #     ij = path[-1]
        #     continue

        path.append(ij)
        # nextcomp = npcomp[ij[0]][ij[1]]
        #
        # nextmini = minimax[-1]
        # while (nextmini - 1) in nextcomp:
        #     nextmini = nextmini - 1
        # minimax.append(nextmini)
        #
        # if (problem == True and minimax[-1] < minprob):
        #     problem = False
            # prohibited = []

        type = npback[ij]
        arcpy.AddMessage('Point added: ' + str(k))
        k += 1

    path.reverse()

    return path

def process_raster(instreams, inIDfield, in_raster, minacc, radius, deviation, demraster, penalty, startpts, endpts,
                   ids, ordids, ordends, ordstarts, lowerleft, cellsize, crs, outstreams):

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

        geometries = get_coordinates(instreams)
        if len(geometries) > 1:
            geometries = geometries[idx]

        arcpy.AddMessage(geometries)

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

            endid = ordends[k]
            enddep = False

            endstr = []
            if endid != -1:
                enddep = True
                endstr = streams[numpy.where(ordids == endid)[0].tolist()[0]]
                for cell in endneigh:
                    if cell in endstr:
                        endneigh = get_neighborhood(cell[0], cell[1], radius, cellsize, ni, nj)
                        break

            startid = ordstarts[k]
            startdep = False

            startstr = []
            if startid != -1:
                startdep = True
                startstr = streams[numpy.where(ordids == startid)[0].tolist()[0]]
                for cell in startneigh:
                    if cell in startstr:
                        startneigh = get_neighborhood(cell[0], cell[1], radius, cellsize, ni, nj)
                        break

            arcpy.AddMessage("ID = " + str(ordids[k]) + ' (' + str(k + 1) + " from " + str(n) + ')')

            stream = []
            extend = False

            if not startdep:
                weight = float('Inf')
                for (i, j) in startneigh:
                    if  inraster[i, j] > minacc:
                        s, e = trace_flow_cells(extinraster, eucs[k,:,:], i, j, minacc, endneigh)
                        ncells = len(e)

                        if False: #ncells > 0 and Utils.frechet_dist(s, geometries[k]) <= deviation:
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

                    if enddep:
                        nl = len(stream)

                        for i in range(nl):
                            if (stream[i] in endstr):
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

                # euc_mask = (euc + 1) * arcpy.sa.Reclassify(euc, "value",
                #                                            arcpy.sa.RemapRange([[0, deviation, penalty],
                #                                                                 [deviation, euc.maximum,
                #                                                                  'NODATA']]))

                euc_mask = arcpy.sa.Reclassify(euc, "value",
                                                           arcpy.sa.RemapRange([[0, deviation, penalty],
                                                                                [deviation, euc.maximum,
                                                                                 'NODATA']]))

                strs = arcpy.sa.Reclassify(in_raster, "value",
                                           arcpy.sa.RemapRange([[0, minacc, 'NODATA'], [minacc, MAXACC, 1]]))

                cost = arcpy.sa.ExtractByMask(strs, euc_mask) #* euc

                # arcpy.Mosaic_management(euc_mask * arcpy.sa.Raster(demraster), cost, 'MINIMUM')
                arcpy.Mosaic_management(euc_mask, cost, 'MINIMUM')

                # arcpy.env.workspace = "X:/DEMGEN/"
                # if arcpy.Exists("cost.tif"):
                #     arcpy.Delete_management("cost.tif")
                #
                # cost.save('cost.tif')

                # backlink = arcpy.sa.CostBackLink(startlyr, cost)
                # costdist = arcpy.sa.CostDistance(startlyr, cost)

                # costpath = arcpy.sa.CostPath(endlyr, cost, backlink)

                path = cost_path(geometries[k], cost, radius, startneigh[0], endneigh[0])

                arcpy.AddMessage('PATH created!')
                arcpy.AddMessage(path)

                # nppath = arcpy.RasterToNumPyArray(costpath, nodata_to_value=-1)
                # npdist = arcpy.RasterToNumPyArray(costdist, nodata_to_value=-1)
                #
                # cellidx = numpy.where(nppath >= 0)
                # cells = numpy.argwhere(nppath >= 0)
                #
                # values = npdist[cellidx]
                #
                # idx = numpy.argsort(values)
                #
                # path = list(map(tuple, cells[idx, :]))

                if startdep:
                    nl = len(path)
                    for i in range(nl):
                        if path[i] not in startstr:
                            nl = i - 1
                            path = path[nl:]
                            break
                if enddep:
                    nl = len(path)
                    for i in range(1, nl):
                        if path[i] in endstr:
                            nl = i + 1
                            break
                    if extend:
                        streams.append(stream + path[1:nl])
                        types.append('Extended Stream')
                    else:
                        streams.append(path[:nl])
                        if startdep:
                            types.append('Braid/Channel')
                        else:
                            types.append('Path')
                else:
                    streams.append(path)
                    if startdep:
                        types.append('Braid/Channel')
                    else:
                        types.append('Path')

            arcpy.SetProgressorPosition(k)

        outraster = None
        nodatavalue = 0 if (min(ordids) > 0) else min(ordids) - 1

        N = len(streams)

        arcpy.AddMessage(N)

        if (N > 0):
            outraster = numpy.full((N, ni, nj), nodatavalue)
            for l in range(N):
                for ncells in range(len(streams[l])):
                    outraster[l, streams[l][ncells][0], streams[l][ncells][1]] = ordids[l]

        arcpy.AddMessage(numpy.amax(outraster))
        arcpy.AddMessage(numpy.amin(outraster))

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

        arcpy.AddMessage(int(arcpy.GetCount_management(result).getOutput(0)))
        # result = set_values(result, 'type', types)

        arcpy.UnsplitLine_management(result, outstreams, ['grid_code', 'type'])

        # ensure right direction
        arcpy.AddMessage('ENSURING RIGHT DIRECTION...')

        # with  arcpy.da.UpdateCursor(outstreams, ["SHAPE@", 'grid_code']) as rows:
        #     for row in rows:
        #         line = row[0].getPart(0)
        #         count_start = [line[0].X, line[0].Y]
        #         id = row[1]
        #
        #         idx = numpy.where(ordids == id)[0].tolist()[0]
        #
        #         hydro_start = startxy[idx]
        #         hydro_end = endxy[idx]
        #
        #         if euc_distance(count_start, hydro_start) > euc_distance(count_start, hydro_end):
        #             row[0] = FlipLine(line)
        #             rows.updateRow(row)

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
    ids = get_values(instreams_crop, inIDfield)

    # Get start and endpoints of rivers
    startpts = 'in_memory/startpts'
    endpts = 'in_memory/endpts'

    arcpy.FeatureVerticesToPoints_management(instreams_crop, startpts, point_location = 'START')
    arcpy.FeatureVerticesToPoints_management(instreams_crop, endpts, point_location = 'END')

    # Calculate distance to determine hierarchy
    end_tbl = 'in_memory/end_tbl'
    arcpy.GenerateNearTable_analysis(endpts, instreams_crop, end_tbl, closest=False, closest_count=2)
    ins = get_values(end_tbl, 'IN_FID')
    nears_end = get_values(end_tbl, 'NEAR_FID')
    dist_end = get_values(end_tbl, 'NEAR_DIST')

    start_tbl = 'in_memory/start_tbl'
    arcpy.GenerateNearTable_analysis(startpts, instreams_crop, start_tbl, closest=False, closest_count=2)
    nears_start = get_values(start_tbl, 'NEAR_FID')
    dist_start = get_values(start_tbl, 'NEAR_DIST')

    # dependent are streams which endpoints are located
    # less or equal to cellsize from another

    ordids = numpy.empty(0, dtype=numpy.int16)
    ordends = numpy.empty(0, dtype=numpy.int16)
    ordstarts = numpy.empty(0, dtype=numpy.int16)

    if len(ids) == 1:
        ordids = ids
        ordends = numpy.asarray([-1])
        ordstarts = numpy.asarray([-1])
    else:

        # remove self-distances
        flt = ins != nears_end
        nears_end = nears_end[flt]
        dist_end = dist_end[flt]

        flt = ins != nears_start
        nears_start = nears_start[flt]
        dist_start = dist_start[flt]

        depends = ids[nears_end - 1]
        depstarts = ids[nears_start - 1]

        depends[dist_end > cellsize] = -1
        depstarts[dist_start > cellsize] = -1

        flt = numpy.logical_and(depends == -1, depstarts == -1)

        arcpy.AddMessage('INDEPENDENT STREAMS: ' + str(ids[flt]))
        arcpy.AddMessage('BRAIDED STREAMS/CHANNELS (parent, braid): ' + str(zip(depstarts[depstarts != -1], ids[depstarts != -1])))

        depids = ids.copy()

        while(len(depends) > 0):
            ord = numpy.logical_and(numpy.logical_not(numpy.in1d(depends, depids)),
                                    numpy.logical_not(numpy.in1d(depstarts, depids)))

            ordids = numpy.append(ordids, depids[ord])
            ordends = numpy.append(ordends, depends[ord])
            ordstarts = numpy.append(ordstarts, depstarts[ord])

            not_ord = numpy.logical_not(ord)

            depids = depids[not_ord]
            depends = depends[not_ord]
            depstarts = depstarts[not_ord]

    process_raster(instreams_crop, inIDfield, inraster, minacc, radius, deviation, demRaster, penalty,
                   startpts, endpts, ids, ordids, ordends, ordstarts, lowerleft, cellsize, crs, outstreams)

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