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

def get_window(npcost, ij, ni, nj, size=3):
    w = int((size - 1) / 2)
    l = range(-w, w + 1)  # calculate kernel indices
    idx = numpy.meshgrid(l, l)  # generate coordinate matrices

    x = idx[0] + ij[0]
    y = idx[1] + ij[1]

    flt_xy = (x >= 0) * (x < ni) * (y >= 0) * (y < nj) #* numpy.logical_or(x != ij[0], y != ij[1])

    x = x[flt_xy]
    y = y[flt_xy]

    costs = numpy.array([npcost[i, j] for i, j in zip(x, y)])
    npdist = costs * (((x - ij[0]) ** 2 + (y - ij[1]) ** 2) ** 0.5)

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

def minimax(arr):
    brr = list(arr)
    m = min(brr)
    mmax = m
    brr.sort()
    for i in range(len(brr)):
        if (m + i) in brr:
            mmax = m + i
        else:
            break
    return mmax

def maximin(arr):
    m = max(arr)
    mmin = m
    for i in range(len(arr)):
        if (m - i) in arr:
            mmin = m - i
        else:
            break
    return mmin

def longgap(arr):
    brr = list(arr)
    m = min(brr)
    gap = 0
    brr.sort()
    for i in range(len(brr)):
        if (m + i) not in brr:
            gap += 1
        else:
            gap = 0

        if gap > 1:
            return True
    return False

def cost_distance(source, npcost, npcomp, destination, nodatavalue=-1):

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

    # visited = []
    calculated = [source]

    ni = npcost.shape[0]
    nj = npcost.shape[1]

    total = (npcost != nodatavalue).sum()

    arcpy.AddMessage('Need to visit ' + str(total) + ' cells')

    npdist = numpy.full((ni, nj), float('Inf'))
    npback = numpy.full((ni, nj), nodatavalue)
    npminimax = numpy.full((ni, nj), nodatavalue).astype(int)
    # npmaximin = numpy.full((ni, nj), nodatavalue).astype(int)

    npdist[source] = 0
    npback[source] = 0
    npminimax[source] = minimax(npcomp[source[0]][source[1]])
    # npmaximin[source] = maximin(npcomp[source[0]][source[1]])

    dead_cells = []
    dead_dist = []

    while len(calculated) > 0:
        newcalc = []
        newdist = []
        incomp_moves = 0
        moves = 0
        for cell in calculated:
            celldist = npdist[cell]
            nb, dist = get_window(npcost, cell, ni, nj)
            cell_minimax = npminimax[cell]
            cell_compat = npcomp[cell[0]][cell[1]]
            deads = 0

            for (ij, d) in zip(nb, dist):
                # if (ij not in visited) and (npcost[ij] != nodatavalue):
                if (npcost[ij] != nodatavalue):

                    ijcomp = npcomp[ij[0]][ij[1]]

                    if len(ijcomp) < 1:
                        if d < 2:
                            deads+=1
                        incomp_moves += 1
                        continue

                    if len(ijcomp.intersection(cell_compat)) == 0:
                        if d < 2:
                            deads += 1
                        incomp_moves += 1
                        continue

                    # if (cell_minimax not in ijcomp) and ((max(cell_compat) != minimax(cell_compat)) or (minimax(ijcomp) != max(ijcomp))):
                    # ijminimax = max(ijcomp)
                    # if cell_minimax in ijcomp:
                    #     ijminimax = minimax(set(range(cell_minimax, max(ijcomp) + 1)).intersection(ijcomp))

                    if min(ijcomp) > cell_minimax:
                        if d < 2:
                            deads += 1
                        incomp_moves += 1
                        continue

                    accum_dist = celldist + d
                    if (accum_dist < npdist[ij]): #or ((npdist[ij] < float('inf')) and (npminimax[ij] < maximin(cell_compat))):
                        npdist[ij] = accum_dist
                        npback[ij] = invback(cell, ij)

                        if cell_minimax in ijcomp:
                            npminimax[ij] = minimax(set(range(cell_minimax, max(ijcomp) + 1)).intersection(ijcomp))
                        else:
                            npminimax[ij] = minimax(ijcomp)
                        if (ij not in newcalc):
                            newcalc.append(ij)
                            newdist.append(accum_dist)
                        else:
                            newdist[newcalc.index(ij)] = accum_dist
                        moves += 1
                    # elif npdist[ij] < float('inf'):
                    #     type = npback[cell]
                    #     while type != 0:
                    #         shift = back[npback[ij]]
                    #         ij = (ij[0] + shift[0], ij[1] + shift[1])
                    #         type = npback[ij]
                    #         if ()


        calculated = [ij for dist, ij in sorted(zip(newdist, newcalc))]

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
    npcost = arcpy.RasterToNumPyArray(cost, nodata_to_value = -1)


    ni = npcost.shape[0]
    nj = npcost.shape[1]

    # Generate compatibility raster
    arcpy.AddMessage('Generating compatibility sets')
    npcomp = [[set() for j in range(nj)] for i in range(ni)] # compatibility raster

    k = 0
    for pnt in coords:
        ip = ni - math.trunc((pnt[1] - miny) / cellsize) - 1
        jp = math.trunc((pnt[0] - minx) / cellsize)
        nb = get_neighborhood(ip, jp, radius, cellsize, ni, nj)
        for (i, j) in nb:
            npcomp[i][j] = npcomp[i][j].union([int(k)])
        k +=1

    arcpy.AddMessage('Calculating distance and backlink rasters')
    npdist, npback = cost_distance(source, npcost, npcomp, destination)

    ij = destination
    path = [ij]
    type = npback[ij]
    k = 0

    arcpy.AddMessage('Tracing')

    while type != 0:
        shift = back[npback[ij]]
        ij = (ij[0] + shift[0], ij[1] + shift[1])

        path.append(ij)
        type = npback[ij]

        k += 1

    path.reverse()

    return path, npdist, npback

def extended_euc(coords, crs, cellsize, start = None, end = None):

    templine = 'in_memory/templine'
    arcpy.CreateFeatureclass_management('in_memory', 'templine', "POLYLINE", spatial_reference=crs)

    cursor = arcpy.da.InsertCursor(templine, ["SHAPE@"])

    points = []
    if start != None:
        points.append(arcpy.Point(start[0], start[1]))
    for xy in coords:
        points.append(arcpy.Point(xy[0], xy[1]))
    if end != None:
        points.append(arcpy.Point(end[0], end[1]))

    array = arcpy.Array(points)
    line = arcpy.Polyline(array)
    cursor.insertRow([line])

    del cursor

    return arcpy.RasterToNumPyArray(arcpy.sa.EucDistance(templine, cell_size=cellsize))


def process_raster(instreams, inIDfield, in_raster, minacc, radius, deviation, demraster, penalty, startpts, endpts,
                   ids, ordids, ordends, ordstarts, lowerleft, cellsize, crs, outstreams, limit):

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
        arcpy.AddMessage(idx)
        if len(geometries) > 1:
            geometries = [geometries[i] for i in idx]

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

            if enddep or startdep:
                xystart = (minx + startneigh[0][1] * cellsize, miny + (ni - startneigh[0][0]) * cellsize) if startdep else None
                xyend = (minx + endneigh[0][1] * cellsize, miny + (ni - endneigh[0][0]) * cellsize) if enddep else None

                eucs[k, :, :] = extended_euc(geometries[k], crs, cellsize, xystart, xyend)

            arcpy.AddMessage("ID = " + str(ordids[k]) + ' (' + str(k + 1) + " from " + str(n) + ')')

            stream = []
            extend = False

            if not startdep:
                weight = float('Inf')
                for (i, j) in startneigh:
                    if  inraster[i, j] > minacc:
                        s, e = trace_flow_cells(extinraster, eucs[k,:,:], i, j, minacc, endneigh)
                        ncells = len(e)
                        if ncells > 0:
                            coords = []
                            for (i, j) in s:
                                coords.append((minx + j * cellsize, miny + (ni - i) * cellsize))
                            dev = Utils.dist_fun[limit](coords, geometries[k])
                            # arcpy.AddMessage(coords)
                            # arcpy.AddMessage(geometries[k])
                            # arcpy.AddMessage(dev)
                            if (dev <= deviation):
                                w = Utils.hausdorff_dist_mod(coords, geometries[k])
                                if (w < weight):
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

                euc_mask = (euc + 1) * arcpy.sa.Reclassify(euc, "value",
                                                           arcpy.sa.RemapRange([[0, deviation, penalty],
                                                                                [deviation, euc.maximum,
                                                                                 'NODATA']]))
                # euc_mask = arcpy.sa.Reclassify(euc, "value",
                #                                            arcpy.sa.RemapRange([[0, deviation, penalty],
                #                                                                 [deviation, euc.maximum,
                #                                                                  'NODATA']]))

                strs = arcpy.sa.Reclassify(in_raster, "value",
                                           arcpy.sa.RemapRange([[0, minacc, 'NODATA'], [minacc, MAXACC, 1]]))

                cost = arcpy.sa.ExtractByMask(strs, euc_mask) #* euc

                arcpy.Mosaic_management(euc_mask * arcpy.sa.Raster(demraster), cost, 'MINIMUM')
                # arcpy.Mosaic_management(euc_mask, cost, 'MINIMUM')

                # cost.save('X:/DEMGEN/cost_new.tif')

                # return

                # arcpy.env.workspace = "X:/DEMGEN/"
                # if arcpy.Exists("cost.tif"):
                #     arcpy.Delete_management("cost.tif")
                #
                # cost.save('cost.tif')

                backlink = arcpy.sa.CostBackLink(startlyr, cost)
                costdist = arcpy.sa.CostDistance(startlyr, cost)

                costpath = arcpy.sa.CostPath(endlyr, cost, backlink)

                nppath = arcpy.RasterToNumPyArray(costpath, nodata_to_value=-1)
                npdist = arcpy.RasterToNumPyArray(costdist, nodata_to_value=-1)

                cellidx = numpy.where(nppath >= 0)
                cells = numpy.argwhere(nppath >= 0)

                values = npdist[cellidx]

                idx = numpy.argsort(values)

                path = list(map(tuple, cells[idx, :]))

                # path, npdist, npback = cost_path(geometries[k], cost, radius, startneigh[0], endneigh[0])
                #
                # ras = arcpy.NumPyArrayToRaster(npdist, lowerleft, cellsize, value_to_nodata=float('Inf'))
                # arcpy.DefineProjection_management(ras, crs)
                # ras.save('X:/DEMGEN/dist.tif')
                #
                # ras = arcpy.NumPyArrayToRaster(npback, lowerleft, cellsize, value_to_nodata=-1)
                # arcpy.DefineProjection_management(ras, crs)
                # ras.save('X:/DEMGEN/back.tif')
                #
                # arcpy.AddMessage('PATH created!')
                # arcpy.AddMessage(path)

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
        arcpy.Densify_edit(outstreams, 'DISTANCE', cellsize)

        arcpy.AddField_management(outstreams, 'frechet_dist', 'FLOAT', field_length=16)
        arcpy.AddField_management(outstreams, 'hausdorff_dist', 'FLOAT', field_length=16)
        arcpy.AddField_management(outstreams, 'dir_hausdorff_dist', 'FLOAT', field_length=16)
        arcpy.AddField_management(outstreams, 'quality', 'TEXT', field_length=16)

        # ensure right direction
        arcpy.AddMessage('ENSURING RIGHT DIRECTION AND ASSESSING THE QUALITY...')

        with  arcpy.da.UpdateCursor(outstreams, ["SHAPE@", 'grid_code', 'frechet_dist', 'hausdorff_dist', 'dir_hausdorff_dist', 'quality']) as rows:
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

                coords = []
                for pnt in row[0].getPart().next():
                    coords.append([pnt.X, pnt.Y])

                row[2] = Utils.frechet_dist(coords, geometries[idx])
                row[3] = Utils.hausdorff_dist(coords, geometries[idx])
                row[4] = Utils.hausdorff_dist_dir(coords, geometries[idx])
                rows.updateRow(row)

                if row[2] <= deviation:
                    row[5] = 'Strong'
                elif row[3] <= deviation:
                    row[5] = 'Regular'
                else:
                    row[5] = 'Weak'
                rows.updateRow(row)
        return

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        raise Exception

def execute(in_streams, inIDfield, inraster, demRaster, outstreams, minacc, penalty, radius, deviation, limit):
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
                   startpts, endpts, ids, ordids, ordends, ordstarts, lowerleft, cellsize, crs, outstreams, limit)

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