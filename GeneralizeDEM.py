# -*- coding: cp1251 -*-
# DEM generalization algorithm
# Important note:  Widen and TraceFlowLines scripts must be in the same directory,
# otherwise Generalize DEM tool will not be able reference them
# 2015-2017 Timofey Samsonov, Lomonosov Moscow State University

import arcpy
import math
import pp
from arcpy.sa import *
import os.path, ExtractStreams, WidenLandforms
__author__ = 'Timofey Samsonov'


def call(oid,
         demdataset,
         marine,
         fishbuffer,
         minacc1,
         minlen1,
         minacc2,
         minlen2,
         is_widen,
         widentype,
         widendist,
         filtersize,
         is_smooth,
         workspace,
         rastertinworkspace):

    i = int(oid)-1
    arcpy.AddMessage(i)
    raster = 'dem0_' + str(i)

    arcpy.AddMessage('Counting')
    N = int(arcpy.GetCount_management(fishbuffer).getOutput(0))

    arcpy.AddMessage('Selecting cell')
    cells = arcpy.da.SearchCursor(fishbuffer, ['SHAPE@', 'OID@'])
    cell = cells.next()
    while cell[1] != oid:
        cell = cells.next()

    arcpy.AddMessage('---')
    arcpy.AddMessage('GENERALIZING DEM ' + str(i + 1) + ' FROM ' + str(N))

    dem0 = arcpy.Raster(rastertinworkspace + '/' + raster)
    dem = dem0

    marine_area = None
    process_marine = False
    if marine:
        marine_area = workspace + "/land" + str(i)
        arcpy.Clip_analysis(marine, cell[0], marine_area)
        if int(arcpy.GetCount_management(marine_area).getOutput(0)) > 0:
            cell_erased = workspace + "/cell_erased" + str(i)
            arcpy.Erase_analysis(cell[0], marine_area, cell_erased)
            dem = ExtractByMask(dem0, cell_erased)
            dem.save(rastertinworkspace + '/' + raster + "_e")
            dem = arcpy.Raster(rastertinworkspace + '/' + raster + "_e")
            process_marine = True

    cellsize = dem.meanCellHeight

    arcpy.AddMessage("PREPROCESSING")

    arcpy.AddMessage("Fill...")
    fill = Fill(dem, "")
    arcpy.AddMessage("Dir...")
    dir = FlowDirection(fill, "", "")
    arcpy.AddMessage("Acc...")
    acc = FlowAccumulation(dir, "", "INTEGER")

    # MAIN STREAMS AND WATERSHEDS
    arcpy.AddMessage("PROCESSING PRIMARY STREAMS AND WATERSHEDS")

    str1_0 = workspace + "/str10"

    arcpy.AddMessage("Extracting primary streams...")

    ExtractStreams.execute(acc, str1_0, minacc1, minlen1)
    str1 = SetNull(str1_0, 1, "value = 0")

    arcpy.AddMessage("Vectorizing streams...")
    streams1 = workspace + "/streams1"
    StreamToFeature(str1, dir, streams1, True)

    arcpy.AddMessage("Deriving endpoints...")
    endpoints1 = workspace + "/endpoints1"
    arcpy.FeatureVerticesToPoints_management(streams1, endpoints1, "END")

    endbuffers1 = workspace + "/endbuffers1"
    radius = 2 * cellsize
    arcpy.AddMessage("Buffering endpoints...")
    arcpy.Buffer_analysis(endpoints1, endbuffers1, radius, "FULL", "ROUND", "NONE", "")

    rendbuffers1 = workspace + "/rendbuffers1"
    # arcpy.AddMessage(rendbuffers1)
    # rendbuffers1 = "X:/Work/Scripts & Tools/MY/DEMGEN/Scratch.gdb/rendbuffers1"
    arcpy.FeatureToRaster_conversion(endbuffers1, "OBJECTID", rendbuffers1, cellsize)

    mask = CreateConstantRaster(0, "INTEGER", cellsize, dem.extent)
    # mask.save(workspace + "/mask")

    # rendbuffers1 = arcpy.Raster(rendbuffers1) + 1

    arcpy.Mosaic_management(mask, rendbuffers1, "MAXIMUM", "FIRST", "", "", "", "0.3", "NONE")

    # str1.save(workspace + "/str1")
    # str1 = workspace + "/str1"

    arcpy.AddMessage("Erasing streams...")

    str1_e = SetNull(rendbuffers1, str1, "value > 0")

    arcpy.AddMessage("Vectorizing erased streams...")
    streams1_e = workspace + "/streams1_e"
    StreamToFeature(str1_e, dir, streams1_e, True)

    arcpy.AddMessage("Deriving erased endpoints...")
    endpoints1_e = workspace + "/endpoints1_e"
    arcpy.FeatureVerticesToPoints_management(streams1_e, endpoints1_e, "END")

    arcpy.AddMessage("Deriving primary watersheds...")
    pour1 = SnapPourPoint(endpoints1_e, acc, cellsize * 1.5, "")

    wsh1 = Watershed(dir, pour1, "")

    arcpy.AddMessage("Vectorizing primary watersheds...")
    watersheds1 = workspace + "/watersheds1"
    arcpy.RasterToPolygon_conversion(wsh1, watersheds1, True, "")

    # SECONDARY STREAMS AND WATERSHEDS

    arcpy.AddMessage("PROCESSING SECONDARY STREAMS AND WATERSHEDS")

    arcpy.AddMessage("Extracting secondary streams...")
    str2_0 = workspace + "/str2"
    ExtractStreams.execute(acc, str2_0, minacc2, minlen2)

    str2 = SetNull(str2_0, 1, "value = 0")
    str2_e = SetNull(str1_0, str2, "value > 0")
    acc_e = SetNull(str1_0, acc, "value > 0")

    arcpy.AddMessage("Vectorizing streams...")
    streams2_e = workspace + "/streams2_e"
    StreamToFeature(str2_e, dir, streams2_e, True)

    arcpy.AddMessage("Deriving endpoints...")
    endpoints2_e = workspace + "/endpoints2_e"
    arcpy.FeatureVerticesToPoints_management(streams2_e, endpoints2_e, "END")

    arcpy.AddMessage("Buffering primary streams...")
    streambuffer = workspace + "/streams1_b"
    arcpy.Buffer_analysis(streams1, streambuffer, radius, "FULL", "ROUND", "NONE", "")

    arcpy.AddMessage("Selecting endpoints...")
    pointslyr = "points"
    arcpy.MakeFeatureLayer_management(endpoints2_e, pointslyr)
    arcpy.SelectLayerByLocation_management(pointslyr, "INTERSECT", streambuffer)

    pourpts2 = arcpy.CreateFeatureclass_management("in_memory", "pourpts2", "POINT", pointslyr, "DISABLED",
                                                   "DISABLED",
                                                   arcpy.Describe(endpoints2_e).spatialReference)
    arcpy.CopyFeatures_management(pointslyr, pourpts2)

    arcpy.AddMessage("Deriving secondary pour pts 1...")
    pour21 = SnapPourPoint(pourpts2, acc_e, cellsize * 1.5, "")

    arcpy.AddMessage("Deriving secondary pour pts 2...")
    pour22 = SnapPourPoint(pourpts2, acc_e, cellsize * 2, "")

    arcpy.AddMessage("Mosaic secondary pour pts...")
    arcpy.Mosaic_management(pour21, pour22, "FIRST", "FIRST", "0", "0", "", "0.3", "NONE")

    arcpy.AddMessage("Deriving secondary watersheds...")
    wsh2 = Watershed(dir, pour22, "")

    arcpy.AddMessage("Vectorizing secondary watersheds...")
    watersheds2 = workspace + "/watersheds2"
    arcpy.RasterToPolygon_conversion(wsh2, watersheds2, True, "")

    arcpy.AddMessage("Interpolating features into 3D...")

    streams1_3d = workspace + "/streams1_3d"

    arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, streams1, streams1_3d)

    watersheds1_3d = workspace + "/watersheds1_3d"
    arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, watersheds1, watersheds1_3d)

    watersheds2_3d = workspace + "/watersheds2_3d"
    arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, watersheds2, watersheds2_3d)

    marine_3d = workspace + "/marine_3d"
    if process_marine:
        arcpy.InterpolateShape_3d(demdataset, marine_area, marine_3d)

    # GENERALIZED TIN SURFACE

    arcpy.AddMessage("DERIVING GENERALIZED SURFACE")

    arcpy.AddMessage("TIN construction...")

    tin = rastertinworkspace + "/tin"
    features = []
    s1 = "'" + streams1_3d + "' Shape.Z " + "hardline"
    w1 = "'" + watersheds1_3d + "' Shape.Z " + "softline"
    w2 = "'" + watersheds2_3d + "' Shape.Z " + "softline"

    features.append(s1)
    features.append(w1)
    features.append(w2)
    if process_marine:
        m2 = "'" + marine_3d + "' Shape.Z " + "hardline"
        features.append(m2)

    featurestring = ';'.join(features)

    arcpy.ddd.CreateTin(tin, "", featurestring, "")

    # GENERALIZED RASTER SURFACE
    arcpy.AddMessage("TIN to raster conversion...")

    rastertin = rastertinworkspace + "/rastertin"
    try:
        arcpy.TinRaster_3d(tin, rastertin, "FLOAT", "NATURAL_NEIGHBORS", "CELLSIZE " + str(cellsize), 1)
    except Exception:
        arcpy.AddMessage("Failed to rasterize TIN using NATURAL_NEIGHBORS method. Switching to linear")
        arcpy.TinRaster_3d(tin, rastertin, "FLOAT", "LINEAR", "CELLSIZE " + str(cellsize), 1)

    # POSTPROCESSING
    arcpy.AddMessage("POSTPROCESSING")

    # Widen valleys and ridges
    widenraster = rastertin
    if is_widen == "true":
        arcpy.AddMessage("Raster widening...")
        widenraster = workspace + "/widenraster"
        WidenLandforms.execute(rastertin, streams1, widendist, filtersize, widenraster, widentype)

    # Smooth DEM
    result = arcpy.Raster(widenraster)
    if is_smooth == "true":
        arcpy.AddMessage("Raster filtering...")
        neighborhood = NbrRectangle(filtersize, filtersize, "CELL")
        result = FocalStatistics(widenraster, neighborhood, "MEAN", "DATA")

    if process_marine:
        arcpy.AddMessage("Masking marine regions...")
        result_erased = ExtractByMask(result, cell_erased)
        arcpy.Mosaic_management(result_erased, rastertin, "FIRST", "FIRST", "", "", "", "0.3", "NONE")
        arcpy.AddMessage("Saving result...")
        res = arcpy.Raster(rastertin)
        res.save(rastertinworkspace + '/gen/dem' + str(i))
    else:
        arcpy.AddMessage("Saving result...")
        result.save(rastertinworkspace + '/gen/dem' + str(i))

    arcpy.Delete_management(pointslyr)

def worker(oid):
    return oid

def execute(demdataset,
            marine,
            output,
            outputcellsize,
            minacc1,
            minlen1,
            minacc2,
            minlen2,
            is_widen,
            widentype,
            widendist,
            filtersize,
            is_smooth,
            is_parallel,
            tilesize):

    arcpy.CheckOutExtension("3D")
    arcpy.CheckOutExtension("Spatial")

    arcpy.env.overwriteOutput = True

    workspace = os.path.dirname(output)

    # raster workspace MUST be a folder, no
    rastertinworkspace = workspace
    n = len(rastertinworkspace)
    if n > 4:
        end = rastertinworkspace[n - 4: n]  # extract last 4 letters
        if end == ".gdb":  # geodatabase
            rastertinworkspace = os.path.dirname(rastertinworkspace)

    arcpy.CreateFolder_management(rastertinworkspace, 'scratch')
    rastertinworkspace += '/scratch'

    arcpy.AddMessage(rastertinworkspace)

    arcpy.CreateFolder_management(rastertinworkspace, 'gen')
    arcpy.CreateFolder_management(rastertinworkspace, 'gencrop')

    arcpy.env.scratchWorkspace = rastertinworkspace
    arcpy.env.workspace = arcpy.env.scratchWorkspace

    demsource = arcpy.Raster(demdataset)

    nrows = int(math.ceil(float(demsource.height) / float(tilesize)))
    ncols = int(math.ceil(float(demsource.width) / float(tilesize)))
    total = nrows * ncols

    cellsize = 0.5 * (demsource.meanCellHeight + demsource.meanCellWidth)

    arcpy.AddMessage(demsource.height)
    arcpy.AddMessage(demsource.width)
    arcpy.AddMessage(cellsize)
    arcpy.AddMessage(nrows)
    arcpy.AddMessage(ncols)

    bufferpixelwidth = math.ceil(max(demsource.width, demsource.height) / (max(nrows, ncols) * 10))
    bufferwidth = bufferpixelwidth * cellsize

    arcpy.AddMessage('Splitting raster into ' + str(nrows) + ' x ' + str(ncols) + ' = ' + str(total) + ' tiles')
    arcpy.AddMessage('Tile overlap will be ' + str(2 * bufferpixelwidth) + ' pixels')

    arcpy.AddMessage('Making fishnet ')
    fishnet = workspace + "/fishnet"
    arcpy.AddMessage('fishnet')
    arcpy.CreateFishnet_management(fishnet,
                                   str(demsource.extent.XMin) + ' ' + str(demsource.extent.YMin),
                                   str(demsource.extent.XMin) + ' ' + str(demsource.extent.YMin + 1),
                                   '', '',
                                   nrows, ncols,
                                   '', '',
                                   demsource.extent, 'POLYGON')
    arcpy.AddMessage('Making buffer ')
    fishbuffer = workspace + "/fishbuffer"
    arcpy.Buffer_analysis(fishnet, fishbuffer, bufferwidth)

    arcpy.AddMessage('Making mask buffer ')
    fishmaskbuffer = workspace + "/fishmaskbuffer"
    arcpy.Buffer_analysis(fishnet, fishmaskbuffer, max(demsource.meanCellHeight, demsource.meanCellWidth))

    arcpy.AddMessage('Splitting raster ')
    arcpy.SplitRaster_management(demdataset,
                                 rastertinworkspace,
                                 'dem',
                                 'POLYGON_FEATURES',
                                 'GRID',
                                 '', '', '', '', '', '', '',
                                 fishbuffer)

    arcpy.AddMessage('Creating cursor')
    rows = arcpy.da.SearchCursor(fishbuffer, 'OID@')
    oids = [row[0] for row in rows]

    # MAIN PROCESSING
    if is_parallel == 'true':
        arcpy.AddMessage('Trying to make multiprocessing')

        prc = int(os.environ["NUMBER_OF_PROCESSORS"]) - 1

        arcpy.AddMessage('Detected ' + str(prc) + ' processors')

        # using Parellel python
        ppservers = ()

        # sets workers to prc
        job_server = pp.Server(prc, ppservers=ppservers)

        jobs = []

        for oid in oids:
            jobs.append(job_server.submit(func=worker, args=(oid,)))

        # Retrieve results of all submited jobs
        for job in jobs:
            print job()
            arcpy.AddMessage(job())

    else:
        for oid in oids:
            try:
                call(oid,
                     demdataset,
                     marine,
                     fishbuffer,
                     minacc1,
                     minlen1,
                     minacc2,
                     minlen2,
                     is_widen,
                     widentype,
                     widendist,
                     filtersize,
                     is_smooth,
                     workspace,
                     rastertinworkspace)
            except Exception:
                arcpy.AddMessage("Failed to generalize dem" + str(oid - 1))

    arcpy.AddMessage("CLIPPING AND MASKING GENERALIZED RASTERS")

    rows = arcpy.da.SearchCursor(fishmaskbuffer, ['SHAPE@', 'OID@'])
    i = 0
    for row in rows:
        dem = arcpy.Raster(rastertinworkspace + '/gen/dem' + str(i))
        dem_clipped = ExtractByMask(dem, row[0])
        dem_clipped.save(rastertinworkspace + '/gencrop/dem' + str(i))
        i += 1

    arcpy.env.workspace = rastertinworkspace + '/gencrop/'

    rasters = arcpy.ListRasters("*", "GRID")

    rasters_str = ';'.join(rasters)

    arcpy.MosaicToNewRaster_management(rasters_str,
                                       os.path.dirname(output),
                                       os.path.basename(output),
                                       "",
                                       "16_BIT_SIGNED",
                                       str(outputcellsize),
                                       "1",
                                       "BLEND",
                                       "FIRST")

    return

if __name__ == '__main__':

    demdataset = arcpy.GetParameterAsText(0)
    marine = arcpy.GetParameterAsText(1)
    output = arcpy.GetParameterAsText(2)
    outputcellsize = float(arcpy.GetParameterAsText(3))
    minacc1 = int(arcpy.GetParameterAsText(4))
    minlen1 = int(arcpy.GetParameterAsText(5))
    minacc2 = int(arcpy.GetParameterAsText(6))
    minlen2 = int(arcpy.GetParameterAsText(7))
    is_widen = arcpy.GetParameterAsText(8)
    widentype = arcpy.GetParameterAsText(9)
    widendist = float(arcpy.GetParameterAsText(10))
    filtersize = int(arcpy.GetParameterAsText(11))
    is_smooth = arcpy.GetParameterAsText(12)
    is_parallel = arcpy.GetParameterAsText(13)
    tilesize = int(arcpy.GetParameterAsText(14))

    execute(demdataset, marine, output, outputcellsize,
            minacc1, minlen1, minacc2, minlen2,
            is_widen, widentype, widendist, filtersize,
            is_smooth, is_parallel, tilesize)







