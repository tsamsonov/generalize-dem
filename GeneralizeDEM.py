# -*- coding: cp1251 -*-
# DEM generalization algorithm
# Important note:  ExtractStreams and WidenLandforms scripts must be in the same directory,
# otherwise Generalize DEM tool will not be able reference them
# 2015-2017 Timofey Samsonov, Lomonosov Moscow State University

import arcpy
import math
import pp
import multiprocessing
import time
import shutil
import sys
import traceback
from arcpy.sa import *
import os.path, ExtractStreams, WidenLandforms
__author__ = 'Timofey Samsonov'

def tin_to_raster(tin, rastertin, rastertype, method, cellsize, factor):
    arcpy.TinRaster_3d(tin, rastertin, rastertype, method, cellsize, factor)
    return

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
         scratchworkspace):
    try:
        arcpy.CheckOutExtension("Spatial")
        arcpy.CheckOutExtension("3D")

        arcpy.CreateFolder_management(scratchworkspace + '/processing', 'dem' + str(oid))

        rastertinworkspace = scratchworkspace + '/processing/dem' + str(oid)

        gdb_name = 'dem' + str(oid) + '.gdb'

        arcpy.CreateFileGDB_management(rastertinworkspace, gdb_name)

        workspace = rastertinworkspace + '/' + gdb_name

        arcpy.env.scratchWorkspace = rastertinworkspace
        arcpy.env.workspace = workspace

        i = int(oid)-1
        raster = 'dem0_' + str(i)

        # arcpy.AddMessage('Counting')
        N = int(arcpy.GetCount_management(fishbuffer).getOutput(0))

        # arcpy.AddMessage('Selecting cell')
        cells = arcpy.da.SearchCursor(fishbuffer, ['SHAPE@', 'OID@'])
        cell = cells.next()
        while cell[1] != oid:
            cell = cells.next()

        arcpy.AddMessage('---')
        arcpy.AddMessage('GENERALIZING DEM ' + str(i + 1) + ' FROM ' + str(N))

        dem0 = arcpy.Raster(scratchworkspace + '/source/' + raster)
        dem = dem0

        marine_area = None
        process_marine = False
        if marine:
            arcpy.AddMessage("Extracting marine area...")
            marine_area = workspace + "/land" + str(i)
            arcpy.Clip_analysis(marine, cell[0], marine_area)
            if int(arcpy.GetCount_management(marine_area).getOutput(0)) > 0:
                cell_erased = workspace + "/cell_erased" + str(i)
                arcpy.Erase_analysis(cell[0], marine_area, cell_erased)

                nareas = int(arcpy.GetCount_management(cell_erased).getOutput(0))

                if nareas == 0:
                    arcpy.AddMessage('NOTHING TO GENERALIZE: The cell is completely in the marine area. Finishing...')
                    return False
                else:
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
        arcpy.FeatureToRaster_conversion(endbuffers1, "OBJECTID", rendbuffers1, cellsize)

        mask = CreateConstantRaster(0, "INTEGER", cellsize, dem.extent)

        arcpy.Mosaic_management(mask, rendbuffers1, "MAXIMUM", "FIRST", "", "", "", "0.3", "NONE")

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

        # pourpts2 = arcpy.CreateFeatureclass_management(workspace, "pourpts2", "POINT", pointslyr, "DISABLED",
        #                                                "DISABLED",
        #                                                arcpy.Describe(endpoints2_e).spatialReference)

        pourpts2 = workspace + "/pourpts2"
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
            # arcpy.AddMessage("Creating process")
            # p = multiprocessing.Process(tin_to_raster,
            #                             args=(tin, rastertin, "FLOAT", "NATURAL_NEIGHBORS", "CELLSIZE " + str(cellsize), 1)
            # )
            #
            # arcpy.AddMessage("Starting process")
            # p.start()
            #
            # # Wait for 30 seconds or until process finishes
            # arcpy.AddMessage("Waiting for finish")
            # p.join(30)
            #
            # # If thread is still active
            # if p.is_alive():
            #     print "Rasterizing TIN using NATURAL_NEIGHBORS takes too long time (> 1 min). It seems that process is freezed"
            #     # Terminate
            #     p.terminate()
            #     p.join()
            #     raise Exception

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
            res.save(scratchworkspace + '/gen/dem' + str(i))
        else:
            arcpy.AddMessage("Saving result...")
            result.save(scratchworkspace + '/gen/dem' + str(i))

        arcpy.AddMessage("CLEANING TEMPORARY DATA")

        arcpy.Delete_management(str1_0)
        arcpy.Delete_management(str2_0)
        arcpy.Delete_management(streams1)
        arcpy.Delete_management(endpoints1)
        arcpy.Delete_management(endbuffers1)
        arcpy.Delete_management(rendbuffers1)
        arcpy.Delete_management(streams1_e)
        arcpy.Delete_management(endpoints1_e)
        arcpy.Delete_management(watersheds1)
        arcpy.Delete_management(streams2_e)
        arcpy.Delete_management(endpoints2_e)
        arcpy.Delete_management(streambuffer)
        arcpy.Delete_management(pointslyr)
        arcpy.Delete_management(pourpts2)
        arcpy.Delete_management(watersheds2)
        arcpy.Delete_management(streams1_3d)
        arcpy.Delete_management(watersheds1_3d)
        arcpy.Delete_management(watersheds2_3d)
        arcpy.Delete_management(tin)
        arcpy.Delete_management(rastertin)
        arcpy.Delete_management(widenraster)

        arcpy.Delete_management(fill)
        arcpy.Delete_management(dir)
        arcpy.Delete_management(acc)
        arcpy.Delete_management(str1)
        arcpy.Delete_management(str2)
        arcpy.Delete_management(mask)
        arcpy.Delete_management(str1_e)
        arcpy.Delete_management(str2_e)
        arcpy.Delete_management(pour1)
        arcpy.Delete_management(pour21)
        arcpy.Delete_management(pour22)
        arcpy.Delete_management(wsh1)
        arcpy.Delete_management(wsh2)
        arcpy.Delete_management(acc_e)

        arcpy.Delete_management(workspace)
        arcpy.Delete_management(rastertinworkspace)

        arcpy.CheckInExtension("3D")
        arcpy.CheckInExtension("Spatial")

    except:
        arcpy.AddMessage("Failed to generalize tile dem" + str(oid - 1))
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)

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

    # ORGANIZE WORKSPACE

    workspace = os.path.dirname(output)
    scratchworkspace = workspace
    n = len(scratchworkspace)
    if n > 4:
        end = scratchworkspace[n - 4: n]  # extract last 4 letters
        if end == ".gdb":  # geodatabase
            scratchworkspace = os.path.dirname(scratchworkspace)

    arcpy.env.workspace = scratchworkspace
    workspaces = arcpy.ListWorkspaces("*", "Folder")
    names = [os.path.basename(w) for w in workspaces]
    i = 0
    name = 'scratch'
    while name in names:
        name = 'scratch' + str(i)
        i += 1

    arcpy.CreateFolder_management(scratchworkspace, name)
    scratchworkspace += '/' + name
    arcpy.env.workspace = scratchworkspace

    arcpy.CreateFolder_management(scratchworkspace, 'gen') # generalized tiles
    arcpy.CreateFolder_management(scratchworkspace, 'gencrop') # generalized and cropped tiles
    arcpy.CreateFolder_management(scratchworkspace, 'processing') # main processing workspace
    arcpy.CreateFolder_management(scratchworkspace, 'source') # source tiles
    arcpy.CreateFileGDB_management(scratchworkspace, 'Scratch.gdb')
    workspace = scratchworkspace + "/Scratch.gdb"


    # SPLIT THE TASK INTO TILES

    demsource = arcpy.Raster(demdataset)

    nrows = int(math.ceil(float(demsource.height) / float(tilesize)))
    ncols = int(math.ceil(float(demsource.width) / float(tilesize)))
    total = nrows * ncols

    cellsize = 0.5 * (demsource.meanCellHeight + demsource.meanCellWidth)

    bufferpixelwidth = math.ceil(max(demsource.width, demsource.height) / (max(nrows, ncols) * 10))
    bufferwidth = bufferpixelwidth * cellsize

    arcpy.AddMessage('Splitting raster into ' + str(nrows) + ' x ' + str(ncols) + ' = ' + str(total) + ' tiles')
    arcpy.AddMessage('Tile overlap will be ' + str(2 * bufferpixelwidth) + ' pixels')

    arcpy.AddMessage('Creating fishnet...')
    fishnet = workspace + "/fishnet"
    arcpy.CreateFishnet_management(fishnet,
                                   str(demsource.extent.XMin) + ' ' + str(demsource.extent.YMin),
                                   str(demsource.extent.XMin) + ' ' + str(demsource.extent.YMin + 1),
                                   '', '',
                                   nrows, ncols,
                                   '', '',
                                   demsource.extent, 'POLYGON')

    arcpy.AddMessage('Creating split buffer...')
    fishbuffer = workspace + "/fishbuffer"
    arcpy.Buffer_analysis(fishnet, fishbuffer, bufferwidth)

    arcpy.AddMessage('Creating mask buffer...')
    fishmaskbuffer = workspace + "/fishmaskbuffer"
    arcpy.Buffer_analysis(fishnet, fishmaskbuffer, max(demsource.meanCellHeight, demsource.meanCellWidth))

    arcpy.AddMessage('Splitting raster...')
    arcpy.SplitRaster_management(demdataset,
                                 scratchworkspace + "/source",
                                 'dem',
                                 'POLYGON_FEATURES',
                                 'GRID',
                                 '', '', '', '', '', '', '',
                                 fishbuffer)

    rows = arcpy.da.SearchCursor(fishbuffer, 'OID@')
    oids = [row[0] for row in rows]

    arcpy.AddMessage('')

    # PERFORM PROCESSING

    if is_parallel == 'true':
        nproc = multiprocessing.cpu_count()

        arcpy.AddMessage('> Trying to make multiprocessing using ' + str(nproc) + ' processor cores')
        arcpy.AddMessage('')

        pool = multiprocessing.Pool(nproc)
        for oid in oids:
            pool.apply_async(call, (oid,
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
                                    scratchworkspace,))
        pool.close()
        pool.join()

    else:
        arcpy.AddMessage('> Processing in sequential mode')
        arcpy.AddMessage('')
        for oid in oids:
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
                 scratchworkspace)

    # FINALIZE

    arcpy.AddMessage("CLIPPING AND MASKING GENERALIZED RASTERS")

    rows = arcpy.da.SearchCursor(fishmaskbuffer, ['SHAPE@', 'OID@'])
    i = 0
    for row in rows:
        raster = scratchworkspace + '/gen/dem' + str(i)
        if arcpy.Exists(raster):
            dem = arcpy.Raster(raster)
            dem_clipped = ExtractByMask(dem, row[0])
            dem_clipped.save(scratchworkspace + '/gencrop/dem' + str(i))
        i += 1

    arcpy.env.workspace = scratchworkspace + '/gencrop/'
    rasters = arcpy.ListRasters("*", "GRID")

    if len(rasters) > 0:
        rasters_str = ';'.join(rasters)
        arcpy.AddMessage('SAVING RESULT')
        arcpy.MosaicToNewRaster_management(rasters_str,
                                           os.path.dirname(output),
                                           os.path.basename(output),
                                           "",
                                           "16_BIT_SIGNED",
                                           str(outputcellsize),
                                           "1",
                                           "BLEND",
                                           "FIRST")
    else:
        arcpy.AddMessage('NOTHING GENERALIZED')

    return

if __name__ == '__main__':
    try:
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
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)