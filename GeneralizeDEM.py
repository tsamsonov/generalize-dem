# -*- coding: cp1251 -*-
# DEM generalization algorithm
# Important note:  all *.py scripts must be in the same directory,
# otherwise Generalize DEM tool will not be able to reference them
# 2015-2019 Timofey Samsonov, Lomonosov Moscow State University

import arcpy
import math
import multiprocessing
import sys
import traceback
from arcpy.sa import *
from itertools import repeat
import os.path
import ExtractStreams, WidenLandforms, CreateFishnet

__author__ = 'Timofey Samsonov'

# Just a crutch for pool.map (Python 2.7)
def call_list(args):
    return call(*args)

def call(oid,
         demdataset,
         flowdir,
         flowacc,
         flines,
         fpolys,
         cliparea,
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
        i = int(oid) - 1
        raster = 'dem' + str(i) + '.tif'
        N = int(arcpy.GetCount_management(fishbuffer).getOutput(0))

        arcpy.AddMessage('\n> GENERALIZING DEM ' + str(i + 1) + ' FROM ' + str(N) + '\n')

        if not arcpy.Exists(scratchworkspace + '/source/' + raster):
            arcpy.AddMessage('\nNOTHING TO GENERALIZE: tile ' + str(i + 1) + ' is empty. Finishing...\n')
            return True

        # arcpy.CheckOutExtension("Spatial")
        # arcpy.CheckOutExtension("3D")

        dem0 = arcpy.Raster(scratchworkspace + '/source/' + raster)
        dem = dem0

        arcpy.env.overwriteOutput = True

        arcpy.env.extent = dem0.extent # Very important!
        arcpy.env.snapRaster = dem0 # Very important!

        processing_folder = scratchworkspace + '/processing'
        tile_folder = 'dem' + str(oid)

        arcpy.env.workspace = processing_folder

        if not arcpy.Exists(tile_folder):
            arcpy.CreateFolder_management(processing_folder, tile_folder)
        rastertinworkspace = processing_folder + '/' + tile_folder

        gdb_name = tile_folder + '.gdb'
        if not arcpy.Exists(tile_folder + '/' + gdb_name):
            arcpy.CreateFileGDB_management(rastertinworkspace, gdb_name)

        workspace = rastertinworkspace + '/' + gdb_name

        arcpy.env.scratchWorkspace = rastertinworkspace
        arcpy.env.workspace = workspace

        cells = arcpy.da.SearchCursor(fishbuffer, ['SHAPE@', 'OID@'])
        cell = cells.next()
        while cell[1] != oid:
            cell = cells.next()

        land_area = None
        process_marine = False
        clip_3d = workspace + "/clip_3d"
        cell_erased = workspace + "/cell_erased" + str(i)

        if cliparea:
            arcpy.AddMessage("Clipping generalization area...")
            land_area = workspace + "/land" + str(i)
            arcpy.Clip_analysis(cliparea, cell[0], land_area)
            if arcpy.Exists(land_area):
                if int(arcpy.GetCount_management(land_area).getOutput(0)) > 0:
                    arcpy.Clip_analysis(cell[0], cliparea, cell_erased)
                    dem = ExtractByMask(dem0, cell_erased)
                    dem.save(rastertinworkspace + '/' + raster + "_e.tif")
                    dem = arcpy.Raster(rastertinworkspace + '/' + raster + "_e.tif")
                    arcpy.InterpolateShape_3d(demdataset, land_area, clip_3d, vertices_only=True)

                    process_marine = True

                    # nareas = 0
                    # if arcpy.Exists(cell_erased):
                    #     nareas = int(arcpy.GetCount_management(cell_erased).getOutput(0))
                    #
                    # if nareas == 0:
                else:
                    arcpy.AddMessage('\nNOTHING TO GENERALIZE: The tile is completely outside clipping area. Finishing...\n')
                    return True
            else:
                arcpy.AddMessage('\nNOTHING TO GENERALIZE: The tile is completely outside clipping area. Finishing...\n')
                return True
                        # arcpy.AddMessage("CLEANING MAIN DATA")
                        # arcpy.Delete_management(land_area)
                        # arcpy.Delete_management(cell_erased)
                        # arcpy.Delete_management(workspace)
                        # arcpy.Delete_management(rastertinworkspace)

                        # arcpy.CheckInExtension("3D")
                        # arcpy.CheckInExtension("Spatial")


        cellsize = dem.meanCellHeight

        arcpy.AddMessage("PREPROCESSING")

        dir = flowdir
        acc = flowacc
        if (dir == None):
            arcpy.AddMessage("Fill...")
            fill = Fill(dem, "")
            arcpy.AddMessage("Dir...")
            dir = FlowDirection(fill, "", "")
            arcpy.AddMessage("Acc...")
            acc = FlowAccumulation(dir, "", "INTEGER")

        # MAIN STREAMS AND WATERSHEDS
        arcpy.AddMessage("PROCESSING PRIMARY STREAMS AND WATERSHEDS")

        str1_0 = rastertinworkspace + "/str10.tif"

        arcpy.AddMessage("Extracting primary streams...")

        ExtractStreams.execute(acc, str1_0, minacc1, minlen1)

        maxstr = int(str(arcpy.GetRasterProperties_management(str1_0, "MAXIMUM")))

        tin = rastertinworkspace + "/tin"
        streams1 = workspace + "/streams1"

        stream_processing = True

        features = [] # to be filled with features for TIN construction

        # If there are any streams extracted
        if maxstr == 1:
            str1 = SetNull(str1_0, 1, "value = 0")

            arcpy.AddMessage("Vectorizing streams...")
            StreamToFeature(str1, dir, streams1, False)

            arcpy.AddMessage("Deriving endpoints...")
            endpoints1 = workspace + "/endpoints1"
            arcpy.FeatureVerticesToPoints_management(streams1, endpoints1, "END")

            # radius = 2 * cellsize

            radius = 1 # number of cells

            rpts1 = workspace + "/rpts1"
            arcpy.PointToRaster_conversion(endpoints1, 'grid_code', rpts1, cellsize=cellsize)

            arcpy.AddMessage("Buffering endpoints...")
            rpts11 = Con(rpts1, 1, 0, "Value>0")
            rendbuffers1 = Expand(rpts11, 2*radius, [1])

            # endbuffers1 = workspace + "/endbuffers1"
            # arcpy.Buffer_analysis(endpoints1, endbuffers1, radius, "FULL", "ROUND", "NONE", "")
            # rendbuffers1 = workspace + "/rendbuffers1"
            # arcpy.FeatureToRaster_conversion(endbuffers1, "OBJECTID", rendbuffers1, cellsize)

            mask = CreateConstantRaster(0, "INTEGER", cellsize, dem.extent)
            arcpy.Mosaic_management(mask, rendbuffers1, "MAXIMUM", "FIRST", "", "", "", "0.3", "NONE")

            arcpy.AddMessage("Erasing streams...")
            str1_e = SetNull(rendbuffers1, str1, "value > 0")

            arcpy.AddMessage("Vectorizing erased streams...")
            streams1_e = workspace + "/streams1_e"
            StreamToFeature(str1_e, dir, streams1_e, False)

            arcpy.AddMessage("Deriving erased endpoints...")
            endpoints1_e = workspace + "/endpoints1_e"
            arcpy.FeatureVerticesToPoints_management(streams1_e, endpoints1_e, "END")

            arcpy.AddMessage("Deriving primary watersheds...")
            pour1 = SnapPourPoint(endpoints1_e, acc, cellsize * 1.5, "")

            wsh1 = Watershed(dir, pour1, "")

            arcpy.AddMessage("Vectorizing primary watersheds...")
            watersheds1 = workspace + "/watersheds1"

            try:
                arcpy.RasterToPolygon_conversion(wsh1, watersheds1, False, "")
            except:
                arcpy.AddMessage("\n> FAILED TO CONVERT PRIMARY WATERSHEDS TO VECTOR - PROBABLY EXTENT ERROR. MASKING...\n")

                domain = workspace + "/domain1"

                arcpy.RasterDomain_3d(wsh1, domain, "POLYGON")

                envpoly = workspace + "/envpoly1"
                arcpy.MinimumBoundingGeometry_management(domain, envpoly,
                                                         "ENVELOPE", "NONE")

                envline = workspace + "/envline1"
                arcpy.PolygonToLine_management(envpoly, envline)

                envpolyr = workspace + "/envpolyr1"
                envliner = workspace + "/envliner1"

                arcpy.PolygonToRaster_conversion(envpoly, "OBJECTID",
                                                 envpolyr, "", "", cellsize)

                arcpy.PolylineToRaster_conversion(envline, "OBJECTID",
                                                 envliner, "", "", cellsize)

                envmaskname = "envmask1"
                envmask = workspace + "/" + envmaskname

                envs = [envpolyr, envliner]

                arcpy.MosaicToNewRaster_management(envs,
                                                   workspace,
                                                   envmaskname,
                                                   "",
                                                   "8_BIT_UNSIGNED",
                                                   str(cellsize),
                                                   "1",
                                                   "SUM",
                                                   "FIRST")

                wsh1_m = SetNull(envmask, wsh1, "VALUE = 2")

                arcpy.RasterToPolygon_conversion(wsh1_m, watersheds1, True, "")

            # SECONDARY STREAMS AND WATERSHEDS

            arcpy.AddMessage("PROCESSING SECONDARY STREAMS AND WATERSHEDS")

            arcpy.AddMessage("Extracting secondary streams...")
            str2_0 = rastertinworkspace + "/str2.tif"
            ExtractStreams.execute(acc, str2_0, minacc2, minlen2)

            str2 = SetNull(str2_0, 1, "value = 0")
            str2_e = SetNull(str1_0, str2, "value > 0")
            acc_e = SetNull(str1_0, acc, "value > 0")

            arcpy.AddMessage("Vectorizing streams...")
            streams2_e = workspace + "/streams2_e"
            StreamToFeature(str2_e, dir, streams2_e, False)

            arcpy.AddMessage("Deriving endpoints...")
            endpoints2_e = workspace + "/endpoints2_e"
            arcpy.FeatureVerticesToPoints_management(streams2_e, endpoints2_e, "END")

            arcpy.AddMessage("Buffering primary streams...")
            streambuffer = Expand(str1, radius, [1])
            # streambuffer = workspace + "/streams1_b"
            # arcpy.Buffer_analysis(streams1, streambuffer, radius, "FULL", "ROUND", "NONE", "")

            arcpy.AddMessage("Selecting endpoints...")
            # pointslyr = "points"
            # arcpy.MakeFeatureLayer_management(endpoints2_e, pointslyr)
            # arcpy.SelectLayerByLocation_management(pointslyr, "INTERSECT", streambuffer)
            # pourpts2 = workspace + "/pourpts2"
            # arcpy.CopyFeatures_management(pointslyr, pourpts2)

            rpts2 = workspace + "/rpts2"
            arcpy.PointToRaster_conversion(endpoints2_e, 'grid_code', rpts2, cellsize=cellsize)
            # rpts21 = Con(rpts2, 1, 0, "Value>0")
            pourpts2 = ExtractByMask(rpts2, streambuffer)


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

            try:
                arcpy.RasterToPolygon_conversion(wsh2, watersheds2, False, "")
            except:
                arcpy.AddMessage("\n> FAILED TO CONVERT SECONDARY WATERSHEDS TO VECTOR - PROBABLY EXTENT ERROR. MASKING...\n")
                domain = workspace + "/domain2"

                arcpy.RasterDomain_3d(wsh2, domain, "POLYGON")

                envpoly = workspace + "/envpoly2"
                arcpy.MinimumBoundingGeometry_management(domain, envpoly,
                                                         "ENVELOPE", "NONE")
                envline = workspace + "/envline2"
                arcpy.PolygonToLine_management(envpoly, envline)

                envpolyr = workspace + "/envpolyr2"
                envliner = workspace + "/envliner2"

                arcpy.PolygonToRaster_conversion(envpoly, "OBJECTID",
                                                 envpolyr, "", "", cellsize)

                arcpy.PolylineToRaster_conversion(envline, "OBJECTID",
                                                  envliner, "", "", cellsize)

                envmaskname = "envmask2"
                envmask = workspace + "/" + envmaskname

                envs = [envpolyr, envliner]

                arcpy.MosaicToNewRaster_management(envs,
                                                   workspace,
                                                   envmaskname,
                                                   "",
                                                   "8_BIT_UNSIGNED",
                                                   str(cellsize),
                                                   "1",
                                                   "SUM",
                                                   "FIRST")

                wsh2_m = SetNull(envmask, wsh2, "VALUE = 2")

                arcpy.RasterToPolygon_conversion(wsh2_m, watersheds2, False, "")

            # PROCESSING NARROW AREAS WITHOUT WATERSHEDS

            arcpy.AddMessage("Processing remaining basins...")
            bsn = Basin(dir)
            basins = workspace + "/basins"
            arcpy.RasterToPolygon_conversion(bsn, basins, False, "")
            basins_e = workspace + "/basins_e"
            arcpy.Erase_analysis(basins, watersheds1, basins_e)

            arcpy.AddMessage("Interpolating features into 3D...")

            streams1_3d = workspace + "/streams1_3d"

            arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, streams1, streams1_3d)

            watersheds1_3d = workspace + "/watersheds1_3d"
            arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, watersheds1, watersheds1_3d)

            watersheds2_3d = workspace + "/watersheds2_3d"
            arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, watersheds2, watersheds2_3d)

            basins_e_3d = workspace + "/basins_3d"
            arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, basins_e, basins_e_3d)

            # GENERALIZED TIN SURFACE

            arcpy.AddMessage("DERIVING GENERALIZED SURFACE")

            arcpy.AddMessage("TIN construction...")


            s1 = "'" + streams1_3d + "' Shape.Z " + "hardline"
            w1 = "'" + watersheds1_3d + "' Shape.Z " + "softline"
            w2 = "'" + watersheds2_3d + "' Shape.Z " + "softline"
            b = "'" + basins_e_3d + "' Shape.Z " + "softline"

            features.append(s1)
            features.append(w1)
            features.append(w2)
            features.append(b)
            if process_marine:
                m2 = "'" + clip_3d + "' Shape.Z " + "hardline"
                features.append(m2)

            # arcpy.AddMessage("CLEANING SUPPLEMENTARY DATA")
            # arcpy.Delete_management(str2_0)
            # arcpy.Delete_management(endpoints1)
            # arcpy.Delete_management(endbuffers1)
            # arcpy.Delete_management(rendbuffers1)
            # arcpy.Delete_management(streams1_e)
            # arcpy.Delete_management(endpoints1_e)
            # arcpy.Delete_management(watersheds1)
            # arcpy.Delete_management(streams2_e)
            # arcpy.Delete_management(endpoints2_e)
            # arcpy.Delete_management(streambuffer)
            # arcpy.Delete_management(pointslyr)
            # arcpy.Delete_management(pourpts2)
            # arcpy.Delete_management(watersheds2)
            # arcpy.Delete_management(streams1_3d)
            # arcpy.Delete_management(watersheds1_3d)
            # arcpy.Delete_management(watersheds2_3d)
            # arcpy.Delete_management(str1)
            # arcpy.Delete_management(str2)
            # arcpy.Delete_management(mask)
            # arcpy.Delete_management(str1_e)
            # arcpy.Delete_management(str2_e)
            # arcpy.Delete_management(pour1)
            # arcpy.Delete_management(pour21)
            # arcpy.Delete_management(pour22)
            # arcpy.Delete_management(wsh1)
            # arcpy.Delete_management(wsh2)
            # arcpy.Delete_management(acc_e)

        else:

            arcpy.AddMessage("NO STREAMS DETECTED. PROCESSING BASINS...")
            bsn = Basin(dir)
            basins = workspace + "/basins"
            arcpy.RasterToPolygon_conversion(bsn, basins, False, "")
            basins_3d = workspace + "/basins_3d"
            arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, basins, basins_3d)
            b = "'" + basins_3d + "' Shape.Z " + "softline"
            features.append(b)

            if process_marine:
                m2 = "'" + clip_3d + "' Shape.Z " + "hardline"
                features.append(m2)

            # Temporarily deprecated in favour of Basin Tool
            # arcpy.AddMessage("\n> NO STREAMS DETECTED. USING TIN DECIMATION TO GENERALIZE RASTER\n")
            # zTol = 0.5*(dem.maximum - dem.minimum)
            # arcpy.RasterTin_3d(dem0, tin, zTol)
            #
            # if process_marine:
            #     features = [[clip_3d, 'Shape', '<None>', 'hardline', True]]
            #     arcpy.EditTin_3d(tin, features)

            stream_processing = False

        if flines:
            flines_3d = workspace + "/flines_3d"
            arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, flines, flines_3d)
            f = "'" + flines_3d + "' Shape.Z " + "softline"
            features.append(f)

        if fpolys:
            fpolys_3d = workspace + "/fpolys_3d"
            arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, fpolys, fpolys_3d)
            f = "'" + fpolys_3d + "' Shape.Z " + "softreplace"
            features.append(f)

        featurestring = ';'.join(features)
        arcpy.ddd.CreateTin(tin, "", featurestring, "")

        # GENERALIZED RASTER SURFACE
        arcpy.AddMessage("TIN to raster conversion...")

        rastertin = rastertinworkspace + "/rastertin.tif"
        try:
            arcpy.TinRaster_3d(tin, rastertin, "FLOAT", "NATURAL_NEIGHBORS", "CELLSIZE " + str(cellsize), 1)
        except:
            arcpy.AddMessage("Failed to rasterize TIN using NATURAL_NEIGHBORS method. Switching to linear")
            arcpy.TinRaster_3d(tin, rastertin, "FLOAT", "LINEAR", "CELLSIZE " + str(cellsize), 1)

        # POSTPROCESSING
        arcpy.AddMessage("POSTPROCESSING")

        # Widen valleys and ridges
        widenraster = rastertin
        if is_widen and stream_processing:
            arcpy.AddMessage("Raster widening...")
            widenraster = rastertinworkspace + "/widenraster.tif"
            WidenLandforms.execute(rastertin, streams1, widendist, filtersize, widenraster, widentype, True)

        # Smooth DEM
        result = arcpy.Raster(widenraster)
        if is_smooth:
            arcpy.AddMessage("Raster filtering...")
            neighborhood = NbrRectangle(filtersize, filtersize, "CELL")
            result = FocalStatistics(widenraster, neighborhood, "MEAN", "DATA")

        if process_marine:
            arcpy.AddMessage("Masking marine regions...")
            try:
                result_erased = ExtractByMask(result, cell_erased)
                arcpy.Mosaic_management(result_erased, rastertin, "FIRST", "FIRST", "", "", "", "0.3", "NONE")
                arcpy.AddMessage("Saving result...")
                res = arcpy.Raster(rastertin)
                res.save(scratchworkspace + '/gen/dem' + str(i) + '.tif')
            except:
                arcpy.AddMessage("Retrying to mask in raster mode...")
                erased_r = rastertinworkspace + "/erased_r.tif"

                arcpy.PolygonToRaster_conversion(cell_erased, "OBJECTID",
                                                 erased_r, cellsize=cellsize)

                result_erased = SetNull(erased_r, result, "value is NULL")
                # result_erased = ExtractByMask(result, erased_r)

                arcpy.Mosaic_management(result_erased, rastertin, "FIRST", "FIRST", "", "", "", "0.3", "NONE")
                arcpy.AddMessage("Saving result...")
                res = arcpy.Raster(rastertin)
                res.save(scratchworkspace + '/gen/dem' + str(i) + '.tif')
        else:
            arcpy.AddMessage("Saving result...")
            result.save(scratchworkspace + '/gen/dem' + str(i) + '.tif')

        # arcpy.AddMessage("CLEANING MAIN DATA")

        # if stream_processing:
        #     arcpy.Delete_management(streams1)
        # if process_marine:
        #     arcpy.Delete_management(dem)
        # if is_widen == "true":
        #     arcpy.Delete_management(widenraster)

        # arcpy.Delete_management(fill)
        # arcpy.Delete_management(dir)
        # arcpy.Delete_management(acc)
        # arcpy.Delete_management(str1_0)
        # arcpy.Delete_management(tin)
        # arcpy.Delete_management(rastertin)
        #
        # arcpy.Delete_management(workspace)
        # arcpy.Delete_management(rastertinworkspace)

        # arcpy.CheckInExtension("3D")
        # arcpy.CheckInExtension("Spatial")

        return True

    except:
        arcpy.AddMessage("FAILED TO GENERALIZE tile dem" + str(oid - 1))
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)

        return False

def execute(demdataset,
            output,
            flowdir,
            flowacc,
            flines,
            fpolys,
            cliparea,
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
            is_tiled,
            tile_size,
            is_parallel,
            num_processes,
            is_continued=False,
            continued_folder=None):

    try:
        # arcpy.CheckOutExtension("3D")
        # arcpy.CheckOutExtension("Spatial")

        arcpy.env.overwriteOutput = True

        # ORGANIZE WORKSPACE

        workspace = os.path.dirname(output)
        scratchworkspace = workspace

        if is_continued:
            scratchworkspace = continued_folder
            arcpy.AddMessage('\n> CONTINUING PREVIOUS PROCESSING')
        else:
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

        demsource = arcpy.Raster(demdataset)

        nrows = int(math.ceil(float(demsource.height) / float(tile_size)))
        ncols = int(math.ceil(float(demsource.width) / float(tile_size)))
        total = nrows * ncols

        arcpy.AddMessage(is_tiled)
        is_tiled = is_tiled and total > 1

        cellsize = 0.5 * (demsource.meanCellHeight + demsource.meanCellWidth)

        bufferpixelwidth = math.ceil(max(demsource.width, demsource.height) / (max(nrows, ncols) * 10))
        bufferwidth = bufferpixelwidth * cellsize
        overlap = bufferwidth * 2

        fishnet = workspace + "/fishnet"
        fishbuffer = workspace + "/fishbuffer"
        fishmaskbuffer = workspace + "/fishmaskbuffer"

        oids = []

        if is_continued:

            arcpy.env.workspace = scratchworkspace + '/gen'
            rasters = arcpy.ListRasters("*", "GRID")

            oids_done = [int(raster[3:])+1 for raster in rasters]

            oids = list(set(oids).symmetric_difference(oids_done))

            arcpy.AddMessage('> Already processed oids = ' + str(oids_done))

            arcpy.AddMessage('> Proccessing oids = ' + str(oids))

        else:

            if is_tiled: # SPLIT THE TASK INTO TILES

                arcpy.AddMessage('Splitting raster into ' + str(nrows) + ' x ' + str(ncols) + ' = ' + str(total) + ' tiles')
                arcpy.AddMessage('Tile overlap will be ' + str(2 * bufferpixelwidth) + ' pixels')

                arcpy.AddMessage('Creating fishnet...')

                CreateFishnet.execute(demdataset, fishnet, nrows, ncols, overlap,
                                      split=True, shrink=cellsize)

                arcpy.AddMessage('Creating split buffer...')

                CreateFishnet.execute(demdataset, fishbuffer, nrows, ncols, overlap,
                                      split=False, shrink=cellsize)

                arcpy.AddMessage('Creating mask buffer...')
                mask_overlap = max(demsource.meanCellHeight, demsource.meanCellWidth)

                CreateFishnet.execute(demdataset, fishmaskbuffer, nrows, ncols, overlap,
                                      split=False, overlap2=mask_overlap, shrink=cellsize)

                arcpy.AddMessage('Splitting raster...')
                arcpy.env.extent = demsource.extent # Very important!
                arcpy.env.snapRaster = demsource # Very important!
                arcpy.SplitRaster_management(demdataset,
                                             scratchworkspace + "/source",
                                             'dem',
                                             split_method='NUMBER_OF_TILES',
                                             format='TIFF',
                                             num_rasters = str(ncols) + ' ' + str(nrows),
                                             overlap=2*bufferpixelwidth)


                rows = arcpy.da.SearchCursor(fishnet, 'OID@')
                oids = [int(row[0]) for row in rows]

            else:
                arcpy.env.extent = demsource.extent  # Very important!
                arcpy.env.snapRaster = demsource  # Very important!
                CreateFishnet.execute(demdataset, fishbuffer, 1, 1) # Just DEM envelope
                # arcpy.AddMessage(float(cellsize))
                # mask = CreateConstantRaster(0, "INTEGER", cellsize, demsource.extent)
                # arcpy.Mosaic_management(demdataset, mask, "FIRST", "FIRST")
                # mask.save(scratchworkspace + "/source/dem0")

                arcpy.CopyRaster_management(demdataset, scratchworkspace + "/source/dem0.tif")
                oids = [1]

        # PERFORM PROCESSING

        jobs = []

        if is_parallel and is_tiled:

            nproc = multiprocessing.cpu_count()

            if 0 < num_processes < 1:
                nproc = int(math.ceil(nproc * num_processes))
            elif num_processes < 0:
                nproc = int(math.floor(nproc - num_processes))
                if nproc < 1:
                    nproc = 1
            elif num_processes >= 1:
                nproc = int(math.ceil(num_processes))

            arcpy.AddMessage('> Trying to make multiprocessing using ' + str(nproc) + ' processor cores')
            arcpy.AddMessage('')

            pool = multiprocessing.Pool(nproc)

            args = zip(oids, repeat(demdataset),
                       repeat(flowdir),
                       repeat(flowacc),
                       repeat(flines),
                       repeat(fpolys),
                       repeat(cliparea),
                       repeat(fishbuffer),
                       repeat(minacc1),
                       repeat(minlen1),
                       repeat(minacc2),
                       repeat(minlen2),
                       repeat(is_widen),
                       repeat(widentype),
                       repeat(widendist),
                       repeat(filtersize),
                       repeat(is_smooth),
                       repeat(scratchworkspace))

            results = pool.map(call_list, args)

            pool.close()
            pool.join()

            falseoids = []
            for state, oid in zip(results, oids):
                if state == False:
                    falseoids.append(oid)

            if len(falseoids) > 0:
                arcpy.AddMessage('\n> FAILED to generalize tiles with OID = ' + str(falseoids) + '\n')

        else:
            arcpy.AddMessage('> Processing in sequential mode')
            arcpy.AddMessage('')
            for oid in oids:
                jobs.append(call(oid,
                                 demdataset,
                                 flowdir,
                                 flowacc,
                                 flines,
                                 fpolys,
                                 cliparea,
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
                                 scratchworkspace))
            falseoids = []
            for state, oid in zip(jobs, oids):
                if state == False:
                    falseoids.append(oid)

            if len(falseoids) > 0:
                arcpy.AddMessage('Failed to generalize tiles with OID = ' + str(falseoids))

        # FINALIZE

        if is_tiled:

            arcpy.AddMessage("CLIPPING AND MASKING GENERALIZED RASTERS")
            arcpy.env.snapRaster = demsource

            rows = arcpy.da.SearchCursor(fishmaskbuffer, ['SHAPE@', 'OID@'])
            i = 0
            for row in rows:
                raster = scratchworkspace + '/gen/dem' + str(i) + '.tif'
                if arcpy.Exists(raster):
                    dem = arcpy.Raster(raster)
                    try:
                        dem_clipped = ExtractByMask(dem, row[0])
                        dem_clipped.save(scratchworkspace + '/gencrop/dem' + str(i) + '.tif')
                    except:
                        arcpy.AddMessage('NOTHING TO EXTRACT: Tile ' + str(i) + ' has significant values only outside its frame')
                i += 1

            arcpy.env.workspace = scratchworkspace + '/gencrop'
            rasters = arcpy.ListRasters("*", "TIFF")

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
        else:
            if arcpy.Exists(scratchworkspace + '/gen/dem0.tif'):
                arcpy.AddMessage('SAVING RESULT')
                arcpy.CopyRaster_management(scratchworkspace + '/gen/dem0.tif', output)
            else:
                arcpy.AddMessage('NOTHING GENERALIZED')

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)

    return

if __name__ == '__main__':

    demdataset = arcpy.GetParameterAsText(0)
    output = arcpy.GetParameterAsText(1)
    flowdir = arcpy.GetParameterAsText(2)
    flowacc = arcpy.GetParameterAsText(3)
    flines = arcpy.GetParameterAsText(4)
    fpolys = arcpy.GetParameterAsText(5)
    cliparea = arcpy.GetParameterAsText(6)
    outputcellsize = float(arcpy.GetParameterAsText(7))
    minacc1 = int(arcpy.GetParameterAsText(8))
    minlen1 = int(arcpy.GetParameterAsText(9))
    minacc2 = int(arcpy.GetParameterAsText(10))
    minlen2 = int(arcpy.GetParameterAsText(11))
    is_widen = True if arcpy.GetParameterAsText(12) == 'true' else False
    widentype = arcpy.GetParameterAsText(13)
    widendist = float(arcpy.GetParameterAsText(14))
    filtersize = int(arcpy.GetParameterAsText(15))
    is_smooth = True if arcpy.GetParameterAsText(16) == 'true' else False
    is_tiled = True if arcpy.GetParameterAsText(17) == 'true' else False
    tile_size = arcpy.GetParameterAsText(18)
    is_parallel = True if arcpy.GetParameterAsText(19) == 'true' else False
    num_processes = int(arcpy.GetParameterAsText(20))
    is_continued = True if arcpy.GetParameterAsText(21) == 'true' else False
    continued_folder = arcpy.GetParameterAsText(22)

    execute(demdataset, output, flowdir, flowacc, flines, fpolys, cliparea,
            outputcellsize, minacc1, minlen1, minacc2, minlen2,
            is_widen, widentype, widendist, filtersize,
            is_smooth, is_tiled, tile_size, is_parallel, num_processes,
            is_continued, continued_folder)