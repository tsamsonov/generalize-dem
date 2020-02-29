# -*- coding: cp1251 -*-
# Automated DEM conflation with reference hydrographic lines
# 2020, Timofey Samsonov, Lomonosov Moscow State University
import os
import arcpy
import sys
import traceback
from arcpy.sa import *
import ScratchWorkspace as SW

__author__ = 'Timofey Samsonov'

def execute(in_raster, in_links, in_area, distance, out_raster):
    desc = arcpy.Describe(in_raster)
    cellsize = arcpy.GetRasterProperties_management(in_raster, "CELLSIZEX").getOutput(0)
    crs = desc.spatialReference
    workspace = os.path.dirname(out_raster)

    scratchworkspace = SW.CreateScratchWorkspace(workspace)

    arcpy.AddMessage("Converting raster to points...")
    dempts = 'in_memory/dempts'
    arcpy.RasterToPoint_conversion(in_raster, dempts)

    arcpy.AddMessage("Preparing points and links...")
    buf = 'in_memory/buf'
    arcpy.Buffer_analysis(in_area, buf, distance, dissolve_option='ALL')
    arcpy.Densify_edit(buf, 'DISTANCE', cellsize)

    identity_links = 'in_memory/idlinks'
    arcpy.FeatureVerticesToPoints_management(buf, identity_links)

    ptslyr = 'ptslyr'
    arcpy.MakeFeatureLayer_management(dempts, ptslyr)
    arcpy.SelectLayerByLocation_management(ptslyr, 'INTERSECT', buf)

    confpts = 'in_memory/confpts'
    arcpy.CopyFeatures_management(ptslyr, confpts)

    arcpy.AddMessage("Rubbersheeting...")
    arcpy.RubbersheetFeatures_edit(confpts, in_links, identity_links, 'NATURAL_NEIGHBOR')

    arcpy.AddMessage("Triangulation...")
    arcpy.SelectLayerByAttribute_management(ptslyr, 'SWITCH_SELECTION')

    features = []
    features.append("'" + confpts + "' grid_code " + "masspoints")
    features.append("'" + ptslyr + "' grid_code " + "masspoints")
    featurestring = ';'.join(features)

    tin = scratchworkspace + '/tin'
    arcpy.CreateTin_3d(tin, crs, featurestring)

    arcpy.AddMessage("Converting to output raster...")
    arcpy.env.extent = desc.extent  # Very important!
    arcpy.env.snapRaster = in_raster

    try:
        arcpy.TinRaster_3d(tin, out_raster, "FLOAT", "NATURAL_NEIGHBORS", "CELLSIZE " + str(cellsize), 1)
    except:
        arcpy.AddMessage("Failed to rasterize TIN using NATURAL_NEIGHBORS method. Switching to linear")
        arcpy.TinRaster_3d(tin, out_raster, "FLOAT", "LINEAR", "CELLSIZE " + str(cellsize), 1)

if __name__ == "__main__":
    try:
        in_raster = arcpy.GetParameterAsText(0)
        in_links = arcpy.GetParameterAsText(1)
        in_area = arcpy.GetParameterAsText(2)
        distance = float(arcpy.GetParameterAsText(3))
        out_raster = arcpy.GetParameterAsText(4)

        execute(in_raster, in_links, in_area, distance, out_raster)
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)