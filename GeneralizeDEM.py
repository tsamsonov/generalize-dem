__author__ = 'Timofey Samsonov'
# -*- coding: cp1251 -*-
# DEM generalization algorithm
# Important note:  Widen and TraceFlowLines scripts must be in the same directory,
# otherwise Generalize DEM tool will not be able reference them
# 2015, Timofey Samsonov, Lomonosov Moscow State University

import arcpy
from arcpy.sa import *
import os.path, TraceFlowLines, Widen

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")

demdataset = arcpy.GetParameterAsText(0)
output = arcpy.GetParameterAsText(1)
outputcellsize = float(arcpy.GetParameterAsText(2))
minacc1 = long(arcpy.GetParameterAsText(3))
minlen1 = long(arcpy.GetParameterAsText(4))
minacc2 = long(arcpy.GetParameterAsText(5))
minlen2 = long(arcpy.GetParameterAsText(6))
widentype = arcpy.GetParameterAsText(7)
widendist = float(arcpy.GetParameterAsText(8))
filtersize = long(arcpy.GetParameterAsText(9))
workspace = arcpy.GetParameterAsText(10)
is_widen = arcpy.GetParameterAsText(11)
is_smooth = arcpy.GetParameterAsText(12)

if len(workspace) == 0 :
    workspace = os.path.dirname(output)

# raster workspace MUST be a folder, no
rastertinworkspace = workspace
n = len(rastertinworkspace)
if n > 4:
    end = rastertinworkspace[n-4 : n] # extract last 4 letters
    if end == ".gdb": # geodatabase
        rastertinworkspace = os.path.dirname(rastertinworkspace)

scriptpath = os.path.realpath(__file__)
toolboxpath = os.path.dirname(scriptpath)

arcpy.env.scratchWorkspace = rastertinworkspace

dem = Raster(demdataset)

cellsize = dem.meanCellHeight

arcpy.AddMessage("\nPREPROCESSING")

arcpy.AddMessage("Fill...")
fill = Fill(dem, "")
arcpy.AddMessage("Dir...")
dir = FlowDirection(fill,"","")
arcpy.AddMessage("Acc...")
acc = FlowAccumulation(dir,"","INTEGER")

# MAIN STREAMS AND WATERSHEDS
arcpy.AddMessage("\nPROCESSING PRIMARY STREAMS AND WATERSHEDS")

str1_0 = workspace + "/str1"
str2_0 = workspace + "/str2"

arcpy.AddMessage("Extracting primary streams...")

TraceFlowLines.execute(acc, str1_0, minacc1, minlen1, rastertinworkspace)
TraceFlowLines.execute(acc, str2_0, minacc2, minlen2, rastertinworkspace)

str1 = SetNull(str1_0, 1, "value = 0")

arcpy.AddMessage("Processing primary streams...")
streams1 = workspace + "/streams1"
StreamToFeature(str1, dir, streams1, True)

endpoints1 = workspace + "/endpoints1"
arcpy.FeatureVerticesToPoints_management(streams1, endpoints1, "END")

endbuffers1 = workspace + "/endbuffers1"
radius = 2*cellsize
arcpy.Buffer_analysis(endpoints1, endbuffers1, radius, "FULL", "ROUND", "NONE", "")
rendbuffers1 = workspace + "/rendbuffers1"
arcpy.FeatureToRaster_conversion(endbuffers1, "OBJECTID", rendbuffers1, cellsize)

mask = CreateConstantRaster(-1, "INTEGER", cellsize, dem.extent)

arcpy.Mosaic_management(mask, rendbuffers1, "MAXIMUM", "FIRST", "", "", "", "0.3", "NONE")

str1_e = SetNull(rendbuffers1, str1, "value >= 0")

streams1_e = workspace + "/streams1_e"
StreamToFeature(str1_e, dir, streams1_e, True)

endpoints1_e = workspace + "/endpoints1_e"
arcpy.FeatureVerticesToPoints_management(streams1_e, endpoints1_e, "END")

arcpy.AddMessage("Deriving primary watersheds...")
pour1 = SnapPourPoint(endpoints1_e,acc,cellsize*1.5,"")

wsh1 = Watershed(dir, pour1, "")

watersheds1 = workspace + "/watersheds1"
arcpy.RasterToPolygon_conversion(wsh1, watersheds1, True, "")

# SECONDARY STREAMS AND WATERSHEDS

arcpy.AddMessage("\nPROCESSING SECONDARY STREAMS AND WATERSHEDS")

arcpy.AddMessage("Extracting secondary streams...")

str2 = SetNull(str2_0, 1, "value = 0")
str2_e = SetNull(str1_0, str2, "value > 0")
acc_e = SetNull(str1_0, acc, "value > 0")

streams2_e = workspace + "/streams2_e"
StreamToFeature(str2_e, dir, streams2_e, True)

endpoints2_e = workspace + "/endpoints2_e"
arcpy.FeatureVerticesToPoints_management(streams2_e, endpoints2_e, "END")

streambuffer = workspace + "/streams1_b"
arcpy.Buffer_analysis(streams1, streambuffer, radius, "FULL", "ROUND", "NONE", "")

pointslyr = "points"
arcpy.MakeFeatureLayer_management(endpoints2_e, pointslyr)
arcpy.SelectLayerByLocation_management(pointslyr,"INTERSECT",streambuffer)

pourpts2 = arcpy.CreateFeatureclass_management("in_memory", "pourpts2", "POINT", pointslyr, "DISABLED", "DISABLED",
                                               arcpy.Describe(endpoints2_e).spatialReference)
arcpy.CopyFeatures_management(pointslyr, pourpts2)

arcpy.AddMessage("Deriving secondary pour pts 1...")
pour21 = SnapPourPoint(pourpts2,acc_e,cellsize*1.5, "")

arcpy.AddMessage("Deriving secondary pour pts 2...")
pour22 = SnapPourPoint(pourpts2,acc_e,cellsize*2, "")

arcpy.AddMessage("Mosaic secondary pour pts...")
arcpy.Mosaic_management(pour21, pour22, "FIRST", "FIRST", "", "", "", "0.3", "NONE")

wsh2 = Watershed(dir, pour22, "")

arcpy.AddMessage("Deriving secondary watersheds...")
watersheds2 = workspace + "/watersheds2"
arcpy.RasterToPolygon_conversion(wsh2, watersheds2, True, "")

arcpy.AddMessage("Interpolating features into 3D...")

streams1_3d = workspace + "/streams1_3d"
arcpy.InterpolateShape_3d(dem, streams1, streams1_3d)

watersheds1_3d = workspace + "/watersheds1_3d"
arcpy.InterpolateShape_3d(dem, watersheds1, watersheds1_3d)

watersheds2_3d = workspace + "/watersheds2_3d"
arcpy.InterpolateShape_3d(dem, watersheds2, watersheds2_3d)

# GENERALIZED TIN SURFACE

arcpy.AddMessage("\nDERIVING GENERALIZED SURFACE")

arcpy.AddMessage("TIN construction...")

tin = rastertinworkspace + "/tin"
features = []
s1 = "'" + streams1_3d + "' Shape.Z " + "hardline"
w1 = "'" + watersheds1_3d + "' Shape.Z " + "softline"
w2 = "'" + watersheds2_3d + "' Shape.Z " + "softline"
features.append(s1)
features.append(w1)
features.append(w2)

featurestring = ';'.join(features)

arcpy.ddd.CreateTin(tin,"",featurestring,"")

# GENERALIZED RASTER SURFACE
arcpy.AddMessage("TIN to raster conversion...")

rastertin = rastertinworkspace + "/rastertin"
arcpy.TinRaster_3d(tin, rastertin, "FLOAT", "NATURAL_NEIGHBORS", "CELLSIZE " + str(cellsize), 1)

# POSTPROCESSING

arcpy.AddMessage("\nPOSTPROCESSING")

# Widen valleys and ridges
widenraster = rastertin
if is_widen == "true":
    arcpy.AddMessage("Raster widening...")
    widenraster = workspace + "/widenraster"
    Widen.execute(rastertin, streams1, widendist, filtersize, widenraster, widentype)

# Smooth DEM
result = Raster(widenraster)
if is_smooth == "true":
    arcpy.AddMessage("Raster filtering...")
    neighborhood = NbrRectangle(filtersize, filtersize, "CELL")
    result = FocalStatistics(widenraster, neighborhood, "MEAN", "DATA")
    
arcpy.AddMessage("Saving result...")
result.save(output)

arcpy.AddMessage("\nCLEANING TEMPORARY DATA\n")

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

arcpy.CheckInExtension("3D")
arcpy.CheckInExtension("Spatial")







