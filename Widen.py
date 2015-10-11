__author__ = 'Timofey Samsonov'
# -*- coding: cp1251 -*-
# Raster DEM widening by modified Leonowicz-Jenny algorithm
# 2014, Timofey Samsonov, Lomonosov Moscow State University
from arcpy.sa import *
import arcpy, sys, traceback, FilterDEM, time, os

def execute(demdataset, streams, distance, windowsize, output, ftype):
    arcpy.env.workspace = "in_memory"
    arcpy.env.snapRaster = demdataset
    #desc = arcpy.Describe(dem + "/Band_1")
    dem = Raster(demdataset)
    cellsize = dem.meanCellHeight

    # Calculate distances from streams
    arcpy.AddMessage("Calculating distances...")
    distances = EucDistance(streams, "", cellsize, "")

    # Derive valleys weights
    arcpy.AddMessage("Calculating valley weights...")
    divdist = Divide(distances, distance)
    divdistminus = Minus(1, divdist)
    w_valleys = Con(divdistminus, 0, divdistminus, "value < 0")

    # Derive ridges weights
    arcpy.AddMessage("Calculating ridge weights...")
    distminus = Minus(distances, distance)
    distminusdiv = Divide(distminus, distance)
    udistminusdiv = Con(distminusdiv, 0, distminusdiv, "value < 0")
    w_ridges = Con(udistminusdiv, 1, udistminusdiv, "value > 1")

    # Derive source dem weights
    arcpy.AddMessage("Calculating source dem weights...")
    w_both = Plus(w_valleys, w_ridges)
    w_dem = Minus(1, w_both)

    # Get valley and ridge values
    arcpy.AddMessage("Filtering elevations...")

    # unique time identifier
    # t = int(round(time.time() * 1000))
    valleys = None
    ridges = None
    
    if(ftype == "Min/Max"):
        neighborhood = NbrRectangle(windowsize, windowsize, "CELL")
        valleys = FocalStatistics(dem, neighborhood, "MINIMUM", "DATA")
        ridges = FocalStatistics(dem, neighborhood, "MAXIMUM", "DATA")
    else:
        workspace = os.path.dirname(output)
        n = len(workspace)
        if(n > 4):
            end = workspace[n-4 : n] # extract last 4 letters
            if(end == ".gdb"): # geodatabase
                workspace = os.path.dirname(workspace)
                
        val = arcpy.env.workspace + "val"
        rig = arcpy.env.workspace + "rig"
        FilterDEM.execute(demdataset, val, windowsize, 1, "Lower Quartile", workspace)
        FilterDEM.execute(demdataset, rig, windowsize, 1, "Upper Quartile", workspace)
        valleys = Raster(val)
        ridges = Raster(rig)
        
    # Calculate weighted values
    arcpy.AddMessage("Calculating weighted values...")
    weightedvalleys = Times(valleys, w_valleys)
    weightedridges = Times(ridges, w_ridges)
    weighteddem = Times(dem, w_dem)

    # Mix values
    arcpy.AddMessage("Mixing values...")
    mixvalleyridges = Plus(weightedvalleys, weightedridges)
    mixall = Plus(weighteddem, mixvalleyridges)

    # Save the result
    mixall.save(output)
    arcpy.AddMessage("Done")

if __name__ == "__main__":
    try:
        demdataset = arcpy.GetParameterAsText(0)
        streams = arcpy.GetParameterAsText(1)
        distance = float(arcpy.GetParameterAsText(2))
        windowsize = long(arcpy.GetParameterAsText(3))
        output = arcpy.GetParameterAsText(4)
        ftype = arcpy.GetParameterAsText(5)

        execute(demdataset, streams, distance, windowsize, output, ftype)

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)






