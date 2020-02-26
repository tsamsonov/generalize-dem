# -*- coding: cp1251 -*-
# Automated DEM conflation with reference hydrographic lines
# 2020, Timofey Samsonov, Lomonosov Moscow State University
import os
import arcpy
import sys
import traceback
from arcpy.sa import *
import CounterpartStreams as CS

__author__ = 'Timofey Samsonov'


def execute(inraster, hydrolines, inidfield, outraster, minacc, radius):

    arcpy.CheckOutExtension("Spatial")
    arcpy.CheckOutExtension("3D")

    workspace = os.path.dirname(outraster)

    arcpy.env.snapRaster = inraster
    dem = arcpy.Raster(inraster)
    cellsize = dem.meanCellHeight

    outstreams = 'in_memory/streams'

    CS.execute(hydrolines, inidfield, inraster, outstreams, minacc, radius)

if __name__ == "__main__":
    try:
        inraster = arcpy.GetParameterAsText(0)
        hydrolines = arcpy.GetParameterAsText(1)
        inidfield = arcpy.GetParameterAsText(2)
        outraster = arcpy.GetParameterAsText(2)
        minacc = float(arcpy.GetParameterAsText(4))
        radius = float(arcpy.GetParameterAsText(3))

        execute(inraster, hydrolines, inidfield, outraster, minacc, radius)
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)