# -*- coding: cp1251 -*-
import arcpy
import sys
import traceback

# 2020, Timofey Samsonov, Lomonosov Moscow State University


def execute(in_rasters, in_subdivision, out_raster):
    return

if __name__ == "__main__":
    try:
        in_rasters = arcpy.GetParameterAsText(0)
        in_subdivision = arcpy.GetParameterAsText(1)
        out_raster = arcpy.GetParameterAsText(2)

        execute(in_rasters, in_subdivision, out_raster)
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)