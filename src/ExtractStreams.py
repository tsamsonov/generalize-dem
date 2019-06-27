import numpy
import sys
import os
import arcpy

# print(os.path.dirname(os.path.abspath(__file__)))
# sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/modules")

import StreamExtractor

def process_raster_cpp(inraster, minacc, minlen):
    nrow = inraster.shape[0]
    ncol = inraster.shape[1]

    outraster = numpy.zeros((nrow, ncol))

    return StreamExtractor.ExtractStreams(inraster, outraster, minacc, minlen)

def execute(inraster, outraster, minacc, minlen):
    global MAXACC

    arcpy.env.overwriteOutput = True
    MAXACC = float(str(arcpy.GetRasterProperties_management(inraster, "MAXIMUM")))

    rasternumpy = arcpy.RasterToNumPyArray(inraster, nodata_to_value = MAXACC + 1)

    # Tracing stream lines
    arcpy.AddMessage("Tracing stream lines...")
    newrasternumpy = process_raster_cpp(rasternumpy, minacc, minlen)

    desc = arcpy.Describe(inraster)
    lowerleft = arcpy.Point(desc.extent.XMin, desc.extent.YMin)
    cellsize = desc.meanCellWidth
    crs = desc.spatialReference

    # Convert python list to ASCII
    arcpy.AddMessage("Writing streams...")
    outinnerraster = arcpy.NumPyArrayToRaster(newrasternumpy, lowerleft, cellsize)
    arcpy.DefineProjection_management(outinnerraster, crs)
    outinnerraster.save(outraster)

if __name__ == "__main__":
    execute("W:/Relief/Default.gdb/acc",
            "W:/Relief/Default.gdb/str",
            200, 20)
