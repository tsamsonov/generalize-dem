# -*- coding: cp1251 -*-
import arcpy
import sys
import traceback
import math
import os
import numpy

shift = []
window = []
noData = -9999
# 2017, Timofey Samsonov, Lomonosov Moscow State University


def sample_window(raster, i, j):
    global window
    global shift
    global noData

    w = 0

    for k in shift:
        ik = i+k
        for l in shift:
            jl = j+l
            if raster[ik, jl] != noData:
                window[w] = raster[ik, jl]
                w += 1
    elems = window[0:w]
    return elems


# nfilt = 0
def calc_lower_quartile(raster, i, j):
    elems = sample_window(raster, i, j)
    elems.sort()
    n = int(math.floor(len(elems)/4))
    return 0.5*(elems[n]+elems[n-1])


# nfilt = 1
def calc_upper_quartile(raster, i, j):
    elems = sample_window(raster, i, j)
    elems.sort(reverse=True)
    n = int(math.floor(len(elems)/4))
    return 0.5*(elems[n]+elems[n-1])


# nfilt = 2
def calc_min(raster, i, j):
    elems = sample_window(raster, i, j)
    value = min(elems)
    return value


# nfilt = 3
def calc_max(raster, i, j):
    elems = sample_window(raster, i, j)
    value = max(elems)
    return value


# nfilt = 4
def calc_mean(raster, i, j):
    elems = sample_window(raster, i, j)
    sum = 0
    for k in range(len(elems)):
        sum += elems[k]
    value = sum/len(elems)
    return value


# nfilt = 5
def calc_median(raster, i, j):
    elems = sample_window(raster, i, j)
    elems.sort()
    n = len(elems)
    if n == 0:
        return noData
    elif n%2 == 0:
        return elems[n//2-1]
    else:
        return (elems[n//2-1] + elems[n//2])*0.5


# Filter selector
filters = {0: calc_lower_quartile,
           1: calc_upper_quartile,
           2: calc_min,
           3: calc_max,
           4: calc_mean,
           5: calc_median
}


def extend_array(array, nx, ny, value):
    ni = array.shape[0]
    nj = array.shape[1]

    extarray = numpy.empty((ni + ny, nj + nx))
    extarray.fill(value)
    for i in range(ni):
        for j in range(nj):
            extarray[i, j] = array[i, j]

    return extarray


def process_raster(inraster, niter, nfilt):
    global noData

    ni = inraster.shape[0]
    nj = inraster.shape[1]

    arcpy.SetProgressor("step", "Processing rows", 0, ni-1, 1)
    raster = inraster

    outraster = numpy.zeros((ni, nj))

    for k in range(niter):
        raster = extend_array(raster, len(shift) - 1, len(shift) - 1, noData)

        arcpy.AddMessage("Iteration " + str(k+1))
        for i in range(ni):
            arcpy.SetProgressorLabel("Processing row " + str(i) + " from " + str(ni))
            for j in range(nj):
                if raster[i, j] == noData:
                    outraster[i, j] = raster[i, j]
                else:
                    outraster[i, j] = filters[nfilt](raster, i, j)
            arcpy.SetProgressorPosition(i)
        raster = outraster
    return outraster


def execute(inraster, outraster, wsize, niter, qtype):
    global window
    global shift
    global noData

    shift = [i for i in range(-math.trunc(wsize / 2), math.trunc(wsize / 2) + 1)]

    corrected_size = 2*math.trunc(wsize / 2) + 1
    
    window = [0 for i in range(corrected_size**2)]

    # Select the appropriate filter number

    nfilt = 0
    if qtype == "Lower Quartile":
        nfilt = 0
    elif qtype == "Upper Quartile":
        nfilt = 1
    elif qtype == "Min":
        nfilt = 2
    elif qtype == "Max":
        nfilt = 3
    elif qtype == "Mean":
        nfilt = 4
    elif qtype == "Median":
        nfilt = 5
    else:
        nfilt = 4

    # Process: RasterToASCII_conversion
    rasternumpy = arcpy.RasterToNumPyArray(inraster)

    # Calculating quartiles
    arcpy.AddMessage("Filtering raster...")
    newrasternumpy = process_raster(rasternumpy, niter, nfilt)

    r = arcpy.Raster(inraster)
    lowerleft = arcpy.Point(r.extent.XMin, r.extent.YMin)
    cellsize = r.meanCellWidth
    crs = r.spatialReference

    # Convert python list to ASCII
    arcpy.AddMessage("Writing streams...")
    outinnerraster = arcpy.NumPyArrayToRaster(newrasternumpy, lowerleft, cellsize)
    arcpy.DefineProjection_management(outinnerraster, crs)
    outinnerraster.save(outraster)

if __name__ == "__main__":
    try:
        # Get input parameters
        inRaster = arcpy.GetParameterAsText(0)
        outRaster = arcpy.GetParameterAsText(1)
        wSize = int(arcpy.GetParameterAsText(2))
        nIter = int(arcpy.GetParameterAsText(3))
        qType = arcpy.GetParameterAsText(4)

        execute(inRaster, outRaster, wSize, nIter, qType)

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        
