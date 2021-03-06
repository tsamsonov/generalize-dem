# -*- coding: cp1251 -*-
# Raster stream network generalization by Leonowicz-Jenny algorithm
# 2017-2019, Timofey Samsonov, Lomonosov Moscow State University
import sys
import arcpy
import numpy
import traceback

can_use_cpp = True

if sys.version_info[:2] == (2, 7):  # ArcGIS for Desktop 10.3+, Python 2.7
    try:
        import StreamExtractor
    except:
        can_use_cpp = False
elif sys.version_info[:2] == (3, 6):  # ArcGIS Pro, Python 3.6
    try:
        import StreamExtractor3 as StreamExtractor
    except:
        can_use_cpp = False
else:
    can_use_cpp = False

MAXACC = 0

def find_up_cell(accraster, i, j):
    w = [[0.70710678, 1, 0.70710678],[1, 1, 1], [0.70710678, 1, 0.70710678]]  # distance weights
    shift = [-1, 0, 1]
    minmax = 4000000000
    a = 0
    kmin = 1
    lmin = 1

    # finding differences in 3x3 neighbourhood

    for k in shift:
        b = 0
        ik = i+k
        for l in shift:
            jl = j+l
            temp = (accraster[i, j] - accraster[ik, jl]) * w[a][b]
            if 0 < temp < minmax:
                minmax = temp 
                kmin = a
                lmin = b      
            b += 1
        a += 1

    iUp = i + shift[kmin]
    jUp = j + shift[lmin]

    return iUp, jUp


def trace_flow_cells(accraster, streamraster, i, j, stream, minacc, minlen):
    acc = accraster[i, j]
    ik = i
    jk = j
    n = 0

    try:
        while n < minlen:
            stream[n] = [ik, jk]
            iup, jup = find_up_cell(accraster, ik, jk)

            acc = accraster[iup, jup]
            if acc < minacc:
                break
            if iup == ik and jup == jk:
                break

            ik = iup
            jk = jup
            n += 1

        if n == minlen:
            for k in range(n):
                streamraster[stream[k][0], stream[k][1]] = 1
            while acc > minacc:
                streamraster[iup, jup] = 1
                iup, jup = find_up_cell(accraster, ik, jk)
                if iup == ik and jup == jk:
                    break
                acc = accraster[iup, jup]
                ik = iup
                jk = jup
                n += 1

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        raise Exception

    return streamraster


def extend_array(array, nx, ny, value):
    
    ni = array.shape[0]
    nj = array.shape[1]

    extarray = numpy.empty((ni+ny, nj+nx))
    extarray.fill(value)
    for i in range(ni):
        for j in range(nj):
            extarray[i, j] = array[i, j]

    return extarray


def process_raster(inraster, minacc, minlen):

    try:
        stream = [[0, 0] for i in range(minlen)]
        ni = inraster.shape[0]
        nj = inraster.shape[1]

        global MAXACC

        outraster = numpy.zeros((ni, nj))
        extinraster = extend_array(inraster, 1, 1, MAXACC + 1)
        arcpy.SetProgressor("step", "Processing rows", 0, ni - 1, 1)
        for i in range(ni):
            arcpy.SetProgressorLabel("Processing row " + str(i) + " from " + str(ni))
            for j in range(nj):
                if  minacc < inraster[i, j] <= MAXACC:
                    trace_flow_cells(extinraster, outraster, i, j, stream, minacc, minlen)
            arcpy.SetProgressorPosition(i)

        return outraster
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        raise Exception

def process_raster_cpp(inraster, minacc, minlen):

    nrow = inraster.shape[0]
    ncol = inraster.shape[1]

    outraster = numpy.zeros((nrow, ncol))

    return StreamExtractor.extract_streams(inraster, outraster, minacc, minlen)

def execute(inraster, outraster, minacc, minlen):
    global MAXACC
    MAXACC = float(str(arcpy.GetRasterProperties_management(inraster, "MAXIMUM")))

    rasternumpy = arcpy.RasterToNumPyArray(inraster, nodata_to_value = MAXACC + 1)

    # Tracing stream lines
    arcpy.AddMessage("Tracing stream lines...")
    newrasternumpy = process_raster_cpp(rasternumpy, minacc, minlen) if can_use_cpp \
                else process_raster(rasternumpy, minacc, minlen)

    desc = arcpy.Describe(inraster)
    lowerleft = arcpy.Point(desc.extent.XMin, desc.extent.YMin)
    cellsize = desc.meanCellWidth
    crs = desc.spatialReference

    arcpy.AddMessage(cellsize)

    # Convert python list to ASCII
    arcpy.AddMessage("Writing streams...")
    outinnerraster = arcpy.NumPyArrayToRaster(newrasternumpy, lowerleft, cellsize, cellsize)
    arcpy.DefineProjection_management(outinnerraster, crs)

    weird_width = outinnerraster.meanCellWidth
    arcpy.AddMessage(weird_width)

    rescaled = 'in_memory/rescaled'
    arcpy.Rescale_management(outinnerraster, rescaled, cellsize/weird_width, cellsize/weird_width)
    arcpy.Shift_management(rescaled, outraster, 0, -2*cellsize, inraster)
    arcpy.AddMessage(arcpy.sa.Raster(outraster).meanCellWidth)

    # outinnerraster.save(outraster)

if __name__ == "__main__":
    try:
        inRaster = arcpy.GetParameterAsText(0)
        outRaster = arcpy.GetParameterAsText(1)
        minAcc = float(arcpy.GetParameterAsText(2))
        minLen = int(arcpy.GetParameterAsText(3))

        execute(inRaster, outRaster, minAcc, minLen)
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)