# -*- coding: cp1251 -*-
# Raster stream network generalization by Leonowicz-Jenny algorithm
# 2017, Timofey Samsonov, Lomonosov Moscow State University
import arcpy, numpy, sys, traceback

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

    arcpy.AddMessage("Streaming...")
    stream = [[0, 0] for i in range(minlen)]
    ni = inraster.shape[0]
    nj = inraster.shape[1]

    global MAXACC

    arcpy.AddMessage("Zeroing...")
    outraster = numpy.zeros((ni, nj))

    arcpy.AddMessage("Extending...")
    extinraster = extend_array(inraster, 1, 1, MAXACC * 10)

    arcpy.AddMessage("Tracing...")

    arcpy.SetProgressor("step", "Processing rows", 0, ni - 1, 1)
    for i in range(ni):
        arcpy.SetProgressorLabel("Processing row " + str(i) + " from " + str(ni))
        for j in range(nj):
            if inraster[i, j] > minacc:
                trace_flow_cells(extinraster, outraster, i, j, stream, minacc, minlen)
        arcpy.SetProgressorPosition(i)

    return outraster


def execute(inraster, outraster, minacc, minlen):
    global MAXACC
    MAXACC = float(str(arcpy.GetRasterProperties_management(inraster, "MAXIMUM")))

    rasternumpy = arcpy.RasterToNumPyArray(inraster)

    # Tracing stream lines
    arcpy.AddMessage("Tracing stream lines...")
    newrasternumpy = process_raster(rasternumpy, minacc, minlen)

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
    try:

        arcpy.AddMessage('Reading parameters')
        # Get input parameters
        inRaster = arcpy.GetParameterAsText(0)
        outRaster = arcpy.GetParameterAsText(1)
        minAcc = float(arcpy.GetParameterAsText(2))
        minLen = long(arcpy.GetParameterAsText(3))

        arcpy.AddMessage('Executing')

        # Execute processing
        execute(inRaster, outRaster, minAcc, minLen)
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
