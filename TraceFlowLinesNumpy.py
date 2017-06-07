# -*- coding: cp1251 -*-
# Raster stream network generalization by Leonowicz-Jenny algorithm
# 2010, Timofey Samsonov, Lomonosov Moscow State University
import arcpy, numpy, sys, traceback, os.path

MAXACC = 0

def findUpCell(accRaster, i, j):
    w = [[0.70710678, 1, 0.70710678],[1, 1, 1], [0.70710678, 1, 0.70710678]] ## distance weights
    shift = [-1, 0, 1]
    minmax = 4000000000
    a = 0
    kmin = 1
    lmin = 1

    ## finding differences in 3x3 neighbourhood

    for k in shift:
        b = 0
        ik = i+k
        for l in shift:
            jl = j+l
            temp = (accRaster[i, j] - accRaster[ik, jl]) * w[a][b]
            if 0 < temp < minmax:
                minmax = temp 
                kmin = a
                lmin = b      
            b+=1
        a+=1 

    iUp = i + shift[kmin]
    jUp = j + shift[lmin]

    return iUp, jUp

def traceFlowCells(accRaster, streamRaster, i, j, stream, minAcc, minLen):
    acc = accRaster[i, j]
    ik = i
    jk = j
    n = 0
    
    while n < minLen:
        stream[n] = [ik, jk]
        iUp, jUp = findUpCell(accRaster, ik, jk)
        acc = accRaster[iUp, jUp]
        if acc < minAcc:
            break
        if iUp == ik and jUp == jk:
            break
        ik = iUp
        jk = jUp
        n+=1
    
    if n == minLen:
        for k in range(n):
            streamRaster[stream[k][0], stream[k][1]] = 1
        while acc > minAcc:
            streamRaster[iUp, jUp] = 1
            iUp, jUp = findUpCell(accRaster, ik, jk)
            if iUp == ik and jUp == jk:
                break
            acc = accRaster[iUp, jUp]
            ik = iUp
            jk = jUp
            n+=1

    return streamRaster 

def extendArray(array, nx, ny, value):
    
    ni = array.shape[0]
    nj = array.shape[1]

    extArray = numpy.empty((ni+ny, nj+nx))
    extArray.fill(value)
    for i in range(ni):
        for j in range(nj):
            extArray[i, j] = array[i, j]

    return extArray
        
def processRaster(inRaster, minAcc, minLen):
    # raster processing here

    arcpy.AddMessage("Streaming...")
    stream = [[0, 0] for i in range(minLen)]
    ni = inRaster.shape[0]
    nj = inRaster.shape[1]

    global MAXACC

    arcpy.AddMessage("Zeroing...")
    outRaster = numpy.zeros((ni, nj))

    arcpy.AddMessage("Extending...")
    extInRaster = extendArray(inRaster, 1, 1, MAXACC * 10)

    arcpy.AddMessage("Tracing...")

    arcpy.SetProgressor("step", "Processing rows", 0, ni - 1, 1)
    for i in range(ni):
        arcpy.SetProgressorLabel("Processing row " + str(i) + " from " + str(ni))
        for j in range(nj):
            if inRaster[i, j] > minAcc:
                traceFlowCells(extInRaster, outRaster, i, j, stream, minAcc, minLen)
        arcpy.SetProgressorPosition(i)

    return outRaster


def execute(inRaster, outRaster, minAcc, minLen, workspace):
    # scratch workspace MUST be a folder, not a geodatabase
    n = len(workspace)
    if(n==0):
        workspace = os.path.dirname(outRaster)
        n = len(workspace)
    if(n > 4):
        end = workspace[n-4 : n] # extract last 4 letters
        if(end == ".gdb"): # geodatabase
            workspace = os.path.dirname(workspace)

    global MAXACC
    MAXACC = float(str(arcpy.GetRasterProperties_management(inRaster, "MAXIMUM")))

    rasterNumpy = arcpy.RasterToNumPyArray(inRaster)

    # Tracing stream lines
    arcpy.AddMessage("Tracing stream lines...")
    newRasterNumpy = processRaster(rasterNumpy, minAcc, minLen)

    r = arcpy.Raster(inRaster)
    lowerLeft = arcpy.Point(r.extent.XMin, r.extent.YMin)
    cellSize = r.meanCellWidth
    crs = r.spatialReference

    # Convert python list to ASCII
    arcpy.AddMessage("Writing streams...")
    outInnerRaster = arcpy.NumPyArrayToRaster(newRasterNumpy, lowerLeft, cellSize)
    arcpy.DefineProjection_management(outInnerRaster, crs)
    outInnerRaster.save(outRaster)

if __name__ == "__main__":
    try:

        # Get input parameters
        inRaster = arcpy.GetParameterAsText(0)
        outRaster = arcpy.GetParameterAsText(1)
        minAcc = int(arcpy.GetParameterAsText(2))
        minLen = int(arcpy.GetParameterAsText(3))
        workspace = arcpy.GetParameterAsText(4)

        # Execute processing
        execute(inRaster, outRaster, minAcc, minLen, workspace)
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        
