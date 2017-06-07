# -*- coding: cp1251 -*-
# Raster stream network generalization by Leonowicz-Jenny algorithm
# 2010, Timofey Samsonov, Lomonosov Moscow State University
import arcpy, sys, traceback, os.path
MAXACC = 0


def findUpCell(accRaster, i, j):
    w = [[0.70710678, 1, 0.70710678],[1, 1, 1], [0.70710678, 1, 0.70710678]] ## distance weights
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
            temp = (accRaster[i][j] - accRaster[ik][jl]) * w[a][b]
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
    acc = accRaster[i][j]
    ik = i
    jk = j
    n = 0
    cell = [i,j]
    
    while n < minLen:
        stream[n] = [ik, jk]
        iUp, jUp = findUpCell(accRaster, ik, jk)
        acc = accRaster[iUp][jUp]
        if acc < minAcc:
            break
        if iUp == ik and jUp == jk:
            break
        ik = iUp
        jk = jUp
        n+=1
    
    if n == minLen:
        for k in range(n):
            streamRaster[stream[k][0]][stream[k][1]] = 1
        while acc > minAcc:
            streamRaster[iUp][jUp] = 1
            iUp, jUp = findUpCell(accRaster, ik, jk)
            if iUp == ik and jUp == jk:
                break
            acc = accRaster[iUp][jUp]
            ik = iUp
            jk = jUp
            n+=1

    return streamRaster 

def extendMatrix(matrix, nx, ny, value):
    
    nullsEnd = [value for i in range(nx)] ## tail for each row
    nullsRow = [value for i in range(len(matrix[0])+nx)] ## additional row template
    nullsMatrix = [nullsRow for i in range(ny)] ## additional rows matrix (bottom)
    for i in range(len(matrix)): 
        matrix[i].extend(nullsEnd) ## Add a tail to each existing row
    matrix.extend(nullsMatrix) # Add new rows to bottom
    return matrix
        
def processRaster(inRaster, outRaster, properties, minAcc, minLen):
    # raster processing here
    arcpy.SetProgressor("step", "Processing rows", 0, int(properties["nrows"])-1, 1)
        
    stream = [[0, 0] for i in range(minLen)]
    ni = int(properties["nrows"])
    nj = int(properties["ncols"])

    global MAXACC
    inRaster = extendMatrix(inRaster, 1, 1, MAXACC*10)

    for i in range(ni):
        arcpy.SetProgressorLabel("Processing row " + str(i) + " from " + properties["nrows"])
        for j in range(nj):
            if inRaster[i][j] > minAcc:
                traceFlowCells(inRaster, outRaster, i, j, stream, minAcc, minLen)
        arcpy.SetProgressorPosition(i)

    return outRaster

def asciiToPythonList(AsciiFile, fieldType):
    #load the raster into a list and this is the output.
    properties = {"nrows":0,
                  "ncols":0,
                  "xllcorner":0,
                  "yllcorner":0,
                  "cellsize":0,
                  "NODATA_value":0
                  }
    rasterAsList = []
    file = open(AsciiFile)
    counter = 0 
    while 1:
        lines = file.readlines(100000)
        if not lines:
            break
        for line in lines:
            counter += 1
            if counter == 1:
                properties["ncols"] = line.split(" ")[-1]
            if counter == 2:
                properties["nrows"] = line.split(" ")[-1]
                arcpy.SetProgressor("step", "Reading rows", 0, int(properties["nrows"])-1, 1)
            if counter == 3:
                properties["xllcorner"] = line.split(" ")[-1]
            if counter == 4:
                properties["yllcorner"] = line.split(" ")[-1]
            if counter == 5:
                properties["cellsize"] = line.split(" ")[-1]
            if counter == 6:
                properties["NODATA_value"] = line.split(" ")[-1]
            if counter > 6:
                rasterAsList.append(line.split(" "))
                arcpy.SetProgressorLabel("Reading row " + str(counter-7) + " from " + properties["nrows"])
                for y in range(int(properties["ncols"])):
                    if fieldType == "INTEGER":
                        rasterAsList[counter-7][y] = int(rasterAsList[counter-7][y])
                    else:
                        rasterAsList[counter-7][y] = float(rasterAsList[counter-7][y])

                ## Delete redundant cell in a row
                while 1:
                    if len(rasterAsList[counter-7]) > int(properties["ncols"]):
                        rasterAsList[counter-7].pop()
                    else:
                        break
                    
                arcpy.SetProgressorPosition(counter-7)
    file.close()

    return rasterAsList, properties

def generateBlankRaster(properties):
    rasterAsList2 = []
    for x in range(int(properties["nrows"])):
        row = [0 for y in range(int(properties["ncols"]))]
        rasterAsList2.append(row)

    return rasterAsList2

def pythonListToASCII(rasterAsList, properties, workspace):
    ModAsciiFile = arcpy.CreateUniqueName("modified.asc", workspace)
    #write out the modified raster to a text file
    file = open(ModAsciiFile, "w")
    file.write("ncols         " + properties["ncols"])
    file.write("nrows         " + properties["nrows"])
    file.write("xllcorner     " + properties["xllcorner"])
    file.write("yllcorner     " + properties["yllcorner"])
    file.write("cellsize      " + properties["cellsize"])
    file.write("NODATA_value  " + properties["NODATA_value"])

    for x in range(int(properties["nrows"])):
        for x in rasterAsList[x]:
            file.write(str(x) + " ")
    file.close()
    
    return ModAsciiFile

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

    fieldType = "INTEGER"
    VALUETYPE = int(str(arcpy.GetRasterProperties_management(inRaster, "VALUETYPE")))
    if 8 < VALUETYPE < 11:
        fieldType = "FLOAT"
    elif VALUETYPE >= 11:
        arcpy.AddMessage("Unsupported field format: " + str(VALUETYPE))
        fieldType = "FLOAT"

    # Process: RasterToASCII_conversion
    arcpy.AddMessage("Reading raster...")
    InAsciiFile = arcpy.CreateUniqueName("temp.asc", workspace)

    arcpy.RasterToASCII_conversion(inRaster, InAsciiFile)

    # Convert ASCII to python list
    rasterAsList, properties = asciiToPythonList(InAsciiFile, fieldType)
    arcpy.AddMessage("Reading completed")
    newRasterList = generateBlankRaster(properties)

    # Tracing stream lines
    arcpy.AddMessage("Tracing stream lines...")
    newRasterList = processRaster(rasterAsList, newRasterList, properties, minAcc, minLen)

    # Convert python list to ASCII
    arcpy.SetProgressor("default", "Finishing...")
    arcpy.AddMessage("Writing streams...")
    ModAsciiFile = pythonListToASCII(newRasterList, properties, workspace)

    # Write ASCII to output raster
    arcpy.ASCIIToRaster_conversion(ModAsciiFile, outRaster, fieldType)
    arcpy.Delete_management(ModAsciiFile)
    arcpy.Delete_management(InAsciiFile)

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
        
