# -*- coding: cp1251 -*-
import arcpy, sys, traceback, math, os
shift = []
window = []
noData = -9999
# 2015, Timofey Samsonov, Lomonosov Moscow State University

def sampleWindow(raster, i, j):
    global window
    global shift
    global noData

    w = 0

    for k in shift:
        ik = i+k
        for l in shift:
            jl = j+l
            if raster[ik][jl] != noData:
                window[w] = raster[ik][jl]
                w+=1
    elems = window[0:w]
    return elems

# nfilt = 0
def calcLowerQuartile(raster, i, j):
    elems = sampleWindow(raster, i, j)
    elems.sort()
    n = int(math.floor(len(elems)/4))
    return 0.5*(elems[n]+elems[n-1])

# nfilt = 1
def calcUpperQuartile(raster, i, j):
    elems = sampleWindow(raster, i, j)
    elems.sort(reverse = True)
    n = int(math.floor(len(elems)/4))
    return 0.5*(elems[n]+elems[n-1])

# nfilt = 2
def calcMin(raster, i, j):
    elems = sampleWindow(raster, i, j)
    value = min(elems)
    return value

# nfilt = 3
def calcMax(raster, i, j):
    elems = sampleWindow(raster, i, j)
    value = max(elems)
    return value

# nfilt = 4
def calcMean(raster, i, j):
    elems = sampleWindow(raster, i, j)
    sum = 0
    for k in range(len(elems)):
        sum += elems[k]
    value = sum/len(elems)
    return value

# nfilt = 5
def calcMedian(raster, i, j):
    elems = sampleWindow(raster, i, j)
    elems.sort()
    n = len(elems)
    if(n==0):
        return noData
    elif(n%2 == 0):
        return elems[n//2-1]
    else:
        return (elems[n//2-1] + elems[n//2])*0.5

# Filter selector
filters = { 0 : calcLowerQuartile,
            1 : calcUpperQuartile,
            2 : calcMin,
            3 : calcMax,
            4 : calcMean,
            5 : calcMedian
}


## this function extends matrix by adding new rows and columns and setting a value for new cells
def extendMatrix(matrix, nx, ny, value):
    nullsEnd = [value for i in range(nx)] ## tail for each row
    nullsRow = [value for i in range(len(matrix[0])+nx)] ## additional row template
    nullsMatrix = [nullsRow for i in range(ny)] ## additional rows matrix (bottom)
    for i in range(len(matrix)): 
        matrix[i].extend(nullsEnd) ## Add a tail to each existing row
    matrix.extend(nullsMatrix) # Add new rows to bottom
    return matrix
        
def processRaster(inRaster, properties, nIter, nFilt):
    global noData

    nrows = int(properties["nrows"])
    arcpy.SetProgressor("step", "Processing rows", 0, nrows-1, 1)
    raster = inRaster

    outRaster = None

    for k in range(nIter):
        raster = extendMatrix(raster, len(shift)-1, len(shift)-1, noData)
        outRaster = generateBlankRaster(properties)
        arcpy.AddMessage("Iteration " + str(k+1))
        for i in range(int(properties["nrows"])):
            arcpy.SetProgressorLabel("Processing row " + str(i) + " from " + str(nrows))
            for j in range(int(properties["ncols"])):
                if raster[i][j] == noData:
                    outRaster[i][j] = raster[i][j]
                else:
                    outRaster[i][j] = filters[nFilt](raster, i, j)
            arcpy.SetProgressorPosition(i)
        raster = outRaster
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
    raster = []
    for x in range(int(properties["nrows"])):
        row = [0 for y in range(int(properties["ncols"]))]
        raster.append(row)
    return raster

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

def execute(inRaster, outRaster, wSize, nIter, qType, workspace):
    global window
    global shift
    global noData

    n = len(workspace)
    if(n==0):
        workspace = os.path.dirname(outRaster)
        n = len(workspace)
    if(n > 4):
        end = workspace[n-4 : n] # extract last 4 letters
        if(end == ".gdb"): # geodatabase
            workspace = os.path.dirname(workspace)

    shift = [i for i in range(-math.trunc(wSize/2), math.trunc(wSize/2)+1)]

    corrected_size = 2*math.trunc(wSize/2)+1
    
    window = [0 for i in range(corrected_size**2)]

    fieldType = "INTEGER"
    VALUETYPE = int(str(arcpy.GetRasterProperties_management(inRaster, "VALUETYPE")))
    if 8 < VALUETYPE < 11:
        fieldType = "FLOAT"
    elif VALUETYPE >= 11:
        arcpy.AddMessage("Unsupported field format: " + str(VALUETYPE) + "\nSetting to FLOAT")
        fieldType = "FLOAT"

    # Select the appropriate filter number
    nFilt = 0
    if(qType == "Lower Quartile"):
        nFilt = 0
    elif(qType == "Upper Quartile"):
        nFilt = 1
    elif(qType == "Min"):
        nFilt = 2
    elif(qType == "Max"):
        nFilt = 3
    elif(qType == "Mean"):
        nFilt = 4
    elif(qType == "Median"):
        nFilt = 5
    else:
        nFilt = 4

    # Process: RasterToASCII_conversion
    arcpy.AddMessage("Reading raster...")
    InAsciiFile = arcpy.CreateUniqueName("temp.asc", workspace)

    arcpy.RasterToASCII_conversion(inRaster, InAsciiFile)

    # Convert ASCII to python list
    rasterAsList, properties = asciiToPythonList(InAsciiFile, fieldType)
    arcpy.AddMessage("Reading completed")

    if(fieldType == "INTEGER"):
        noData = int(properties["NODATA_value"])
    else:
        noData = float(properties["NODATA_value"])


    # Calculating quartiles
    arcpy.AddMessage("Filtering raster...")
    newRasterList = processRaster(rasterAsList, properties, nIter, nFilt)

    # Convert python list to ASCII
    arcpy.SetProgressor("default", "Finishing...")
    arcpy.AddMessage("Writing output raster...")
    ModAsciiFile = pythonListToASCII(newRasterList, properties, workspace)

    # Write ASCII to output raster
    arcpy.ASCIIToRaster_conversion(ModAsciiFile, outRaster)
    arcpy.Delete_management(InAsciiFile)

if __name__ == "__main__":
    try:
        # Get input parameters
        inRaster = arcpy.GetParameterAsText(0)
        outRaster = arcpy.GetParameterAsText(1)
        wSize = int(arcpy.GetParameterAsText(2))
        nIter = int(arcpy.GetParameterAsText(3))
        qType = arcpy.GetParameterAsText(4)
        workspace = arcpy.GetParameterAsText(5)

        execute(inRaster, outRaster, wSize, nIter, qType, workspace)

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        
