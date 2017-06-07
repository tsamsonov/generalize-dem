# -*- coding: cp1251 -*-
import arcpy, sys, traceback, math, os
shift = []
window = []
noData = -9999
# 2015, Timofey Samsonov, Lomonosov Moscow State University


def sample_window(raster, i, j):
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
def calc_lower_quartile(raster, i, j):
    elems = sample_window(raster, i, j)
    elems.sort()
    n = int(math.floor(len(elems)/4))
    return 0.5*(elems[n]+elems[n-1])


# nfilt = 1
def calc_upper_quartile(raster, i, j):
    elems = sample_window(raster, i, j)
    elems.sort(reverse = True)
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
    if(n==0):
        return noData
    elif(n%2 == 0):
        return elems[n//2-1]
    else:
        return (elems[n//2-1] + elems[n//2])*0.5


# Filter selector
filters = { 0 : calc_lower_quartile,
            1 : calc_upper_quartile,
            2 : calc_min,
            3 : calc_max,
            4 : calc_mean,
            5 : calc_median
}


# this function extends matrix by adding new rows and columns and setting a value for new cells
def extend_matrix(matrix, nx, ny, value):
    nullsend = [value for i in range(nx)]  # tail for each row
    nullsrow = [value for i in range(len(matrix[0])+nx)]  # additional row template
    nullsmatrix = [nullsrow for i in range(ny)]  # additional rows matrix (bottom)
    for i in range(len(matrix)): 
        matrix[i].extend(nullsend)  # Add a tail to each existing row
    matrix.extend(nullsmatrix)  # Add new rows to bottom
    return matrix

        
def process_raster(inraster, properties, niter, nfilt):
    global noData

    nrows = int(properties["nrows"])
    arcpy.SetProgressor("step", "Processing rows", 0, nrows-1, 1)
    raster = inraster

    outraster = None

    for k in range(niter):
        raster = extend_matrix(raster, len(shift) - 1, len(shift) - 1, noData)
        outraster = generate_blank_raster(properties)
        arcpy.AddMessage("Iteration " + str(k+1))
        for i in range(int(properties["nrows"])):
            arcpy.SetProgressorLabel("Processing row " + str(i) + " from " + str(nrows))
            for j in range(int(properties["ncols"])):
                if raster[i][j] == noData:
                    outraster[i][j] = raster[i][j]
                else:
                    outraster[i][j] = filters[nfilt](raster, i, j)
            arcpy.SetProgressorPosition(i)
        raster = outraster
    return outraster


def ascii_to_python_list(AsciiFile, fieldType):

    # load the raster into a list and this is the output
    properties = {"nrows":0,
                  "ncols":0,
                  "xllcorner":0,
                  "yllcorner":0,
                  "cellsize":0,
                  "NODATA_value":0
                  }
    rasteraslist = []
    infile = open(AsciiFile)
    counter = 0 
    while 1:
        lines = infile.readlines(100000)
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
                rasteraslist.append(line.split(" "))
                arcpy.SetProgressorLabel("Reading row " + str(counter-7) + " from " + properties["nrows"])
                for y in range(int(properties["ncols"])):
                    if fieldType == "INTEGER":
                        rasteraslist[counter-7][y] = int(rasteraslist[counter-7][y])
                    else:
                        rasteraslist[counter-7][y] = float(rasteraslist[counter-7][y])

                # Delete redundant cell in a row
                while 1:
                    if len(rasteraslist[counter-7]) > int(properties["ncols"]):
                        rasteraslist[counter-7].pop()
                    else:
                        break
                    
                arcpy.SetProgressorPosition(counter-7)
    infile.close()
    
    return rasteraslist, properties


def generate_blank_raster(properties):
    raster = []
    for x in range(int(properties["nrows"])):
        row = [0 for y in range(int(properties["ncols"]))]
        raster.append(row)
    return raster


def python_list_to_ascii(rasteraslist, properties, workspace):

    modasciifile = arcpy.CreateUniqueName("modified.asc", workspace)

    # write out the modified raster to a text file
    infile = open(modasciifile, "w")
    infile.write("ncols         " + properties["ncols"])
    infile.write("nrows         " + properties["nrows"])
    infile.write("xllcorner     " + properties["xllcorner"])
    infile.write("yllcorner     " + properties["yllcorner"])
    infile.write("cellsize      " + properties["cellsize"])
    infile.write("NODATA_value  " + properties["NODATA_value"])

    for x in range(int(properties["nrows"])):
        for x in rasteraslist[x]:
            infile.write(str(x) + " ")
    infile.close()
    
    return modasciifile


def execute(inraster, outraster, wsize, niter, qtype, workspace):
    global window
    global shift
    global noData

    n = len(workspace)
    if n == 0:
        workspace = os.path.dirname(outraster)
        n = len(workspace)
    if n > 4:
        end = workspace[n-4 : n] # extract last 4 letters
        if end == ".gdb": # geodatabase
            workspace = os.path.dirname(workspace)

    shift = [i for i in range(-math.trunc(wsize / 2), math.trunc(wsize / 2) + 1)]

    corrected_size = 2*math.trunc(wsize / 2) + 1
    
    window = [0 for i in range(corrected_size**2)]

    fieldtype = "INTEGER"

    valuetype = int(str(arcpy.GetRasterProperties_management(inraster, "VALUETYPE")))
    if 8 < valuetype < 11:
        fieldtype = "FLOAT"
    elif valuetype >= 11:
        arcpy.AddMessage("Unsupported field format: " + str(valuetype) + "\nSetting to FLOAT")
        fieldtype = "FLOAT"

    # Select the appropriate filter number
    nfilt = 0
    if(qtype == "Lower Quartile"):
        nfilt = 0
    elif(qtype == "Upper Quartile"):
        nfilt = 1
    elif(qtype == "Min"):
        nfilt = 2
    elif(qtype == "Max"):
        nfilt = 3
    elif(qtype == "Mean"):
        nfilt = 4
    elif(qtype == "Median"):
        nfilt = 5
    else:
        nfilt = 4

    # Process: RasterToASCII_conversion
    arcpy.AddMessage("Reading raster...")
    inasciifile = arcpy.CreateUniqueName("temp.asc", workspace)

    arcpy.RasterToASCII_conversion(inraster, inasciifile)

    # Convert ASCII to python list
    rasteraslist, properties = ascii_to_python_list(inasciifile, fieldtype)
    arcpy.AddMessage("Reading completed")

    if fieldtype == "INTEGER":
        noData = int(properties["NODATA_value"])
    else:
        noData = float(properties["NODATA_value"])

    # Calculating quartiles
    arcpy.AddMessage("Filtering raster...")
    newrasterlist = process_raster(rasteraslist, properties, niter, nfilt)

    # Convert python list to ASCII
    arcpy.SetProgressor("default", "Finishing...")
    arcpy.AddMessage("Writing output raster...")
    modasciifile = python_list_to_ascii(newrasterlist, properties, workspace)

    # Write ASCII to output raster
    arcpy.ASCIIToRaster_conversion(modasciifile, outraster)
    arcpy.Delete_management(inasciifile)

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
        
