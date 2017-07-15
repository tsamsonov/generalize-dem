import os
import sys
import arcpy
import traceback

def execute(template, output, nrows, ncols, overlap, split = False, overlap2 = 0):

    desc = arcpy.Describe(template)

    adjustment = 0
    # if desc.dataType == "RasterLayer" or desc.dataType == "RasterDataset":
    #     adjustment = desc.meanCellWidth

    xmin = desc.extent.XMin
    xmax = desc.extent.XMax - adjustment
    ymin = desc.extent.YMin + adjustment
    ymax = desc.extent.YMax

    Lx = xmax - xmin
    Ly = ymax - ymin

    wx = (Lx + (ncols - 1) * overlap) / ncols
    wy = (Ly + (nrows - 1) * overlap) / nrows

    xcoords = []
    ycoords = []

    jdx = range(nrows)
    idx = range(ncols)

    if split:
        wx0 = wx - 0.5 * overlap  # width of the first and last cell
        wx1 = wx - overlap  # width of other cells

        wy0 = wy - 0.5 * overlap  # height of the first and last cell
        wy1 = wy - overlap  # height of other cells

        x = xmin
        xcoords.append(xmin)
        for i in range(ncols):
            if i == 0 or i == ncols - 1:
                x += wx0
                xcoords.append(x)
            else:
                x += wx1
                xcoords.append(x)

        y = ymin
        ycoords.append(ymin)
        for j in range(nrows):
            if j == 0 or j == nrows - 1:
                y += wy0
                ycoords.append(y)
            else:
                y += wy1
                ycoords.append(y)
    else:
        if overlap2 == 0:
            x = xmin
            for i in range(ncols):
                xcoords.append(x)
                x += wx
                xcoords.append(x)
                x -= overlap

            y = ymin
            for j in range(nrows):
                ycoords.append(y)
                y += wy
                ycoords.append(y)
                y -= overlap
        else:
            wx0 = wx - 0.5 * overlap + 0.5 * overlap2
            wx1 = wx - overlap + overlap2
            x = xmin
            for i in range(ncols):
                xcoords.append(x)
                if i == 0 or i == ncols - 1:
                    x += wx0
                else:
                    x += wx1
                xcoords.append(x)
                x -= overlap2

            wy0 = wy - 0.5 * overlap + 0.5 * overlap2
            wy1 = wy - overlap + overlap2
            y = ymin
            for j in range(nrows):
                ycoords.append(y)
                if j == 0 or j == nrows - 1:
                    y += wy0
                else:
                    y += wy1
                ycoords.append(y)
                y -= overlap2

        jdx = range(0, 2 * nrows, 2)
        idx = range(0, 2 * ncols, 2)


    workspace = os.path.dirname(output)
    name = os.path.basename(output)

    arcpy.CreateFeatureclass_management(workspace, name, "POLYGON", spatial_reference=template)

    fishnet = workspace + '/' + name

    cursor = arcpy.da.InsertCursor(fishnet, ["SHAPE@"])

    for j in jdx:
        for i in idx:
            points = [arcpy.Point(xcoords[i], ycoords[j]),
                      arcpy.Point(xcoords[i + 1], ycoords[j]),
                      arcpy.Point(xcoords[i + 1], ycoords[j + 1]),
                      arcpy.Point(xcoords[i], ycoords[j + 1])]
            array = arcpy.Array(points)
            polygon = arcpy.Polygon(array)
            cursor.insertRow([polygon])

    del cursor

    return

if __name__ == 'main':
    template = arcpy.GetParameterAsText(0)
    output = arcpy.GetParameterAsText(1)
    nrows = int(arcpy.GetParameterAsText(2))
    ncols = int(arcpy.GetParameterAsText(3))
    overlap = float(arcpy.GetParameterAsText(4))
    overlap2 = float(arcpy.GetParameterAsText(5))
    split = arcpy.GetParameterAsText(6)

    try:
        execute(template, output, nrows, ncols, overlap, split, overlap2)

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        print("Processing failed")