# -*- coding: cp1251 -*-
import arcpy
import sys
import traceback
import Utils
import os

# 2020, Timofey Samsonov, Lomonosov Moscow State University


def execute(in_rasters, in_subdivision, in_field, in_overlap, out_raster):

    value_table = arcpy.ValueTable(2)
    value_table.loadFromString(in_rasters)

    ids = sorted(Utils.get_values(in_subdivision, in_field))

    sublyr = 'sublyr'
    arcpy.MakeFeatureLayer_management(in_subdivision, sublyr)

    # wsp = Utils.CreateScratchWorkspace(os.path.dirname(out_raster))

    dems = []

    for i in range(0, value_table.rowCount):
        dem = value_table.getValue(i, 0)
        iscrop = True if value_table.getValue(i, 1) == 'Yes' else False

        arcpy.AddMessage(in_field + ' = "' + ids[i] + '"')

        arcpy.SelectLayerByAttribute_management(sublyr, 'NEW_SELECTION', in_field + " = '" + ids[i] + "'")

        mask = sublyr
        if (not iscrop):
            mask = 'in_memory/mask'
            arcpy.Buffer_analysis(sublyr, mask, 0.5 * in_overlap)

        dems.append(arcpy.sa.ExtractByMask(dem, mask))

    arcpy.MosaicToNewRaster_management(dems,
                                       os.path.dirname(out_raster),
                                       os.path.basename(out_raster),
                                       pixel_type = '32_BIT_FLOAT',
                                       number_of_bands = 1,
                                       mosaic_method = 'BLEND')


    return

if __name__ == "__main__":
    try:
        in_rasters = arcpy.GetParameterAsText(0)
        in_subdivision = arcpy.GetParameterAsText(1)
        in_field = arcpy.GetParameterAsText(2)
        in_overlap = arcpy.GetParameterAsText(3)
        out_raster = arcpy.GetParameterAsText(4)

        execute(in_rasters, in_subdivision, in_field, in_overlap, out_raster)
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)