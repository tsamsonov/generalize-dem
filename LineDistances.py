import os
import sys
import arcpy
import traceback
import math
import numpy
import Utils
from scipy.spatial.distance import cdist

def execute(in_hydrolines, hydro_field, in_counterparts, count_field, deviation, out_table):

    coords_hydro = Utils.get_coordinates(in_hydrolines)
    coords_count = Utils.get_coordinates(in_counterparts)
    id_hydro = Utils.get_values(in_hydrolines, hydro_field)
    id_count = Utils.get_values(in_counterparts, count_field)

    arcpy.CreateTable_management(os.path.dirname(out_table), os.path.basename(out_table))

    arcpy.AddField_management(out_table, 'ID', 'FLOAT', field_length=16)
    arcpy.AddField_management(out_table, 'frechet', 'FLOAT', field_length=16)
    arcpy.AddField_management(out_table, 'hausdorff', 'FLOAT', field_length=16)
    arcpy.AddField_management(out_table, 'hausdorff_forw', 'FLOAT', field_length=16)
    arcpy.AddField_management(out_table, 'hausdorff_back', 'FLOAT', field_length=16)
    arcpy.AddField_management(out_table, 'quality', 'TEXT', field_length=16)

    N = len(id_hydro)
    idx = []
    for i in range(N):
        idx.append(numpy.where(id_count == id_hydro[i])[0])

    coords_count = [coords_count[i] for i in idx]

    insertcursor = arcpy.da.InsertCursor(out_table, ['ID', 'frechet', 'hausdorff', 'hausdorff_forw', 'hausdorff_back', 'quality'])

    for i in range(N):
        id = id_hydro[i]
        arcpy.AddMessage('ID = ' + str(id))

        frechet = Utils.frechet_dist(coords_count[i], coords_hydro[i])
        haus = Utils.hausdorff_dist(coords_count[i], coords_hydro[i])
        haus_forw = Utils.hausdorff_dist_dir(coords_count[i], coords_hydro[i])
        haus_back = Utils.hausdorff_dist_dir(coords_hydro[i], coords_count[i])

        quality = 'Unknown'

        if frechet <= deviation:
            quality = 'Strong'
        elif haus <= deviation:
            quality = 'Regular'
        elif haus_forw <= deviation:
            quality = 'Weak'

        insertcursor.insertRow([id, frechet, haus, haus_forw, haus_back, quality])

    return

if __name__ == 'main':
    in_hydrolines = arcpy.GetParameterAsText(0)
    hydro_field = arcpy.GetParameterAsText(1)
    in_counterparts = arcpy.GetParameterAsText(2)
    count_field = int(arcpy.GetParameterAsText(3))
    deviation = float(arcpy.GetParameterAsText(4))
    out_table = arcpy.GetParameterAsText(5)

    try:
        execute(in_hydrolines, hydro_field, in_counterparts, count_field, deviation, out_table)

    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        print("Processing failed")