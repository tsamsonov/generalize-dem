# -*- coding: cp1251 -*-
# DEM generalization algorithm
# 2011-2017 Timofey Samsonov, Lomonosov Moscow State University

import sys
import arcpy
import traceback

import FilterDEM as FD
import CreateFishnet as CF
import GeneralizeDEM as GD
import ExtractStreams as ES
import WidenLandforms as WL

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Generalize DEM"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [CreateFishnet, ExtractStreams, FilterDEM, WidenLandforms, GeneralizeDEM]

class CreateFishnet(object):

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Fishnet"
        self.description = ""
        self.canRunInBackground = True

    def getParameterInfo(self):

        template = arcpy.Parameter(
            displayName="Template layer",
            name="template",
            datatype="GPLayer",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output fishnet polygon feature class",
            name="output",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        nrows = arcpy.Parameter(
            displayName="Number of columns",
            name="nrows",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        nrows.value = 2

        ncols = arcpy.Parameter(
            displayName="Number of rows",
            name="ncols",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        ncols.value = 2

        overlap = arcpy.Parameter(
            displayName="Overlap between tiles (in template units)",
            name="overlap",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        overlap.value = 0

        overlap2 = arcpy.Parameter(
            displayName="Second-order overlap between tiles (in template units)",
            name="overlap2",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        overlap2.value = 0

        split = arcpy.Parameter(
            displayName="Split tiles at the middle of the overlap",
            name="split",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [template, output, nrows, ncols, overlap, overlap2, split]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        template = parameters[0].valueAsText
        output = parameters[1].valueAsText
        nrows = int(parameters[2].valueAsText)
        ncols = int(parameters[3].valueAsText)
        overlap = float(parameters[4].valueAsText)
        overlap2 = float(parameters[5].valueAsText)
        split = parameters[6].valueAsText


        CF.execute(template, output, nrows, ncols, overlap, split, overlap2)

        return

class ExtractStreams(object):

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Streams"
        self.description = ""
        self.canRunInBackground = True

    def getParameterInfo(self):

        in_raster = arcpy.Parameter(
            displayName="Input flow accumulation raster",
            name="in_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        out_raster = arcpy.Parameter(
            displayName="Output stream raster",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        min_acc = arcpy.Parameter(
            displayName="Minimum flow accumulation",
            name="min_acc",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        min_len = arcpy.Parameter(
            displayName="Minimum stream length (in cells)",
            name="min_len",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        params = [in_raster, out_raster, min_acc, min_len]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        inraster = parameters[0].valueAsText
        outraster = parameters[1].valueAsText
        minacc = float(parameters[2].valueAsText)
        minlen = int(parameters[3].valueAsText)

        ES.execute(inraster, outraster, minacc, minlen)

        return


class FilterDEM(object):

    def __init__(self):
        self.label = "Filter DEM"
        self.description = ""
        self.canRunInBackground = True
        self.wsize = 3

    def getParameterInfo(self):
        in_raster = arcpy.Parameter(
            displayName="Input raster DEM",
            name="in_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        out_raster = arcpy.Parameter(
            displayName="Output raster DEM",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        wsize = arcpy.Parameter(
            displayName="Filter size (in cells)",
            name="min_len",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        wsize.value = 3

        niter = arcpy.Parameter(
            displayName="Number of iterations",
            name="min_acc",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        niter.value = 1

        ftype = arcpy.Parameter(
            displayName="Statistics type",
            name="ftype",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        ftype.value = 'Min'
        ftype.filter.list = ['Min', 'Max', 'Mean', 'Median', 'Upper Quartile', 'Lower Quartile']

        params = [in_raster, out_raster, wsize, niter, ftype]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        try:
            # Get input parameters
            in_raster = parameters[0].valueAsText
            out_raster = parameters[1].valueAsText
            wsize = int(parameters[2].valueAsText)
            niter = int(parameters[3].valueAsText)
            ftype = parameters[4].valueAsText

            FD.execute(in_raster, out_raster, wsize, niter, ftype)
        except:
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                    str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
            arcpy.AddError(pymsg)
        return


class WidenLandforms(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Widen Landforms"
        self.description = ""
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""
        inraster = arcpy.Parameter(
            displayName="Input raster DEM",
            name="inraster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input streams feature layer",
            name="streams",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        outraster = arcpy.Parameter(
            displayName="Output raster DEM",
            name="outraster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        distance = arcpy.Parameter(
            displayName="Widening distance (in DEM projection units)",
            name="distance",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        distance.value = 1000

        windowsize = arcpy.Parameter(
            displayName="Filter size (in cells)",
            name="windowsize",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        windowsize.value = 3

        ftype = arcpy.Parameter(
            displayName="Widening statistics",
            name="ftype",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        ftype.value = 'Min/Max'
        ftype.filter.list = ['Min/Max', 'Lower/Upper Quartile']

        params = [inraster, streams, outraster, distance, windowsize, ftype]
        return params

    def isLicensed(self):
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False  # tool cannot be executed

        return True  # tool can be executed

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        demdataset = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        distance = float(parameters[3].valueAsText)
        windowsize = int(parameters[4].valueAsText)
        ftype = parameters[5].valueAsText

        WL.execute(demdataset, streams, distance, windowsize, output, ftype)

        return


class GeneralizeDEM(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Generalize DEM"
        self.description = ""
        self.canRunInBackground = True

    def getParameterInfo(self):

        demdataset = arcpy.Parameter(
            displayName="Input raster DEM",
            name="demdataset",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        demdataset.category = '1. Input and output'

        marine = arcpy.Parameter(
            displayName="Marine area polygon feature layer",
            name="marine",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        marine.category = '1. Input and output'

        output = arcpy.Parameter(
            displayName="Output raster DEM",
            name="output",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")
        output.category = '1. Input and output'

        outputcellsize = arcpy.Parameter(
            displayName="Output cell size",
            name="outputcellsize",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        outputcellsize.category = '2. Main parameters'
        outputcellsize.value = 1000

        minacc1 = arcpy.Parameter(
            displayName="Minimum primary flow accumulation",
            name="minacc1",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minacc1.category = '2. Main parameters'
        minacc1.value = 40

        minlen1 = arcpy.Parameter(
            displayName="Minimum primary flow length",
            name="minlen1",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minlen1.category = '2. Main parameters'
        minlen1.value = 40

        minacc2 = arcpy.Parameter(
            displayName="Minimum secondary flow accumulation",
            name="minacc2",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minacc2.category = '2. Main parameters'
        minacc2.value = 20

        minlen2 = arcpy.Parameter(
            displayName="Minimum secondary flow length",
            name="minlen2",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minlen2.category = '2. Main parameters'
        minlen2.value = 10

        is_widen = arcpy.Parameter(
            displayName="Widen Landforms",
            name="is_widen",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_widen.category = '3. Widening and smoothing'
        is_widen.value = 'true'

        widentype = arcpy.Parameter(
            displayName="Widening method",
            name="widentype",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        widentype.category = '3. Widening and smoothing'
        widentype.filter.list = ['Min/Max', 'Lower/Upper Quartile']
        widentype.value = 'Min/Max'

        widendist = arcpy.Parameter(
            displayName="Widening distance",
            name="widendist",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        widendist.category = '3. Widening and smoothing'
        widendist.value = 8000

        filtersize = arcpy.Parameter(
            displayName="Widening filter size",
            name="filtersize",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        filtersize.category = '3. Widening and smoothing'
        filtersize.value = 3

        is_smooth = arcpy.Parameter(
            displayName="Smooth DEM",
            name="is_smooth",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_smooth.category = '3. Widening and smoothing'
        is_smooth.value = 'false'

        is_tiled = arcpy.Parameter(
            displayName="Tile the source DEM",
            name="is_tiled",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_tiled.category = '4. Tiling and parallel processing'
        is_tiled.value = 'true'

        tile_size = arcpy.Parameter(
            displayName="Tile size",
            name="tile_size",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        tile_size.category = '4. Tiling and parallel processing'
        tile_size.value = 2048

        is_parallel = arcpy.Parameter(
            displayName="Parallel processing",
            name="is_parallel",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_parallel.category = '4. Tiling and parallel processing'
        is_parallel.value = 'false'

        num_processes = arcpy.Parameter(
            displayName="Number of processes",
            name="num_processes",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        num_processes.category = '4. Tiling and parallel processing'
        num_processes.value = 0

        params = [demdataset, marine, output, outputcellsize, minacc1, minlen1, minacc2, minlen2,
                  is_widen, widentype, widendist, filtersize, is_smooth, is_tiled, tile_size, is_parallel, num_processes]
        return params

    def isLicensed(self):
        try:
            if arcpy.CheckExtension("Spatial") != "Available" or arcpy.CheckExtension("3D") != "Available":
                raise Exception
        except Exception:
            return False  # tool cannot be executed
        return True  # tool can be executed

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        demdataset = parameters[0].valueAsText
        marine = parameters[1].valueAsText
        output = parameters[2].valueAsText
        outputcellsize = float(parameters[3].valueAsText)
        minacc1 = int(parameters[4].valueAsText)
        minlen1 = int(parameters[5].valueAsText)
        minacc2 = int(parameters[6].valueAsText)
        minlen2 = int(parameters[7].valueAsText)
        is_widen = parameters[8].valueAsText
        widentype = parameters[9].valueAsText
        widendist = float(parameters[10].valueAsText)
        filtersize = int(parameters[11].valueAsText)
        is_smooth = parameters[12].valueAsText
        is_tiled = parameters[13].valueAsText
        tile_size = int(parameters[14].valueAsText)
        is_parallel = parameters[15].valueAsText
        num_processes = float(parameters[16].valueAsText)

        GD.execute(demdataset,
                              marine,
                              output,
                              outputcellsize,
                              minacc1,
                              minlen1,
                              minacc2,
                              minlen2,
                              is_widen,
                              widentype,
                              widendist,
                              filtersize,
                              is_smooth,
                              tile_size,
                              num_processes,
                              is_parallel)

        return