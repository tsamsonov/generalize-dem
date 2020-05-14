# -*- coding: cp1251 -*-
# DEM generalization algorithm
# 2011-2017 Timofey Samsonov, Lomonosov Moscow State University

import sys
import arcpy
import traceback

import FilterDEM as FD
import CarveDEM as CD
import MosaicDEM as MD
import ConflationLinks as CL
import CreateFishnet as CF
import ConflateDEMbyLinks as CB
import LineDistances as LD
import GeneralizeDEM as GD
import ExtractStreams as ES
import CounterpartStreams as CS
import WidenLandforms as WL

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Generalize DEM"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [CreateFishnet, CarveDEM, CalculateLineDistances, ExtractStreams, CounterpartStreams, GenerateConflationLinks, FilterDEM, MosaicDEM, WidenLandforms, GeneralizeDEM, ConflateDEMbyLinks]

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

class CounterpartStreams(object):

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Counterpart Streams"
        self.description = ""
        self.canRunInBackground = True

    def getParameterInfo(self):

        in_streams = arcpy.Parameter(
            displayName="Input reference hydrographic lines",
            name="in_streams",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        in_field = arcpy.Parameter(
            displayName="Hydrographic line ID field",
            name="in_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        in_field.filter.list = ['Short', 'Long']
        in_field.parameterDependencies = [in_streams.name]

        in_raster = arcpy.Parameter(
            displayName="Input flow accumulation raster",
            name="in_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        dem_raster = arcpy.Parameter(
            displayName="Input elevation raster",
            name="dem_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        out_streams = arcpy.Parameter(
            displayName="Output counterpart streams feature class",
            name="out_streams",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        min_acc = arcpy.Parameter(
            displayName="Minimum flow accumulation",
            name="min_acc",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        min_acc.value = 10

        penalty = arcpy.Parameter(
            displayName="Offstream penalty",
            name="penalty",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        penalty.value = 30

        radius = arcpy.Parameter(
            displayName="Catch radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        deviation = arcpy.Parameter(
            displayName="Maximum deviation",
            name="deviation",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        limit = arcpy.Parameter(
            displayName="Deviation distance metric (flowline counterparts only)",
            name="limit",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        limit.value = 'DIRECTED HAUSDORFF'
        limit.filter.list = ['DIRECTED HAUSDORFF', 'HAUSDORFF', 'FRECHET']

        params = [in_streams, in_field, in_raster, dem_raster, out_streams, min_acc, penalty, radius, deviation, limit]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        instreams = parameters[0].valueAsText
        inidfield= parameters[1].valueAsText
        inraster = parameters[2].valueAsText
        demraster = parameters[3].valueAsText
        outstreams = parameters[4].valueAsText
        minacc = float(parameters[5].valueAsText)
        penalty = int(parameters[6].valueAsText)
        radius = float(parameters[7].valueAsText)
        deviation = float(parameters[8].valueAsText)
        limit = parameters[9].valueAsText

        CS.execute(instreams, inidfield, inraster, demraster, outstreams, minacc, penalty, radius, deviation, limit)

        return

class CarveDEM(object):

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Carve DEM Along Streams"
        self.description = ""
        self.canRunInBackground = True

    def getParameterInfo(self):

        in_raster = arcpy.Parameter(
            displayName="Input raster DEM",
            name="in_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        in_streams = arcpy.Parameter(
            displayName="Input reference hydrographic lines",
            name="in_streams",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        in_field = arcpy.Parameter(
            displayName="Hydrographic line ID field",
            name="in_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        in_field.filter.list = ['Short', 'Long']
        in_field.parameterDependencies = [in_streams.name]

        out_raster = arcpy.Parameter(
            displayName="Output raster DEM",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        params = [in_raster, in_streams, in_field, out_raster]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        in_raster = parameters[0].valueAsText
        in_streams = parameters[1].valueAsText
        in_field = parameters[2].valueAsText
        out_raster = parameters[3].valueAsText

        CD.execute(in_raster, in_streams, in_field, out_raster)

        return

class MosaicDEM(object):

    def __init__(self):
        self.label = "Mosaic DEM"
        self.description = ""
        self.canRunInBackground = True
        self.wsize = 3

    def getParameterInfo(self):
        in_rasters = arcpy.Parameter(
            displayName="Input raster DEMs",
            name="in_rasters",
            datatype="GPValueTable",
            parameterType="Required",
            direction="Input")

        in_rasters.columns = [['GPRasterLayer', 'Raster'], ['String', 'Crop?']]
        in_rasters.filters[1].type = 'ValueList'
        in_rasters.filters[1].list = ['No', 'Yes']

        in_subdivision = arcpy.Parameter(
            displayName="Input mosaic polygons",
            name="in_subdivision",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        in_field = arcpy.Parameter(
            displayName="Order field",
            name="in_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        in_field.parameterDependencies = [in_subdivision.name]

        in_overlap = arcpy.Parameter(
            displayName="Overlap distance",
            name="in_overlap",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        out_raster = arcpy.Parameter(
            displayName="Input raster DEMs",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        params = [in_rasters, in_subdivision, in_field, in_overlap, out_raster]
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
            in_rasters = parameters[0].valueAsText
            in_subdivision = parameters[1].valueAsText
            in_field = parameters[2].valueAsText
            in_overlap = float(parameters[3].valueAsText)
            out_raster = parameters[4].valueAsText

            MD.execute(in_rasters, in_subdivision, in_field, in_overlap, out_raster)
        except:
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                    str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
            arcpy.AddError(pymsg)
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

        ridges = arcpy.Parameter(
            displayName="Widen ridges?",
            name="ridges",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [inraster, streams, outraster, distance, windowsize, ftype, ridges]
        return params

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed

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
        ridges = parameters[6].valueAsText

        WL.execute(demdataset, streams, distance, windowsize, output, ftype, ridges)

        return

class ConflateDEMbyLinks(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Conflate DEM by Links"
        self.description = ""
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""
        in_raster = arcpy.Parameter(
            displayName="Input raster DEM",
            name="in_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        in_links = arcpy.Parameter(
            displayName="Input conflation links",
            name="in_links",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        in_area = arcpy.Parameter(
            displayName="Input conflation area",
            name="in_area",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        distance = arcpy.Parameter(
            displayName="Conflation distance",
            name="distance",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        out_raster = arcpy.Parameter(
            displayName="Output raster DEM",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        params = [in_raster, in_links, in_area, distance, out_raster]
        return params

    def isLicensed(self):

        return True  # tool can be executed

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        in_raster = parameters[0].valueAsText
        in_links = parameters[1].valueAsText
        in_area = parameters[2].valueAsText
        distance = float(parameters[3].valueAsText)
        out_raster = parameters[4].valueAsText

        CB.execute(in_raster, in_links, in_area, distance, out_raster)

        return

class GenerateConflationLinks(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Generate Conflation Links"
        self.description = ""
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""

        in_hydrolines = arcpy.Parameter(
            displayName="Input reference hydrographic lines",
            name="in_hydrolines",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        hydro_field = arcpy.Parameter(
            displayName="Hydrographic line ID field",
            name="hydro_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        hydro_field.filter.list = ['Short', 'Long']
        hydro_field.parameterDependencies = [in_hydrolines.name]

        in_counterparts = arcpy.Parameter(
            displayName="Input counterpart streams",
            name="in_counterparts",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        count_field = arcpy.Parameter(
            displayName="Couinterpart line ID field",
            name="count_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        count_field.filter.list = ['Short', 'Long']
        count_field.parameterDependencies = [in_counterparts.name]

        out_links = arcpy.Parameter(
            displayName="Output conflation links",
            name="out_links",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        out_area = arcpy.Parameter(
            displayName="Output conflation area",
            name="out_area",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Output")

        params = [in_hydrolines, hydro_field, in_counterparts, count_field, out_links, out_area]
        return params

    def isLicensed(self):

        return True  # tool can be executed

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        in_hydrolines = parameters[0].valueAsText
        hydro_field = parameters[1].valueAsText
        in_counterparts = parameters[2].valueAsText
        count_field = parameters[3].valueAsText
        out_links = parameters[4].valueAsText
        out_area = parameters[5].valueAsText

        CL.execute(in_hydrolines, hydro_field, in_counterparts, count_field, out_links, out_area)

        return

class CalculateLineDistances(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Calculate Distances Between Lines"
        self.description = "The tool calculates Directed Hausdorff, Hausdorff and Frechet distances between corresponding lines in two layers"
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""

        in_hydrolines = arcpy.Parameter(
            displayName="Input reference hydrographic lines",
            name="in_hydrolines",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        hydro_field = arcpy.Parameter(
            displayName="Hydrographic line ID field",
            name="hydro_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        hydro_field.filter.list = ['Short', 'Long']
        hydro_field.parameterDependencies = [in_hydrolines.name]

        in_counterparts = arcpy.Parameter(
            displayName="Input counterpart streams",
            name="in_counterparts",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        count_field = arcpy.Parameter(
            displayName="Couinterpart line ID field",
            name="count_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        count_field.filter.list = ['Short', 'Long']
        count_field.parameterDependencies = [in_counterparts.name]

        deviation = arcpy.Parameter(
            displayName="Maximum deviation",
            name="deviation",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        out_table = arcpy.Parameter(
            displayName="Output table",
            name="out_table",
            datatype="DETable",
            parameterType="Required",
            direction="Output")

        params = [in_hydrolines, hydro_field, in_counterparts, count_field, deviation, out_table]
        return params

    def isLicensed(self):

        return True  # tool can be executed

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        in_hydrolines = parameters[0].valueAsText
        hydro_field = parameters[1].valueAsText
        in_counterparts = parameters[2].valueAsText
        count_field = parameters[3].valueAsText
        deviation = float(parameters[4].valueAsText)
        out_table = parameters[5].valueAsText

        LD.execute(in_hydrolines, hydro_field, in_counterparts, count_field, deviation, out_table)

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
        demdataset.category = '1. Input and output DEM'

        output = arcpy.Parameter(
            displayName="Output raster DEM",
            name="output",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")
        output.category = '1. Input and output DEM'

        flowdir = arcpy.Parameter(
            displayName="Input flow direction raster",
            name="flowdir",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")
        flowdir.category = '2. Supplementary input'

        flowacc = arcpy.Parameter(
            displayName="Input flow accumulation raster",
            name="flowacc",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")
        flowacc.category = '2. Supplementary input'

        flines = arcpy.Parameter(
            displayName="Input linear features",
            name="flines",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        flines.category = '2. Supplementary input'

        fpolys = arcpy.Parameter(
            displayName="Input polygonal features",
            name="fpolys",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        fpolys.category = '2. Supplementary input'

        cliparea = arcpy.Parameter(
            displayName="Clip area",
            name="cliparea",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        cliparea.category = '2. Supplementary input'

        outputcellsize = arcpy.Parameter(
            displayName="Output cell size",
            name="outputcellsize",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        outputcellsize.category = '3. Main parameters'
        outputcellsize.value = 1000

        minacc1 = arcpy.Parameter(
            displayName="Minimum primary flow accumulation",
            name="minacc1",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minacc1.category = '3. Main parameters'
        minacc1.value = 40

        minlen1 = arcpy.Parameter(
            displayName="Minimum primary flow length",
            name="minlen1",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minlen1.category = '3. Main parameters'
        minlen1.value = 40

        minacc2 = arcpy.Parameter(
            displayName="Minimum secondary flow accumulation",
            name="minacc2",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minacc2.category = '3. Main parameters'
        minacc2.value = 20

        minlen2 = arcpy.Parameter(
            displayName="Minimum secondary flow length",
            name="minlen2",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minlen2.category = '3. Main parameters'
        minlen2.value = 10

        is_widen = arcpy.Parameter(
            displayName="Widen Landforms",
            name="is_widen",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_widen.category = '4. Widening and smoothing'
        is_widen.value = 'true'

        widentype = arcpy.Parameter(
            displayName="Widening method",
            name="widentype",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        widentype.category = '4. Widening and smoothing'
        widentype.filter.list = ['Min/Max', 'Lower/Upper Quartile']
        widentype.value = 'Min/Max'

        widendist = arcpy.Parameter(
            displayName="Widening distance",
            name="widendist",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        widendist.category = '4. Widening and smoothing'
        widendist.value = 8000

        filtersize = arcpy.Parameter(
            displayName="Widening filter size",
            name="filtersize",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        filtersize.category = '4. Widening and smoothing'
        filtersize.value = 3

        is_smooth = arcpy.Parameter(
            displayName="Smooth DEM",
            name="is_smooth",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_smooth.category = '4. Widening and smoothing'
        is_smooth.value = 'false'

        is_tiled = arcpy.Parameter(
            displayName="Tile the source DEM",
            name="is_tiled",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_tiled.category = '5. Tiling and parallel processing'
        is_tiled.value = 'false'

        tile_size = arcpy.Parameter(
            displayName="Tile size",
            name="tile_size",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        tile_size.category = '5. Tiling and parallel processing'
        tile_size.value = 2048

        is_parallel = arcpy.Parameter(
            displayName="Parallel processing",
            name="is_parallel",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_parallel.category = '5. Tiling and parallel processing'
        is_parallel.value = 'false'

        num_processes = arcpy.Parameter(
            displayName="Number of processes",
            name="num_processes",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        num_processes.category = '5. Tiling and parallel processing'
        num_processes.value = 0

        is_continued = arcpy.Parameter(
            displayName="Continue previous processing",
            name="is_continued",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_continued.category = '6. Continue previous processing'
        is_continued.value = 'false'

        continued_folder = arcpy.Parameter(
            displayName="Scratch folder from previous processing",
            name="continued_folder",
            datatype="DEFolder",
            parameterType="Optional",
            direction="Input")
        continued_folder.category = '5. Continue previous processing'

        params = [demdataset, output, flowdir, flowacc, flines, fpolys, cliparea,
                  outputcellsize, minacc1, minlen1, minacc2, minlen2,
                  is_widen, widentype, widendist, filtersize, is_smooth, is_tiled, tile_size,
                  is_parallel, num_processes, is_continued, continued_folder]
        return params

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available" or arcpy.CheckExtension("3D") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed
        return True  # tool can be executed

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        demdataset = parameters[0].valueAsText
        output = parameters[1].valueAsText
        flowdir = parameters[2].valueAsText
        flowacc = parameters[3].valueAsText
        flines = parameters[4].valueAsText
        fpolys = parameters[5].valueAsText
        cliparea = parameters[6].valueAsText
        outputcellsize = float(parameters[7].valueAsText)
        minacc1 = int(parameters[8].valueAsText)
        minlen1 = int(parameters[9].valueAsText)
        minacc2 = int(parameters[10].valueAsText)
        minlen2 = int(parameters[11].valueAsText)
        is_widen = True if parameters[12].valueAsText == 'true' else False
        widentype = parameters[13].valueAsText
        widendist = float(parameters[14].valueAsText)
        filtersize = int(parameters[15].valueAsText)
        is_smooth = True if parameters[16].valueAsText == 'true' else False
        is_tiled = True if parameters[17].valueAsText == 'true' else False
        tile_size = int(parameters[18].valueAsText)
        is_parallel = True if parameters[19].valueAsText == 'true' else False
        num_processes = float(parameters[20].valueAsText)
        is_continued = True if parameters[21].valueAsText == 'true' else False
        continued_folder = parameters[22].valueAsText

        GD.execute(demdataset,
                   output,
                   flowdir,
                   flowacc,
                   flines,
                   fpolys,
                   cliparea,
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
                   is_tiled,
                   tile_size,
                   is_parallel,
                   num_processes,
                   is_continued,
                   continued_folder)

        return