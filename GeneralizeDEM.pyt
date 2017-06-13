# -*- coding: cp1251 -*-
# DEM generalization algorithm
# 2011-2017 Timofey Samsonov, Lomonosov Moscow State University

import arcpy
import numpy
import sys
import math
import traceback
import os.path
import multiprocessing
from arcpy.sa import *
from Test import execute

class Tools:
    @staticmethod
    def extend_array(array, nx, ny, value):
        ni = array.shape[0]
        nj = array.shape[1]

        extarray = numpy.empty((ni + ny, nj + nx))
        extarray.fill(value)
        for i in range(ni):
            for j in range(nj):
                extarray[i, j] = array[i, j]
        return extarray


class Filter:

    def __init__(self, wsize, nodata):

        corrected_size = 2 * math.trunc(wsize / 2) + 1

        self.shift = [i for i in range(-math.trunc(wsize / 2), math.trunc(wsize / 2) + 1)]
        self.window = [0 for i in range(corrected_size ** 2)]
        self.nodata = nodata

        self.filters = {'Min': self.calc_min,
                        'Max': self.calc_max,
                        'Mean': self.calc_mean,
                        'Median': self.calc_median,
                        'Lower Quartile': self.calc_lower_quartile,
                        'Upper Quartile': self.calc_upper_quartile
        }

    def sample_window(self, raster, i, j):
        w = 0
        for k in self.shift:
            ik = i + k
            for l in self.shift:
                jl = j + l
                if raster[ik, jl] != self.nodata:
                    self.window[w] = raster[ik, jl]
                    w += 1
        elems = self.window[0:w]
        return elems

    # nfilt = 0
    def calc_lower_quartile(self, raster, i, j):
        elems = self.sample_window(raster, i, j)
        elems.sort()
        n = int(math.floor(len(elems) / 4))
        return 0.5 * (elems[n] + elems[n - 1])

    # nfilt = 1
    def calc_upper_quartile(self, raster, i, j):
        elems = self.sample_window(raster, i, j)
        elems.sort(reverse=True)
        n = int(math.floor(len(elems) / 4))
        return 0.5 * (elems[n] + elems[n - 1])

    # nfilt = 2
    def calc_min(self, raster, i, j):
        elems = self.sample_window(raster, i, j)
        value = min(elems)
        return value

    # nfilt = 3
    def calc_max(self, raster, i, j):
        elems = self.sample_window(raster, i, j)
        value = max(elems)
        return value

    # nfilt = 4
    def calc_mean(self, raster, i, j):
        elems = self.sample_window(raster, i, j)
        sum = 0
        for k in range(len(elems)):
            sum += elems[k]
        value = sum / len(elems)
        return value

    # nfilt = 5
    def calc_median(self, raster, i, j):
        elems = self.sample_window(raster, i, j)
        elems.sort()
        n = len(elems)
        if n == 0:
            return self.nodata
        elif n % 2 == 0:
            return elems[n // 2 - 1]
        else:
            return (elems[n // 2 - 1] + elems[n // 2]) * 0.5

    def filter(self, raster, i, j, ftype):
        return self.filters[ftype](raster, i, j)

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Generalize DEM"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [ExtractStreams, FilterDEM, WidenLandforms, GeneralizeDEM]


class ExtractStreams(object):

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Streams"
        self.description = ""
        self.canRunInBackground = True
        self.MAXACC = 0

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

    def find_up_cell(self, accraster, i, j):
        w = [[0.70710678, 1, 0.70710678], [1, 1, 1], [0.70710678, 1, 0.70710678]]  # distance weights
        shift = [-1, 0, 1]
        minmax = 4000000000
        a = 0
        kmin = 1
        lmin = 1

        # finding differences in 3x3 neighbourhood

        for k in shift:
            b = 0
            ik = i + k
            for l in shift:
                jl = j + l
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

    def trace_flow_cells(self, accraster, streamraster, i, j, stream, minacc, minlen):
        acc = accraster[i, j]
        ik = i
        jk = j
        n = 0

        while n < minlen:
            stream[n] = [ik, jk]
            iup, jup = self.find_up_cell(accraster, ik, jk)
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
                iup, jup = self.find_up_cell(accraster, ik, jk)
                if iup == ik and jup == jk:
                    break
                acc = accraster[iup, jup]
                ik = iup
                jk = jup
                n += 1

        return streamraster

    def process_raster(self, inraster, minacc, minlen):

        arcpy.AddMessage("Streaming...")
        stream = [[0, 0] for i in range(minlen)]
        ni = inraster.shape[0]
        nj = inraster.shape[1]

        arcpy.AddMessage("Zeroing...")
        outraster = numpy.zeros((ni, nj))

        arcpy.AddMessage("Extending...")
        extinraster = Tools.extend_array(inraster, 1, 1, self.MAXACC * 10)

        arcpy.AddMessage("Tracing...")

        arcpy.SetProgressor("step", "Processing rows", 0, ni - 1, 1)
        for i in range(ni):
            arcpy.SetProgressorLabel("Processing row " + str(i) + " from " + str(ni))
            for j in range(nj):
                if inraster[i, j] > minacc:
                    self.trace_flow_cells(extinraster, outraster, i, j, stream, minacc, minlen)
            arcpy.SetProgressorPosition(i)

        return outraster

    def call(self, inraster, outraster, minacc, minlen):

        self.MAXACC = float(str(arcpy.GetRasterProperties_management(inraster, "MAXIMUM")))

        rasternumpy = arcpy.RasterToNumPyArray(inraster)

        # Tracing stream lines
        arcpy.AddMessage("Tracing stream lines...")
        newrasternumpy = self.process_raster(rasternumpy, minacc, minlen)

        desc = arcpy.Describe(inraster)
        lowerleft = arcpy.Point(desc.extent.XMin, desc.extent.YMin)
        cellsize = desc.meanCellWidth
        crs = desc.spatialReference

        arcpy.AddMessage("Writing streams...")
        outinnerraster = arcpy.NumPyArrayToRaster(newrasternumpy, lowerleft, cellsize)
        arcpy.DefineProjection_management(outinnerraster, crs)
        outinnerraster.save(outraster)

    def execute(self, parameters, messages):

        inraster = parameters[0].valueAsText
        outraster = parameters[1].valueAsText
        minacc = float(parameters[2].valueAsText)
        minlen = int(parameters[3].valueAsText)

        self.call(inraster, outraster, minacc, minlen)

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

    def process_raster(self, inraster, niter, ftype, wsize, nodata):

        ni = inraster.shape[0]
        nj = inraster.shape[1]

        arcpy.SetProgressor("step", "Processing rows", 0, ni - 1, 1)
        raster = inraster

        outraster = numpy.zeros((ni, nj))

        # Filter selector
        flt = Filter(wsize, nodata)

        for k in range(niter):
            raster = Tools.extend_array(raster, wsize - 1, wsize - 1, nodata)

            arcpy.AddMessage("Iteration " + str(k + 1))
            for i in range(ni):
                arcpy.SetProgressorLabel("Processing row " + str(i) + " from " + str(ni))
                for j in range(nj):
                    if raster[i, j] == nodata:
                        outraster[i, j] = raster[i, j]
                    else:
                        outraster[i, j] = flt.filter(raster, i, j, ftype)
                arcpy.SetProgressorPosition(i)
            raster = outraster
        return outraster

    def call(self, inraster, outraster, wsize, niter, ftype):

        self.wsize = wsize

        r = arcpy.Raster(inraster)
        lowerleft = arcpy.Point(r.extent.XMin, r.extent.YMin)
        cellsize = r.meanCellWidth
        crs = r.spatialReference

        arcpy.AddMessage("Filtering raster...")
        rasternumpy = arcpy.RasterToNumPyArray(inraster)
        newrasternumpy = self.process_raster(rasternumpy, niter, ftype, wsize, r.noDataValue)

        arcpy.AddMessage("Writing output...")
        outinnerraster = arcpy.NumPyArrayToRaster(newrasternumpy, lowerleft, cellsize)
        arcpy.DefineProjection_management(outinnerraster, crs)
        outinnerraster.save(outraster)

    def execute(self, parameters, messages):
        try:
            # Get input parameters
            in_raster = parameters[0].valueAsText
            out_raster = parameters[1].valueAsText
            wsize = int(parameters[2].valueAsText)
            niter = int(parameters[3].valueAsText)
            ftype = parameters[4].valueAsText

            self.call(in_raster, out_raster, wsize, niter, ftype)
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

    def call(self, demdataset, streams, distance, windowsize, output, ftype):
        arcpy.env.workspace = "in_memory"
        arcpy.env.snapRaster = demdataset
        dem = arcpy.Raster(demdataset)
        cellsize = dem.meanCellHeight

        # Calculate distances from streams
        arcpy.AddMessage("Calculating distances...")
        distances = EucDistance(streams, "", cellsize, "")

        # Derive valleys weights
        arcpy.AddMessage("Calculating valley weights...")
        divdist = Divide(distances, distance)
        divdistminus = Minus(1, divdist)
        w_valleys = Con(divdistminus, 0, divdistminus, "value < 0")

        # Derive ridges weights
        arcpy.AddMessage("Calculating ridge weights...")
        distminus = Minus(distances, distance)
        distminusdiv = Divide(distminus, distance)
        udistminusdiv = Con(distminusdiv, 0, distminusdiv, "value < 0")
        w_ridges = Con(udistminusdiv, 1, udistminusdiv, "value > 1")

        # Derive source dem weights
        arcpy.AddMessage("Calculating source dem weights...")
        w_both = Plus(w_valleys, w_ridges)
        w_dem = Minus(1, w_both)

        # Get valley and ridge values
        arcpy.AddMessage("Filtering elevations...")
        if ftype == "Min/Max":
            neighborhood = NbrRectangle(windowsize, windowsize, "CELL")
            valleys = FocalStatistics(dem, neighborhood, "MINIMUM", "DATA")
            ridges = FocalStatistics(dem, neighborhood, "MAXIMUM", "DATA")
        else:
            val = arcpy.env.workspace + "val"
            rig = arcpy.env.workspace + "rig"
            FilterDEM.call(demdataset, val, windowsize, 1, "Lower Quartile")
            FilterDEM.call(demdataset, rig, windowsize, 1, "Upper Quartile")
            valleys = arcpy.Raster(val)
            ridges = arcpy.Raster(rig)

        # Calculate weighted values
        arcpy.AddMessage("Calculating weighted values...")
        weightedvalleys = Times(valleys, w_valleys)
        weightedridges = Times(ridges, w_ridges)
        weighteddem = Times(dem, w_dem)

        # Mix values
        arcpy.AddMessage("Mixing values...")
        mixvalleyridges = Plus(weightedvalleys, weightedridges)
        mixall = Plus(weighteddem, mixvalleyridges)

        # Save the result
        mixall.save(output)

    def execute(self, parameters, messages):

        demdataset = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        distance = float(parameters[3].valueAsText)
        windowsize = int(parameters[4].valueAsText)
        ftype = parameters[5].valueAsText

        self.call(demdataset, streams, distance, windowsize, output, ftype)

        return


class GeneralizeDEM(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Generalize DEM"
        self.description = ""
        self.canRunInBackground = True
        self.demdataset = None
        self.rastertinworkspace = None
        self.fishbuffer = None
        self.marine = None
        self.workspace = None
        self.cell = None
        self.minacc1 = None
        self.minlen1 = None
        self.minacc2 = None
        self.minlen2 = None
        self.is_widen = None
        self.is_smooth = None
        self.widendist = None
        self.filtersize = None
        self.widentype = None
        self.i = None
        self.N = None

    def getParameterInfo(self):

        demdataset = arcpy.Parameter(
            displayName="Input raster DEM",
            name="demdataset",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        demdataset.category = 'Input and output'

        marine = arcpy.Parameter(
            displayName="Marine area polygon feature layer",
            name="marine",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        marine.category = 'Input and output'

        output = arcpy.Parameter(
            displayName="Output raster DEM",
            name="output",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")
        output.category = 'Input and output'

        outputcellsize = arcpy.Parameter(
            displayName="Output cell size",
            name="outputcellsize",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        outputcellsize.category = 'Main parameters'
        outputcellsize.value = 1000

        minacc1 = arcpy.Parameter(
            displayName="Minimum primary flow accumulation",
            name="minacc1",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minacc1.category = 'Main parameters'
        minacc1.value = 40

        minlen1 = arcpy.Parameter(
            displayName="Minimum primary flow length",
            name="minlen1",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minlen1.category = 'Main parameters'
        minlen1.value = 40

        minacc2 = arcpy.Parameter(
            displayName="Minimum secondary flow accumulation",
            name="minacc2",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minacc2.category = 'Main parameters'
        minacc2.value = 20

        minlen2 = arcpy.Parameter(
            displayName="Minimum secondary flow length",
            name="minlen2",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        minlen2.category = 'Main parameters'
        minlen2.value = 10

        is_widen = arcpy.Parameter(
            displayName="Widen Landforms",
            name="is_widen",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_widen.category = 'Widening and smoothing'
        is_widen.value = 'true'

        widentype = arcpy.Parameter(
            displayName="Widening method",
            name="widentype",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        widentype.category = 'Widening and smoothing'
        widentype.filter.list = ['Min/Max', 'Lower/Upper Quartile']
        widentype.value = 'Min/Max'

        widendist = arcpy.Parameter(
            displayName="Widening distance",
            name="widendist",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        widendist.category = 'Widening and smoothing'
        widendist.value = 8000

        filtersize = arcpy.Parameter(
            displayName="Widening filter size",
            name="filtersize",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        filtersize.category = 'Widening and smoothing'
        filtersize.value = 3

        is_smooth = arcpy.Parameter(
            displayName="Smooth DEM",
            name="is_smooth",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_smooth.category = 'Widening and smoothing'
        is_smooth.value = 'false'

        is_parallel = arcpy.Parameter(
            displayName="Parallel processing",
            name="is_parallel",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        is_parallel.category = 'Parallel processing and tiling'
        is_parallel.value = 'false'

        tile_size = arcpy.Parameter(
            displayName="Tile size",
            name="tile_size",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        tile_size.category = 'Parallel processing and tiling'
        tile_size.value = 2048

        params = [demdataset, marine, output, outputcellsize, minacc1, minlen1, minacc2, minlen2,
                  is_widen, widentype, widendist, filtersize, is_smooth, is_parallel, tile_size]
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

    def call(self, oid):
        try:
            i = int(oid)-1
            arcpy.AddMessage(i)
            raster = 'dem' + str(i)

            arcpy.AddMessage('Counting')
            N = int(arcpy.GetCount_management(self.fishbuffer).getOutput(0))

            arcpy.AddMessage('Selecting cell')
            cells = arcpy.da.SearchCursor(self.fishbuffer, ['SHAPE@', 'OID@'])
            cell = cells.next()
            while cell[1] != oid:
                cell = cells.next()

            arcpy.AddMessage('---')
            arcpy.AddMessage('GENERALIZING DEM ' + str(i + 1) + ' FROM ' + str(N))

            dem0 = arcpy.Raster(self.rastertinworkspace + '/' + raster)
            dem = dem0

            marine_area = None
            process_marine = False
            if self.marine:
                marine_area = self.workspace + "/land" + str(i)
                arcpy.Clip_analysis(self.marine, cell[0], marine_area)
                if int(arcpy.GetCount_management(marine_area).getOutput(0)) > 0:
                    cell_erased = self.workspace + "/cell_erased" + str(i)
                    arcpy.Erase_analysis(cell[0], marine_area, cell_erased)
                    dem = ExtractByMask(dem0, cell_erased)
                    dem.save(self.rastertinworkspace + '/' + raster + "_e")
                    dem = arcpy.Raster(self.rastertinworkspace + '/' + raster + "_e")
                    process_marine = True

            cellsize = dem.meanCellHeight

            arcpy.AddMessage("PREPROCESSING")

            arcpy.AddMessage("Fill...")
            fill = Fill(dem, "")
            arcpy.AddMessage("Dir...")
            dir = FlowDirection(fill, "", "")
            # dir.save(rastertinworkspace + "/dir")
            arcpy.AddMessage("Acc...")
            acc = FlowAccumulation(dir, "", "INTEGER")
            # acc.save(rastertinworkspace + "/acc")

            # MAIN STREAMS AND WATERSHEDS
            arcpy.AddMessage("PROCESSING PRIMARY STREAMS AND WATERSHEDS")

            str1_0 = self.workspace + "/str1"

            arcpy.AddMessage("Extracting primary streams...")

            stream_extractor = ExtractStreams()
            stream_extractor.call(acc, str1_0, self.minacc1, self.minlen1)
            str1 = SetNull(str1_0, 1, "value = 0")

            arcpy.AddMessage("Vectorizing streams...")
            streams1 = self.workspace + "/streams1"
            StreamToFeature(str1, dir, streams1, True)

            arcpy.AddMessage("Deriving endpoints...")
            endpoints1 = self.workspace + "/endpoints1"
            arcpy.FeatureVerticesToPoints_management(streams1, endpoints1, "END")

            endbuffers1 = self.workspace + "/endbuffers1"
            radius = 2 * cellsize
            arcpy.AddMessage("Buffering endpoints...")
            arcpy.Buffer_analysis(endpoints1, endbuffers1, radius, "FULL", "ROUND", "NONE", "")
            rendbuffers1 = self.workspace + "/rendbuffers1"

            arcpy.FeatureToRaster_conversion(endbuffers1, "OBJECTID", rendbuffers1, cellsize)

            mask = CreateConstantRaster(-1, "INTEGER", cellsize, dem.extent)

            arcpy.Mosaic_management(mask, rendbuffers1, "MAXIMUM", "FIRST", "", "", "", "0.3", "NONE")

            arcpy.AddMessage("Erasing streams...")
            str1_e = SetNull(rendbuffers1, str1, "value >= 0")

            arcpy.AddMessage("Vectorizing erased streams...")
            streams1_e = self.workspace + "/streams1_e"
            StreamToFeature(str1_e, dir, streams1_e, True)

            arcpy.AddMessage("Deriving erased endpoints...")
            endpoints1_e = self.workspace + "/endpoints1_e"
            arcpy.FeatureVerticesToPoints_management(streams1_e, endpoints1_e, "END")

            arcpy.AddMessage("Deriving primary watersheds...")
            pour1 = SnapPourPoint(endpoints1_e, acc, cellsize * 1.5, "")
            # pour1.save(rastertinworkspace + "/pour1")

            wsh1 = Watershed(dir, pour1, "")

            arcpy.AddMessage("Vectorizing primary watersheds...")
            watersheds1 = self.workspace + "/watersheds1"
            arcpy.RasterToPolygon_conversion(wsh1, watersheds1, True, "")

            # SECONDARY STREAMS AND WATERSHEDS

            arcpy.AddMessage("PROCESSING SECONDARY STREAMS AND WATERSHEDS")

            arcpy.AddMessage("Extracting secondary streams...")
            str2_0 = self.workspace + "/str2"
            stream_extractor.call(acc, str2_0, self.minacc2, self.minlen2)

            str2 = SetNull(str2_0, 1, "value = 0")
            str2_e = SetNull(str1_0, str2, "value > 0")
            acc_e = SetNull(str1_0, acc, "value > 0")

            arcpy.AddMessage("Vectorizing streams...")
            streams2_e = self.workspace + "/streams2_e"
            StreamToFeature(str2_e, dir, streams2_e, True)

            arcpy.AddMessage("Deriving endpoints...")
            endpoints2_e = self.workspace + "/endpoints2_e"
            arcpy.FeatureVerticesToPoints_management(streams2_e, endpoints2_e, "END")

            arcpy.AddMessage("Buffering primary streams...")
            streambuffer = self.workspace + "/streams1_b"
            arcpy.Buffer_analysis(streams1, streambuffer, radius, "FULL", "ROUND", "NONE", "")

            arcpy.AddMessage("Selecting endpoints...")
            pointslyr = "points"
            arcpy.MakeFeatureLayer_management(endpoints2_e, pointslyr)
            arcpy.SelectLayerByLocation_management(pointslyr, "INTERSECT", streambuffer)

            pourpts2 = arcpy.CreateFeatureclass_management("in_memory", "pourpts2", "POINT", pointslyr, "DISABLED",
                                                           "DISABLED",
                                                           arcpy.Describe(endpoints2_e).spatialReference)
            arcpy.CopyFeatures_management(pointslyr, pourpts2)

            arcpy.AddMessage("Deriving secondary pour pts 1...")
            pour21 = SnapPourPoint(pourpts2, acc_e, cellsize * 1.5, "")

            arcpy.AddMessage("Deriving secondary pour pts 2...")
            pour22 = SnapPourPoint(pourpts2, acc_e, cellsize * 2, "")

            arcpy.AddMessage("Mosaic secondary pour pts...")
            arcpy.Mosaic_management(pour21, pour22, "FIRST", "FIRST", "0", "0", "", "0.3", "NONE")

            # pour22.save(rastertinworkspace + "/pour22")

            arcpy.AddMessage("Deriving secondary watersheds...")
            wsh2 = Watershed(dir, pour22, "")
            # wsh2.save(rastertinworkspace + "/wsh2")

            arcpy.AddMessage("Vectorizing secondary watersheds...")
            watersheds2 = self.workspace + "/watersheds2"
            arcpy.RasterToPolygon_conversion(wsh2, watersheds2, True, "")

            arcpy.AddMessage("Interpolating features into 3D...")

            streams1_3d = self.workspace + "/streams1_3d"

            arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, streams1, streams1_3d)

            watersheds1_3d = self.workspace + "/watersheds1_3d"
            arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, watersheds1, watersheds1_3d)

            watersheds2_3d = self.workspace + "/watersheds2_3d"
            arcpy.InterpolateShape_3d(dem0.path + '/' + dem0.name, watersheds2, watersheds2_3d)

            marine_3d = self.workspace + "/marine_3d"
            if process_marine:
                arcpy.InterpolateShape_3d(self.demdataset, marine_area, marine_3d)

            # GENERALIZED TIN SURFACE

            arcpy.AddMessage("DERIVING GENERALIZED SURFACE")

            arcpy.AddMessage("TIN construction...")

            tin = self.rastertinworkspace + "/tin"
            features = []
            s1 = "'" + streams1_3d + "' Shape.Z " + "hardline"
            w1 = "'" + watersheds1_3d + "' Shape.Z " + "softline"
            w2 = "'" + watersheds2_3d + "' Shape.Z " + "softline"

            features.append(s1)
            features.append(w1)
            features.append(w2)
            if process_marine:
                m2 = "'" + marine_3d + "' Shape.Z " + "hardline"
                features.append(m2)

            featurestring = ';'.join(features)

            arcpy.ddd.CreateTin(tin, "", featurestring, "")

            # GENERALIZED RASTER SURFACE
            arcpy.AddMessage("TIN to raster conversion...")

            rastertin = self.rastertinworkspace + "/rastertin"
            try:
                arcpy.TinRaster_3d(tin, rastertin, "FLOAT", "NATURAL_NEIGHBORS", "CELLSIZE " + str(cellsize), 1)
            except Exception:
                arcpy.AddMessage("Failed to rasterize TIN using NATURAL_NEIGHBORS method. Switching to linear")
                arcpy.TinRaster_3d(tin, rastertin, "FLOAT", "LINEAR", "CELLSIZE " + str(cellsize), 1)
            # POSTPROCESSING

            arcpy.AddMessage("POSTPROCESSING")

            # Widen valleys and ridges
            widenraster = rastertin
            if self.is_widen == "true":
                arcpy.AddMessage("Raster widening...")
                widenraster = self.workspace + "/widenraster"
                dem_widener = WidenLandforms()
                dem_widener.call(rastertin, streams1, self.widendist, self.filtersize, widenraster, self.widentype)

            # Smooth DEM
            result = arcpy.Raster(widenraster)
            if self.is_smooth == "true":
                arcpy.AddMessage("Raster filtering...")
                neighborhood = NbrRectangle(self.filtersize, self.filtersize, "CELL")
                result = FocalStatistics(widenraster, neighborhood, "MEAN", "DATA")

            if process_marine:
                arcpy.AddMessage("Masking marine regions...")
                result_erased = ExtractByMask(result, cell_erased)
                arcpy.Mosaic_management(result_erased, rastertin, "FIRST", "FIRST", "", "", "", "0.3", "NONE")
                arcpy.AddMessage("Saving result...")
                res = arcpy.Raster(rastertin)
                res.save(self.rastertinworkspace + '/gen/dem' + str(i))
            else:
                arcpy.AddMessage("Saving result...")
                result.save(self.rastertinworkspace + '/gen/dem' + str(i))

            arcpy.Delete_management(pointslyr)

            self.i += 1

            if self.i == self.N:
                arcpy.AddMessage("---")
                arcpy.AddMessage("CLEANING TEMPORARY DATA")

                # arcpy.Delete_management(str1_0)
                # arcpy.Delete_management(str2_0)
                # arcpy.Delete_management(streams1)
                # arcpy.Delete_management(endpoints1)
                # arcpy.Delete_management(endbuffers1)
                # arcpy.Delete_management(rendbuffers1)
                # arcpy.Delete_management(streams1_e)
                # arcpy.Delete_management(endpoints1_e)
                # arcpy.Delete_management(watersheds1)
                # arcpy.Delete_management(streams2_e)
                # arcpy.Delete_management(endpoints2_e)
                # arcpy.Delete_management(streambuffer)
                # arcpy.Delete_management(pourpts2)
                # arcpy.Delete_management(watersheds2)
                # arcpy.Delete_management(streams1_3d)
                # arcpy.Delete_management(watersheds1_3d)
                # arcpy.Delete_management(watersheds2_3d)
                # arcpy.Delete_management(tin)
                # arcpy.Delete_management(rastertin)
                # arcpy.Delete_management(widenraster)
                #
                # arcpy.Delete_management(fill)
                # arcpy.Delete_management(dir)
                # arcpy.Delete_management(acc)
                # arcpy.Delete_management(str1)
                # arcpy.Delete_management(str2)
                # arcpy.Delete_management(mask)
                # arcpy.Delete_management(str1_e)
                # arcpy.Delete_management(str2_e)
                # arcpy.Delete_management(pour1)
                # arcpy.Delete_management(pour21)
                # arcpy.Delete_management(pour22)
                # arcpy.Delete_management(wsh1)
                # arcpy.Delete_management(wsh2)
                # arcpy.Delete_management(acc_e)

        except Exception:
            arcpy.AddMessage("Failed to generalize " + raster)

    def execute(self, parameters, messages):
        self.demdataset = parameters[0].valueAsText
        self.marine = parameters[1].valueAsText
        output = parameters[2].valueAsText
        outputcellsize = float(parameters[3].valueAsText)
        self.minacc1 = int(parameters[4].valueAsText)
        self.minlen1 = int(parameters[5].valueAsText)
        self.minacc2 = int(parameters[6].valueAsText)
        self.minlen2 = int(parameters[7].valueAsText)
        self.is_widen = parameters[8].valueAsText
        self.widentype = parameters[9].valueAsText
        self.widendist = float(parameters[10].valueAsText)
        self.filtersize = int(parameters[11].valueAsText)
        self.is_smooth = parameters[12].valueAsText
        is_parallel = parameters[13].valueAsText
        tilesize = int(parameters[14].valueAsText)

        self.workspace = os.path.dirname(output)

        # raster workspace MUST be a folder, no
        self.rastertinworkspace = self.workspace
        n = len(self.rastertinworkspace)
        if n > 4:
            end = self.rastertinworkspace[n - 4: n]  # extract last 4 letters
            if end == ".gdb":  # geodatabase
                self.rastertinworkspace = os.path.dirname(self.rastertinworkspace)

        arcpy.AddMessage(self.rastertinworkspace)

        arcpy.CreateFolder_management(self.rastertinworkspace, 'scratch')
        self.rastertinworkspace += '/scratch'

        arcpy.CreateFolder_management(self.rastertinworkspace, 'gen')
        arcpy.CreateFolder_management(self.rastertinworkspace, 'gencrop')

        arcpy.env.scratchWorkspace = self.rastertinworkspace
        arcpy.env.workspace = arcpy.env.scratchWorkspace

        demsource = arcpy.Raster(self.demdataset)

        nrows = math.ceil(demsource.height / tilesize)
        ncols = math.ceil(demsource.width / tilesize)
        total = nrows * ncols

        cellsize = 0.5 * (demsource.meanCellHeight + demsource.meanCellWidth)

        bufferpixelwidth = math.ceil(max(demsource.width, demsource.height) / (max(nrows, ncols) * 10))
        bufferwidth = bufferpixelwidth * cellsize

        arcpy.AddMessage('Splitting raster into ' + str(nrows) + ' x ' + str(ncols) + ' = ' + str(total) + ' tiles')
        arcpy.AddMessage('Tile overlap will be ' + str(2 * bufferpixelwidth) + ' pixels')

        fishnet = self.workspace + "/fishnet"
        arcpy.CreateFishnet_management(fishnet,
                                       str(demsource.extent.XMin) + ' ' + str(demsource.extent.YMin),
                                       str(demsource.extent.XMin) + ' ' + str(demsource.extent.YMin + 1),
                                       '', '',
                                       nrows, ncols,
                                       '', '',
                                       demsource.extent, 'POLYGON')

        self.fishbuffer = self.workspace + "/fishbuffer"
        arcpy.Buffer_analysis(fishnet, self.fishbuffer, bufferwidth)

        fishmaskbuffer = self.workspace + "/fishmaskbuffer"
        arcpy.Buffer_analysis(fishnet, fishmaskbuffer, max(demsource.meanCellHeight, demsource.meanCellWidth))

        arcpy.SplitRaster_management(self.demdataset,
                                     self.rastertinworkspace,
                                     'dem',
                                     'POLYGON_FEATURES',
                                     'GRID',
                                     '', '', '', '', '', '', '',
                                     self.fishbuffer)

        # rasters = arcpy.ListRasters("*", "GRID")

        rows = arcpy.da.SearchCursor(self.fishbuffer, 'OID@')
        oids = [row[0] for row in rows]

        # MAIN PROCESSING
        if is_parallel == 'true':

            arcpy.AddMessage('Trying to make multiprocessing')

            execute(*oids)

        else:
            for oid in oids:
                self.call(oid)

        arcpy.AddMessage("CLIPPING ANS MASKING GENERALIZED RASTERS")

        rows = arcpy.da.SearchCursor(fishmaskbuffer, ['SHAPE@', 'OID@'])
        i = 0
        for row in rows:
            dem = arcpy.Raster(self.rastertinworkspace + '/gen/dem' + str(i))
            dem_clipped = ExtractByMask(dem, row[0])
            dem_clipped.save(self.rastertinworkspace + '/gencrop/dem' + str(i))
            i += 1

        arcpy.env.workspace = self.rastertinworkspace + '/gencrop/'

        rasters = arcpy.ListRasters("*", "GRID")

        rasters_str = ';'.join(rasters)

        arcpy.MosaicToNewRaster_management(rasters_str,
                                           os.path.dirname(output),
                                           os.path.basename(output),
                                           "",
                                           "16_BIT_SIGNED",
                                           str(outputcellsize),
                                           "1",
                                           "BLEND",
                                           "FIRST")

        return