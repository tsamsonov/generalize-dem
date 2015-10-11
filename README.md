# Small-Scale-Terrain-Generalization
DEM Generalization toolbox has been developed for small-scale terrain generalization (production scales 1:1 000 000 and smaller). Resulting DEM is best suited for hypsometric terrain representation (classified image in GIS software).

You can find an example DEM in dem.gdb geodatabase (please unzip gdb before using). Source dem is called simply "dem" and generalized is called "dem_gen". Generalization parameters for Generalize DEM tool can be found in Parameters_Example.png screenshot.

The toolbox consists of four tools/scripts: Generalize DEM, Widen Landforms, Flow Accumulation to Streams, Filter DEM. The associated Python scripts must be in the same folder as the toolbox itself.

You can find an extensive explanation of the tool parameters directly in the tool dialog.

GENERALIZE DEM tool performs structural generalization of raster digital elevation model for small-scale mapping (1:500 000 and smaller). The underlying algorithm workflow is presented in the figure below. It is based on the idea that the generalized terrain surface can be reconstructed from skeletal lines such as streams and watersheds. The derived skeleton includes three types of lines: streams representing the axes of negative forms; watershed boundaries representing the axes of positive forms; and secondary watersheds representing the slopes of negative and positive terrain features. Streams are derived using Flow Accumulation to Streams tool (see this toolbox) using the minimal stream length that is representative for the resulting scale. Watersheds are simply derived between the constructed streams. The secondary watersheds are derived for the supplementary secondary stream layer which is derived by Flow Accumulation to Streams tool too, but with minimal stream length that is representative for the source DEM. After skeletal lines are derived, they are triangulated, and the resulting TIN surface is rasterized. Resulting raster is then post-processed by Widen Valleys and Hills tool (see this toolbox) to derive clearer representation of negative and positive terrainn features.

For more details and case studies see the following paper:

Samsonov T. Multiscale hypsometric mapping // Advances in Cartography and GIScience, Vol. 1. — Lecture Notes in Geoinformation and Cartography. — Springer: Heidelberg, Germany, 2011. — P. 497–520. http://dx.doi.org/10.1007/978-3-642-19143-5_28

Note that this tool is not suitable for generalization of a detailed DEM in large- and middle-scale mapping (larger than 1:500 000). For these scales try using Terrain Equalizer by Bernhard Jenny (http://www.terraincartography.com/terrainequalizer/).

It is not recommended to use this toolbox for generalization of high-resolution models at large and middle scales, because this scale range presumes different principles of terrain generalization.

WIDEN LANDFORMS tool is intended to facilitate terrain image readability in small scales. Valley wideining is a standard operation that is used in small-scale terrain generalization and is applied at the post-processing stage. If you use contours as representation technique, then contours would spread apart from the streams and the valley bottom would become more prominent in the image. The same procedure is often applied to the positive landforms - hills and ridges. This is why the tool is called Widen Landforms. 

At first, the algorithm produces 2 derivative rasters using MIN and MAX filters correspondingly. Then SOURCE, MIN and MAX rasters are added together with weights set according to the euclidian distance from the nearest stream to the current pixel.

This tool is a modification of DEM generalization algorithm introduced in this article:

Anna M. Leonowicz, Bernhard Jenny, Lorenz Hurni, Automatic generation of hypsometric layers for small-scale maps, Computers & Geosciences, Volume 35, Issue 10, October 2009, Pages 2074-2083, ISSN 0098-3004, http://dx.doi.org/10.1016/j.cageo.2008.12.012

The difference is that original methodology used quartile filters instead of min/max and exxagerated valleys only.

FLOW ACCUMULATION TO STREAMS tool extracts and generalizes stream network from flow accumulation raster. It has two boundary conditions: minimum flowacc value that defines how much uphill stream sources can be and minimum length which allows omission of short streams. This tool is an alrtenative to standard stream network extraction workflow (also described in Esri's Spatial Analyst help) which is based on flow accumulation only. Instead of applying simple condition (Con("Value" > 100), for example) this script traces upstream from every pixel until the value is below the flowacc treshold specified. Then, if the length of the traced path is larger than minimum length specified, all traced pixels are marked as streams. This approach works much better than simple condition, because it does not clip upper reaches.The script is optimized not to trace already traced paths. 

For more details of the algorithm please read this article:

Anna M. Leonowicz, Bernhard Jenny, Lorenz Hurni, Automatic generation of hypsometric layers for small-scale maps, Computers & Geosciences, Volume 35, Issue 10, October 2009, Pages 2074-2083, ISSN 0098-3004, http://dx.doi.org/10.1016/j.cageo.2008.12.012

FILTER DEM tool implements focal statistics filtering for raster datasets. It includes several standard filters including mean, median, min, max and also low and high quartile filters that can be useful in DEM generalization. Another useful feature is the possibility to make several iterations

2015, Timofey Samsonov, Lomonosov MSU Faculty of Geography
