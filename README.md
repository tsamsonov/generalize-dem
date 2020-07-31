[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3813731.svg)](https://doi.org/10.5281/zenodo.3813731) [![Python 3.6](https://img.shields.io/badge/python-3.6-red.svg)](https://www.python.org/downloads/release/python-360/) [![Python 2.7](https://img.shields.io/badge/python-2.7-orange.svg)](https://www.python.org/downloads/release/python-270/) [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%203%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)

# Generalization and conflation of digital elevation models

**generalize-dem** toolbox contains tools for structural generalization and conflation of digital elevation models.

## Requirements

You need **ArcGIS Pro** (Python 3.6) or **ArcGIS for Desktop 10.3+** (Python 2.7) with *Spatial Analyst* and *3D Analyst* extension modules to use the toolbox.

## Installation

Download the latest [release](https://github.com/tsamsonov/generalize-dem/releases) and extract the contents of ZIP archive. You should see `generalize-dem.pyt` Python toolbox in the *Catalog* window pane of **ArcGIS Pro** or **ArcMap**:

![toolbox](img/toolbox.png)

## Usage

Ten tools are contained in the toolbox:

1. **Сalculate Distances Between Lines** tool calculates Directed Hausdorff, Hausdorff and Frechét distances between correspoding lines in two feature classes.

2. **Carve DEM along streams** tool changes cell values in DEM so that elevation is monotonoysly decreasing along given lines.

3. **Conflate DEM by Links** tool performes rubbersheeting of DEM so that it becomes spatially adjusted with a given set of reference hydrographic lines.

4. **Create Fishnet** tool generates rectangular grids with overlapping or non-overlapping cells.

5. **Extract Counterpart Streams** tool finds path in drainage network that are similar to given reference hydrographic lines.

6. **Extract Streams** tool traces the streams using flow accumulation and length criteria.

7. **Filter DEM** tool performs filtering of DEM.

8. **Generalize DEM** tool performs structural generalization of raster digital elevation model

9. **Generate Conflation Links** tool generates conflation links between counterpart streams and reference hydrographic lines.

10. **Widen Landforms** tool allows widening of negative and positive terrain features, which can be very effective to improve visual analysis of generalized DEM.

## Example data

You can find an example DEM in `dem.gdb` geodatabase (please unzip gdb before using). Source dem is called simply `dem` and generalized is called `dem_gen`.

## Acknowledgements

[pybind11](https://github.com/pybind/pybind11/) library is used to compile stream extraction module.

## Further reading

A detailed description of the method can be found in the following papers:

Samsonov, T. (2011) *Multiscale Hypsometric Mapping*, in Ruas, A. (ed.) **Advances in Cartography and GIScience. Vol.1: Selection from ICC-2011, Paris**. Berlin, Heidelberg: Springer (Lecture Notes in Geoinformation and Cartography), pp. 497–520. DOI: [10.1007/978-3-642-19143-5_28](http://dx.doi.org/10.1007/978-3-642-19143-5_28)

Samsonov, T. (2020) *Automated Conflation of Digital Elevation Model with Reference Hydrographic Lines*, **ISPRS International Journal of Geo-Information**, 9(5), 334. DOI: [10.3390/ijgi9050334](http://dx.doi.org/10.3390/ijgi9050334).

## Citation & Copyright

To cite the software use the following reference:

Samsonov, T. (2020). *Generalize DEM: ArcGIS Python toolbox for automated structural generalization and conflation of digital elevation models*. **Zenodo**. DOI: [10.5281/zenodo.3813731](https://doi.org/10.5281/zenodo.3813731)

## Funding 

Development of this software is funded by Russian Science Foundation (RSF) grant No 19-77-00071

© 2010-2020, Timofey Samsonov, Lomonosov MSU Faculty of Geography.