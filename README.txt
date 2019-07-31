The code available in this directory is released under license GPL Version 3 
(or higher). See LICENCE_gpl-3.0.txt for more infomration.


##########
# README #
##########

This directory contain code and data as supplementary material to the research article:

Castagna, A., Johnson, B. C., Voss, K., Dierssen, H. M., Patrick, H., Germer, T., Sabbe, 
K., & Vyverman, W. 2019. Uncertainty in global downwelling plane irradiance estimates from 
sintered polytetrafluoroethylene plaque radiance measurements. Applied Optics XX, XXX.

A basic knowledge of the R programing language is assumed. To load the data, results and 
functions, is necessary to run:
load('restore.RData')

To manipulate data and use the provided plot functions for visualization of BRDF and sky 
radiance, it is necessary an R version >= 3.3.0 and packages abind, sp and raster. Those 
can be installed from an R session with the command:

install.packages(c('sp', 'raster', 'abind'))

To run the simulations, additionally the R packages rgeos, rgdal and akima are necessary.
To install all required packages, use:

install.packages(c('sp', 'raster', 'rgeos', 'rgdal', 'akima', 'abind'))

In Windows, it is recommended to install first the free and open source Quantum GIS 
software (https://www.qgis.org/en/site/), that will provide a compiled binary of GDAL, 
necessary for the rgdal package. This package is used only to read the vector shadow masks 
from disk.

Parts of the code for the simulations are written in C (BRDF functions). The functions 
were extracted from the SCATMECH C++ library. Complementary compiled C programs were 
compiled in Debian GNU/Linux in an AMD64 architecture, using gcc version 6.3.0 and may not 
work properly in Windows. C source code is provided for new compilation as necessary. The 
header in each C source code file give instructions for compilation in GNU/Linux with gcc.

Note that to reduce data volume, only the results of the analysis are provided directly.
The sky radiances distributions and the BRDFs used however, can be recreated by running 
the simulation code. Complete calculation can take up to three weeks in a common desktop 
computer. To allow exploration and test of the code, sky radiances for continental clean
aerosol, with 95% humidity and aerosol optical thickness at 550 nm of 0.1 are provided. 
For the same reason, the pre-computed BRDF for the level plaques is also provided.

All source codes were saved with a Windows type line ending (\r\n).

###################
# Directory tree: #
###################

ROOT
|- ANC: ancillary data
|   |- OPAC_INPUT: OPAC input files for the aerosol models
|   |- OPAC_OUT: OPAC output files for the aerosol models
|   |- PTFE_BRDF_COEF: zernike coefficients for PTFE BRDF interpolation, including the new
|                      coefficients for the grey plaque at 532 nm
|
|- BRDF_C: C code for the PTFE BRDF Zernike polynomial interpolation
|- FIELD: Data from the inter-comparison experiment
|- PLOT: plots of results, in pdf format
|- RES: analysis results in R binary data format
|- SIM: simulation data
|   |- BRDF: pre-computed BRDF for both plaques and view geometries, both level and tilt
|   |    |- TILT:
|   |         |- GEOMETRY_1: pre-computed BRDF for both plaques under nadir view with tilt
|   |         |    |- 03
|   |         |    |- 06
|   |         |    |- 09
|   |         |    |- 12
|   |         |
|   |         |- GEOMETRY_2: pre-computed BRDF for both plaques under side view with tilt
|   |              |- 03
|   |              |- 06
|   |              |- 09
|   |              |- 12
|   |
|   |- SHADOW: shadow vector for each geometry
|   |- SKYRAD: pre-computed sky radiance distribution and Es at surface
|
|- SRC: R source code for simulation and analysis

######################
# Measurements data: #
######################

Two sets of measurements are made available: (1) grey plaque BRDF and (2) Irradiance from 
field inter-comparison experiments. The grey plaque BRDF data were obtained and 
interpolated from a 10 % reflective PTFE sample at 532 nm using the methods described in 
Germer (2017). The Zernike expansion coefficients based on the BRDF data are included in 
ANC/PTFE_BRDF_COEF/PTFE.grey10.532.coef.csv, Those are suitable for use by the 
Zernike_Expansion_BRDF_Model in the SCATMECH library.

Irradiance data estimated with the plaque method and measured with the irradiance sensors
are provided in the FIELD directory. Those are available in two equivalent sets for nadir
view and for side view. Eg measurements for the irradiance sensors are averaged for the 
duration of the plaque measurements. Data has units of W/m2/nm. The position of the 
stations is provided as a KML file.

####################
# Simulation data: #
####################

The main data provided is the Eg_plaque_estimate / Eg_real ratio for BRDF effect alone, 
shadow effect alone, BRDF + shadow, tilt effect alone and BRDF + shadow + tilt. Data is 
organized as multidimensional arrays and is saved in R binary format. It will be already
loaded for an R session initiated at the root folder (R session is saved as binary data 
hidden file and automatically loaded upon start). The datasets can also be loaded in any 
R session with the function 'load':

load('full/path.Rdata')

The ratio data is saved in the RES directory (aka, results). The general naming convention 
is 'res_effect_component.Rdata'. 'effect' refers to BRDF if 'brdf', shadow if 'shdw', 
tilt if 'tilt' or combined affects if 'all'. 'component' refers to Eg if 'tot', Es if 
'dir' and Ed if 'dif'. When the lambertian assumption is used, '_pi' is added to the file 
and object name. The common dimensions of the data are:

- Wavelength [66]: from 350 to 1000 nm in steps of 10 nm.
- Sun zenith angle [10]: Overcast, then from 0 to 80 degrees, in steps of 10 degrees.
- Plaque [2]: Plaque type, white and grey.
- View geometry [2]: Nadir and 40 degrees.
- Sky model [24]: Sky model dependent on aerosol model.
- Zenith tilt [4]: Tilt from local zenith from 3 to 12 degrees, in steps of 3 degrees.
- Azimuth tilt [8]: Tilt from Sun direction from 0 to 315 degrees, in steps of 45 degrees.

The files and their content and dimensions are described next:

res_brdf_tot.Rdata, res_brdf_tot_pi.Rdata : 
	Ratios for the BRDF effect alone. Array with dimensions:
	(1) wavelength [66]
	(2) sun zenith angle [10] 
	(3) plaque [2] 
	(4) geometry [2] 
	(5) sky model [24] 

res_shdw_tot_nbrdf.Rdata, res_shdw_dir_nbrdf.Rdata, res_shdw_dif_nbrdf.Rdata :
	Ratios of the shadow effect allone. Note that since BRDF effects are not 
	considered, there is no need to calculate per plaque or per retrieval equation. 
	Array with dimensions:
	(1) wavelength [66]
	(2) sun zenith angle [10] 
	(3) geometry [2] 
	(4) sky model [24]

res_shdw_tot.Rdata, res_shdw_dir.Rdata, res_shdw_dif.Rdata, res_shdw_tot_pi.Rdata :
	Ratios for the BRDF and shadow effects (interaction). Array with dimensions:
	(1) wavelength [66]
	(2) sun zenith angle [10] 
	(3) plaque [2] 
	(4) geometry [2] 
	(5) sky model [24] 

res_tilt_tot.Rdata, res_tilt_dir.Rdata, res_tilt_dif.Rdata, res_tilt_tot_pi.Rdata :
	Ratios of the tilt effect alone. Note that since BRDF and shadow effects are not
	included, there is no need to calculate per plaque, per geometry or per retrieval
	equation. Dimensions are:
	(1) zenith tilt [4]
	(2) azimuth tilt [8]
	(3) wavelength [66]
	(4) sun zenith angle [10] 
	(5) sky model [24]

res_all_tot.Rdata, res_all_dir.Rdata, res_all_dif.Rdata, res_all_tot_pi.Rdata :
	Ratios for the BRDF, shadow and tilt effects (interaction). Array with dimensions:
	(1) zenith tilt [4]
	(2) azimuth tilt [8]
	(3) wavelength [66]
	(4) sun zenith angle [10] 
	(5) plaque [2]
	(6) geometry [2]
	(7) sky model [24]

res_brdfg.Rdata :
	The weighted average BRDF for the total incident light (BRDFg), when considering 
	only shadow effects (no tilt). Array with dimensions:
	(1) wavelength [66]
	(2) sun zenith angle [10] 
	(3) plaque [2] 
	(4) geometry [2] 
	(5) sky model [24]

crosscal_shdw_tot_pi.Rdata and crosscal_all_tot_pi.Rdata :
	Specific data to simulate the effect of cross-calibration of the white and grey
	plaque before deployment. The first considers only BRDF and shadow effects and the
	second all the effects. Dimensions are equal to res_shdw_tot.Rdata for the first 
	and res_all_tot.Rdata for the second. Only calculated with equation 5B as was done
	with field measurements.

Other data made available is the sky radiance distribution and pre-computed BRDF for each
plaque, view geometry and tilt condition (from 0 to 12 degrees, all azimuths). Sky 
radiance is provided in long data format, in text files. BRDFs are provided in R binary 
format. Little guidance will be provided on how to manipulate those data since they are 
not the primary result of the study. Note also that in order to reduce data volume, the 
sky radiance distributions and the BRDFs are not directly provided and must be recreated 
by running the provided simulation code. Interested users can inspect the code for 
examples on how to import and use this data, after they are generated.

####################
# Data extraction: #
####################

A source code with examples for data manipulation is provided in root directory. It show 
examples of plot functions ('plot_brdf' and 'plot_skyrad') and how to extract error 
statistics for the ratios results.

###################
# Simulation run: #
###################

To re-run or run new simulations (e.g., specific conditions during field deployment) an R 
session should be started at the root directory (as a R working directory). R source code 
names start with 'sim' if provide base function for simulation of sky radiance, tilt and 
BRDF. The source code 'sim_01_setup.R' will load all required libraries and load all 
'sim' source code to the R session. The file 'sim_09_production.R' has the code that uses
those functions to carry out the sky radiance simulations.

The files starting with 'analysis' show the sequence of analysis steps, from defining 
fixed grid summation approximation to integral to pre-calculating the BRDFs and finally to
calculate the Eg_plaque_estimate / Eg_real ratios and extract statistics of errors. The 
source code for each component of the analysis can be used as example and changed as 
necessary (e.g, analysis_03_brdf.R, analysis_04_shadow.R, analysis_05_tilt.R, 
analysis_06_all.R). Of particular interest might be the code in 
'analysis_08_field_condition.R', that shows the analysis of errors after cross-calibration
with a white plaque. The code can be easily adapted for other cross-calibration 
conditions.


