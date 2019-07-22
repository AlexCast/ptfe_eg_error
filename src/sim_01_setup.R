
###################
# Load libraries: #
###################

library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(akima)
library(abind)
source("src/sim_02_read_opac.R")
source("src/sim_03_geometrical_functions.R")
source("src/sim_04_tilt.R")
source("src/sim_05_surface_reflection.R")
source("src/sim_06_atmosphere.R")
source("src/sim_07_sky_radiance_model.R")
source("src/sim_08_plaque_brdf.R")

##################
# Set constants: #
##################

crs.wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
crs.ste <- "+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
crs.ort <- "+proj=ortho +lat_0=90 +datum=WGS84 +units=m +no_defs"

#########################################
# Create upper hemisphere Bounding box: #
#########################################

bnds <- cbind(x = c(-180, -180, 180, 180, -180), y = c(90, 0, 0, 90, 90))
uhbb <- SpatialPolygons(list(Polygons(list(Polygon(bnds)), "1")),
                 proj4string=CRS(crs.wgs))

##############################
# Create zenith angle lines: #
##############################

LinesList <- list(Line(cbind(x = -180:180, y = rep(0, length(-180:180)))), 
                  Line(cbind(x = -180:180, y = rep(50, length(-180:180)))),
                  Line(cbind(x = -180:180, y = rep(75, length(-180:180)))))
sl        <- Lines(LinesList, ID = 1)
zlines    <- SpatialLines(list(sl), proj4string = CRS(crs.wgs))
zlin_ste  <- spTransform(zlines, crs.ste)
zlin_ort  <- spTransform(zlines, crs.ort)
