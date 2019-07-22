# Eg estimation from PTFE plaques
#
# Analysis part 1 - Setup shadows
#
# This script set up shadow geometries for nadir view and 40º at azimuths of 90 degrees. 
# Shadows shapes are calculated for when the plaque is hold at upper chest/neck high, at 
# 40 cm from the body and with the sensor at least at 30 cm distance from the surface.
#
# Nadir view shadow is coded here as 'shadow_1' and at 40º zenith as 'shadow_2'. Those are 
# SpatialPolygons as defined in the sp package. Those objects are used with the raster
# package as masks to set to NA the radiance grid values for the directions the shadows
# intersect. 
#
# The procedure below manually creates the shadows vertices using geometry. The shadow is 
# drawn for both view angles at an azimuth of 90 degrees. They are densified and corrected 
# for projection distortions. The code below is a bit complex but of little importance as 
# the shadow spatial polygon is the only final product. It is expected that tilt of the 
# shadow as a spatial polygon would be challenging, since vertices could cross (specially 
# in nadir view) with crossing of line segments. The procedure that will be used is to 
# first apply the shadow mask to the untilted radiance distribution, and then tilt it. 
# Because the tilt expose the plaque view to the lower hemisphere, shadow is calculated 
# below 90 degrees zenith too.
#
# To load the created shadows at the recommended side view, do:
#
# library(sp)
# shadow_1 <- readOGR("sim/shadow", "shadow1_090")
# shadow_2 <- readOGR("sim/shadow", "shadow2_090")
#
# Version: 2018-07-08
#
# Alexandre Castagna
# alexandre.castagna@ugent.be
# Protistology and Aquatic Ecology
# Biology Department, Faculty of Sciences, Gent University
# Krijgslaan, 281, S8, 3rd floor, C-wing, Gent 9000, Belgium - BE
#

#
# Load libraries:
#

library(sp)
library(raster)

#
# Compute shadow manually:
#

dir.create("sim/shadow")

##############################
# Set up of the base shadow: #
##############################

# 40º off-nadir view:

x1 <- cos(rad(0:360)) * 10.62								# Head
y1 <- (sin(rad(0:360)) * 10.62 * 1.4) + 21.992

x2 <- c(-x1[57], seq(-x1[57], x1[57], length.out = 50), x1[57])				# Neck
y2 <- c(y1[305], rep(0, 50), y1[305])

x3 <- cos(rad(0:360)) * 4.5								# Hand
y3 <- (sin(rad(0:360)) * 4.5) + 50

x4 <- c(2.5, 20, seq(17, 27, length.out = 50), 30, 4.4) 				# Arm
y4 <- c(46.25, 10, rep(0, 50), 10, 49)

x5 <- c(-27, -25, -23, seq(-20, 20, 1), 23, seq(27, -27, -1))   			# Torso
y5 <- c(0, -10, -13, rep(-20, 41), -10, rep(0, 55))

plot(x1, y1, asp = 1, ylim = c(-20, 55), xlim = c(-30, 30), type = "l")
lines(x2, y2)
lines(x3, y3)
lines(x4, y4)
lines(x5, y5)

xf1 <- c(x1[1:237], x2, x1[305:361])
yf1 <- c(y1[1:237], y2, y1[305:361])
xf2 <- c(x3[1:304], x4, x3[348:361])
yf2 <- c(y3[1:304], y4, y3[348:361])
xf3 <- x5
yf3 <- y5

# Make spatial:
xf1 <- xf1 / sin(rad(90 - yf1)) # Sine distortion with elevation...
xf2 <- xf2 / sin(rad(90 - yf2))
xf3 <- xf3 / sin(rad(90 - abs(yf3)))

xf1 <- xf1 + 90
xf2 <- xf2 + 90
xf3 <- xf3 + 90

Sr1   <- Polygon(cbind(xf1, yf1))
Sr2   <- Polygon(cbind(xf2 ,yf2))
Sr3   <- Polygon(cbind(xf3 ,yf3))
Srs1  <- Polygons(list(Sr1), "s1")
Srs2  <- Polygons(list(Sr2), "s2")
Srs3  <- Polygons(list(Sr3), "s3")
geom2 <- SpatialPolygons(list(Srs1, Srs2, Srs3), 1:3, proj4string = CRS(crs.wgs))

plot(spTransform(geom2[1:2], crs.ste), col = "black")

# Didn't work quite well in the end, when masked raster is projected it show some
# distortions where a simple line should be present. Correction for it is a bit 
# cumbersome, but... :

geom2ste <- spTransform(geom2, crs.ste)
temp <- geom2ste

temp <- temp[2]@polygons[[1]]@Polygons[[1]]@coords
data <- as.data.frame(temp[15:16, ])
colnames(data) <- c("x", "y")
tempxa <- seq(data[1, 1], data[2, 1], length.out = 50)
tempya <- predict(lm(y~x, data), data.frame(x = tempxa))
data <- as.data.frame(temp[16:17, ])
colnames(data) <- c("x", "y")
tempxb <- seq(data[1, 1], data[2, 1], length.out = 50)
tempyb <- predict(lm(y~x, data), data.frame(x = tempxb))
data <- as.data.frame(temp[66:67, ])
colnames(data) <- c("x", "y")
tempxc <- seq(data[1, 1], data[2, 1], length.out = 50)
tempyc <- predict(lm(y~x, data), data.frame(x = tempxc))
data <- as.data.frame(temp[67:68, ])
colnames(data) <- c("x", "y")
tempxd <- seq(data[1, 1], data[2, 1], length.out = 50)
tempyd <- predict(lm(y~x, data), data.frame(x = tempxd))
temp   <- rbind(temp[1:15, ], cbind(tempxa, tempya)[-1, ], cbind(tempxb, tempyb)[-1, ], 
		temp[18:66, ], cbind(tempxc, tempyc)[-1, ], cbind(tempxd, tempyd)[-1, ], 
		temp[69:372, ])
rownames(temp) <- NULL
colnames(temp) <- c("x", "y")
Sr2p   <- Polygon(temp)
Srs2p  <- Polygons(list(Sr2p), "s2")
Srs4p  <- Polygons(list(Sr2p), "s4")
Srs2p.sp <- SpatialPolygons(list(Srs2p, Srs4p), 1:2, proj4string = CRS(crs.ste))
Srs2p.sp <- spTransform(Srs2p.sp, crs.wgs)
temp  <- Srs2p.sp[1]@polygons[[1]]@Polygons[[1]]@coords
Sr2   <- Polygon(temp)
Srs2  <- Polygons(list(Sr2), "s2")

geom2 <- SpatialPolygons(list(Srs1, Srs2, Srs3), 1:3, proj4string = CRS(crs.wgs))
plot(spTransform(geom2[1:2], crs.ste), col = "black")


# Nadir view:

x3 <- cos(rad(0:360)) * 4.5								# Hand
y3 <- (sin(rad(0:360)) * 4.5) + 90

x4 <- c(2.5, 20, seq(17, 27, length.out = 50), 30, 4.4)					# Arm
y4 <- c(86.25, 30, rep(0, 50), 30, 89)

plot(x1, y1, asp = 1, ylim = c(-20, 95), xlim = c(-30, 30), type = "l")
lines(x2, y2)
lines(x3, y3)
lines(x4, y4)
lines(x5, y5)

xf1 <- c(x1[1:237], x2, x1[305:361])
yf1 <- c(y1[1:237], y2, y1[305:361])
xf2 <- c(x3[1:304], x4, x3[348:361])
yf2 <- c(y3[1:304], y4, y3[348:361])
xf3 <- x5
yf3 <- y5

# Make spatial:
xf1 <- xf1 / sin(rad(90 - yf1)) # Sine distortion with elevation...
xf2 <- xf2 / sin(rad(90 - yf2))
xf3 <- xf3 / sin(rad(90 - abs(yf3)))

xf2[yf2 > 90] <- seq(180, -180, length.out = 179)
yf2[yf2 > 90] <- 90

xf2[xf2 > 180] <- 180
xf2[xf2 < -180] <- -180

Sr1    <- Polygon(cbind(xf1, yf1))
Sr2    <- Polygon(cbind(xf2 ,yf2))
Srs1   <- Polygons(list(Sr1), "s1")
Srs2   <- Polygons(list(Sr2), "s2")
geom00 <- SpatialPolygons(list(Srs1, Srs2), 1:2, proj4string = CRS(crs.wgs))

temp   <- spTransform(geom00, crs.ste)
coordt <- temp@polygons[[2]]@Polygons[[1]]@coords

tempc <- coordt[67:16, ]

x3 <- cos(rad(0:360))*0.38e6
y3 <- sin(rad(0:360))*0.38e6

coordt[1:91, 1] <- x3[1:91]
coordt[1:91, 2] <- y3[1:91]
coordt[92:93, 1] <- c(0, x3[91])
coordt[92:93, 2] <- c(0, y3[91])
coordt[94:306, 1] <- x3[92:304]
coordt[94:306, 2] <- y3[92:304]
coordt[307:358, ] <- tempc
coordt[359, ] <- c(x3[1], y3[1])
coordt <- coordt[1:359, ]

Sr1    <- Polygon(temp@polygons[[1]]@Polygons[[1]]@coords)
Sr2    <- Polygon(coordt)
Srs1   <- Polygons(list(Sr1), "s1")
Srs2   <- Polygons(list(Sr2), "s2")
geom00 <- SpatialPolygons(list(Srs1, Srs2), 1:2, proj4string = CRS(crs.ste))
plot(geom00)
geom00 <- spTransform(geom00, crs.wgs)
coordt <- geom00@polygons[[2]]@Polygons[[1]]@coords

t1 <- coordt[1:54, ]
t1[, 1] <- t1[, 1] + 90
t2 <- cbind(seq(123, -180, -1), rep(coordt[1, 2], 304))
t3 <- cbind(c(-180, 180, 180), c(90, 90, coordt[1, 2]))
coordt <- rbind(t1, t2, t3)

temp1 <- geom00@polygons[[1]]@Polygons[[1]]@coords
temp1[, 1] <- temp1[, 1] + 90
Sr1   <- Polygon(temp1)
Sr2   <- Polygon(coordt)
xf3   <- xf3 + 90
Sr3   <- Polygon(cbind(xf3 ,yf3))
Srs1  <- Polygons(list(Sr1), "s1")
Srs2  <- Polygons(list(Sr2), "s2")
Srs3  <- Polygons(list(Sr3), "s3")
geom1 <- SpatialPolygons(list(Srs1, Srs2, Srs3), 1:3, proj4string = CRS(crs.wgs))

plot(spTransform(geom1[1:2], crs.ste), col = "black")

# Didn't work quite well in the end, when masked raster is projected it show some
# distortions where a simple line should be present, as observed for side view. Correction
# is:

geom1ste <- spTransform(geom1, crs.ste)
temp <- geom1ste

temp <- temp[2]@polygons[[1]]@Polygons[[1]]@coords
data <- as.data.frame(temp[1:2, ])
colnames(data) <- c("x", "y")
tempxa <- seq(data[1, 1], data[2, 1], length.out = 50)
tempya <- predict(lm(y~x, data), data.frame(x = tempxa))
data <- as.data.frame(temp[2:3, ])
colnames(data) <- c("x", "y")
tempxb <- seq(data[1, 1], data[2, 1], length.out = 50)
tempyb <- predict(lm(y~x, data), data.frame(x = tempxb))
data <- as.data.frame(temp[52:53, ])
colnames(data) <- c("x", "y")
tempxc <- seq(data[1, 1], data[2, 1], length.out = 50)
tempyc <- predict(lm(y~x, data), data.frame(x = tempxc))
data <- as.data.frame(temp[53:54, ])
colnames(data) <- c("x", "y")
tempxd <- seq(data[1, 1], data[2, 1], length.out = 50)
tempyd <- predict(lm(y~x, data), data.frame(x = tempxd))
temp   <- rbind(cbind(tempxa, tempya), cbind(tempxb, tempyb)[-1, ], temp[4:52, ], 
		cbind(tempxc, tempyc)[-1, ], cbind(tempxd, tempyd)[-1, ], temp[55:361, ])
rownames(temp) <- NULL
colnames(temp) <- c("x", "y")
Sr2p   <- Polygon(temp)
Srs2p  <- Polygons(list(Sr2p), "s2")
Srs4p  <- Polygons(list(Sr2p), "s4")
Srs2p.sp <- SpatialPolygons(list(Srs2p, Srs4p), 1:2, proj4string = CRS(crs.ste))
Srs2p.sp <- spTransform(Srs2p.sp, crs.wgs)
temp  <- Srs2p.sp[1]@polygons[[1]]@Polygons[[1]]@coords
temp[551:552, 1] <- c(-180,180)
Sr2   <- Polygon(temp)
Srs2  <- Polygons(list(Sr2), "s2")

geom1 <- SpatialPolygons(list(Srs1, Srs2, Srs3), 1:3, proj4string = CRS(crs.wgs))
plot(spTransform(geom1[1:2], crs.ste), col = "black")

shadow_1 <- geom1
shadow_2 <- geom2

temp <- SpatialPolygonsDataFrame(shadow_1, data.frame(name = 1:3), match.ID = F)
writeOGR(temp, "sim/shadow", "shadow1_090", driver="ESRI Shapefile", overwrite = T)

temp <- SpatialPolygonsDataFrame(shadow_2, data.frame(name = 1:3), match.ID = F)
writeOGR(temp, "sim/shadow", "shadow2_090", driver="ESRI Shapefile", overwrite = T)


#
# Plots:
#

library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
	"\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
	"\\documentclass[18pt]{article}")

tikz("plot/shadow_40_90.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	par(mar = c(1, 1, 1, 1))
	plot(zlin_ste, col = "grey45", lwd = 4, xlim = c(-18525000, 18525000), ylim = c(-18525000, 18525000))
	lines(c(0, 0), c(-12350000, 12350000), lwd = 4, col = "grey45")
	lines(c(-12350000, 12350000), c(0, 0), lwd = 4, col = "grey45")
	plot(spTransform(geom2[1:2], crs.ste), col = "black", add = TRUE)
	text(x = 9.6E6, y = -9.6E6, labels = '90$^{\\circ}$', col = 'grey45', cex = 1.5)
	text(x = 4.0E6, y = -4.0E6, labels = '40$^{\\circ}$', col = 'grey45', cex = 1.5)
	text(x = 1.95E6, y = -1.95E6, labels = '15$^{\\circ}$', col = 'grey45', cex = 1.5)
	text(x = 0.3E6, y = -13.2E6, labels = '0$^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
	text(x = -13.8E6, y = 0, labels = '270$^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
	text(x = 0.3E6, y = 13.1E6, labels = '180$^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
	text(x = 13.6E6, y = 0, labels = '90$^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
dev.off()
tools::texi2pdf("plot/shadow_40_90.tex")
file.rename("shadow_40_90.pdf", "plot/shadow_40_90.pdf")
file.remove("shadow_40_90.aux", "shadow_40_90.log")

tikz("plot/shadow_00_90.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	par(mar = c(1, 1, 1, 1))
	plot(zlin_ste, col = "grey45", lwd = 4, xlim = c(-18525000, 18525000), ylim = c(-18525000, 18525000))
	lines(c(0, 0), c(-12350000, 12350000), lwd = 4, col = "grey45")
	lines(c(-12350000, 12350000), c(0, 0), lwd = 4, col = "grey45")
	plot(spTransform(geom1[1:2], crs.ste), col = "black", add = TRUE)
	text(x = 9.6E6, y = -9.6E6, labels = '90$^{\\circ}$', col = 'grey45', cex = 1.5)
	text(x = 4.0E6, y = -4.0E6, labels = '40$^{\\circ}$', col = 'grey45', cex = 1.5)
	text(x = 1.95E6, y = -1.95E6, labels = '15$^{\\circ}$', col = 'grey45', cex = 1.5)
	text(x = 0.3E6, y = -13.2E6, labels = '0$^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
	text(x = -13.8E6, y = 0, labels = '270$^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
	text(x = 0.3E6, y = 13.1E6, labels = '180$^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
	text(x = 13.6E6, y = 0, labels = '90$^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
dev.off()
tools::texi2pdf("plot/shadow_00_90.tex")
file.rename("shadow_00_90.pdf", "plot/shadow_00_90.pdf")
file.remove("shadow_00_90.aux", "shadow_00_90.log")

