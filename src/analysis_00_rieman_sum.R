# Eg estimation from PTFE plaques
#
# Analysis part 0 - Define integration procedure
#
# Adaptive multidimensional integration is available in R, but in some of the 
# simulations it failed to converge, specially under tilt. Maybe I could have solved with 
# special tweaks in the control parameters, but since this integration is simple, a fixed 
# integration grid is defined, based on the two dimensional trapezoidal rule. Comparisons 
# show very good agreement between the adaptive and fixed grid procedures.
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
# Define integration parameters:
#

dx  <- (rad(90) - rad(0))  / 720							# Step in theta (1/8th of a degree)
dy  <- (rad(360) - rad(0)) / 2880							# Step in phi (1/8th of a degree)
sx  <- seq(rad(0), rad(90),  dx)							# Theta coordinate
sy  <- seq(rad(0), rad(360), dy)							# Phi coordinate
gxy <- expand.grid(sx, sy)								# Combined coordinates
fx  <- abs(cos(gxy[, 1])) * sin(gxy[, 1]) * dx * dy					# Factor to multiply interpolated matrix. It includes dx, 
											# dy and the cosine and sine weights of the solid angle integration.

#
# The correct way to compute the integral by the trapezoidal rule is shown below, tdata is 
# the sky radiance matrix at the original 1x1 degree resolution. Bicubic interpolation is 
# performed first to get values at the integration grid points. The interpolated radiance 
# values are then multiplied by the step, cosine and sine factor (fx). Data structure is 
# changed to a bi-dimensional array and the trapezoidal rule in two dimension is applied. 
# The argument z in the bicubic function is the only value to be changed: is the input 
# data.
# Example:
#
# library(akima)
# intp <- bicubic(x = rad(0:90), y = rad(0:360), z = tdata, x0 = gxy[, 1], y0 = gxy[, 2])$z
# intp <- intp * fx
# dim(intp) <- c(length(sx), length(sy))
# intp <- intp[-length(sx), ] + intp[-1, ]
# intp <- intp[, -length(sy)] + intp[, -1]
# sum(intp / 4)
#

# Tests:
#
# library(akima)
# library(cubature)
#
# par.fls <- list.files(path = "sim", pattern = "simpar", full.names = T)
# dif.fls <- list.files(path = "sim", pattern = "diffuse", full.names = T)
# dir.fls <- list.files(path = "sim", pattern = "direct", full.names = T)
# sky.fls <- list.files(path = "sim", pattern = "sky", full.names = T)
#
# i = 6
# sim_par  <- read_simpar(par.fls[i])
# dire     <- read.table(dir.fls[i], header = T)
# dife     <- read.table(dif.fls[i], header = T)
# skydata  <- read.table(sky.fls[i], header = T)
# skydata  <- as.matrix(skydata[, -c(1:2), drop = F])
# l = 21
# dim(skydata)   <- c(length(sim_par$zenith), length(sim_par$razimuth), length(sim_par$lambda))
# sdata <- skydata[,,l]
# sdata <- cbind(sdata, sdata[, 180:1])
#
# a) Using cubature:
#
# l_plaque.v <- function(x, data) {
#	lsbrdf  <- bicubic(x = rad(0:90), y = rad(0:360), z = data, x0 = x[1,], 
#                     y0 = x[2,])$z
#	x       <- rbind(lsbrdf, x)
#	res     <- matrix(apply(x, 2, function(z) z[1] * abs(cos(z[2])) * 
#                    sin(z[2])), ncol = ncol(x))
#	return(res)
# }
#
# cuba  <- hcubature(f = l_plaque.v, lowerLimit = c(0, 0), upperLimit = c(pi/2, pi*2), data = sdata, vectorInterface = TRUE)$integral
#
# b) Using fixed grid:
#
# intp  <- bicubic(x = rad(0:90), y = rad(0:360), z = sdata, x0 = gxy[, 1], y0 = gxy[, 2])$z
# intp  <- intp * fx
# dim(intp) <- c(length(sx), length(sy))
# intp  <- intp[-length(sx), ] + intp[-1, ]
# intp  <- intp[, -length(sy)] + intp[, -1]
# fgrid <- sum(intp / 4)
#
# c) Compare:
#
# ((fgrid - cuba) / cuba) * 100
#

