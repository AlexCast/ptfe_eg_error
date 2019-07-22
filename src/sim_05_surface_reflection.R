#
# Surface reflection
#
# Uses a measured water leaving reflectance for coastal waters in Long Island 
# Sound and a Fresnel reflection for a flat surface, assuming the diffuse 
# component as isotropic fraction, 20% of the total Ed.
#
# Version: 2018-05-31
#
# Alexandre Castagna
# alexandre.castagna@ugent.be
# Protistology and Aquatic Ecology
# Biology Department, Faculty of Sciences, Gent University
# Krijgslaan, 281, S8, 3rd floor, C-wing, Gent 9000, Belgium - BE
#
# Data list:
# coastal      - L27  - Coastal water leaving lambertian equivalent reflectance (unitless)
# lake         - L34  - Lake water leaving lambertian equivalent reflectance (unitless)
# n_w          - L56  - Fresh water refractive index (real part)
# n_sw         - L57  - Marine water refractive index (real part)
# fref_w_diff  - L76  - Fresnel reflectance for fresh water under uniform illumination
# fref_Sw_diff - L77  - Fresnel reflectance for marine water under cardioidal illumination
#
# Functions list:
# surfacelight - L174 - Compute the (flat) surface radiance distribution for a 
#                       given sky radiance distribution. 
#
# Read water-leaving reflectances cases:
#

coastal <- read.table("anc/wl_reflectance_coastal.dat", header = T)
coastal[, 2][coastal[, 2] < 0] <- 0
temp <- as.numeric(filter(coastal[, 2], rep(1/7, 7), "convolution"))
coastal[4:606, 2] <- temp[4:606]
coastal[, 2] <- coastal[, 2] * pi
coastal[650:651, 2] <- 0

lake <- read.table("anc/wl_reflectance_lake.dat", header = T)
lake[, 2][lake[, 2] < 0] <- 0
temp <- as.numeric(filter(lake[, 2], rep(1/7, 7), "convolution"))
lake[4:606, 2] <- temp[4:606]
lake[, 2] <- lake[, 2] * pi
lake[645:651, 2] <- 0

png("plot/rho_wl.png", width = 1000, height = 1000, res = 200)
par(mar = c(5, 5, 3, 2))
xlab = expression(lambda~~(nm))
ylab = expression(italic(rho[wl])~~(unitless))
plot(coastal, ylim = c(0, 0.035), type = "l", xlab = xlab, ylab = ylab)
lines(lake, lty = 2)
legend("topright", c("Coastal", "Lake"), lty = 1:2, bty = "n")
dev.off()

#
# Read refractive indexes generated with WOPP:
#

n_w  <- read.table("anc/nw_15deg_00sal.dat", skip = 2)[, 1:2]
n_sw <- read.table("anc/nw_15deg_34sal.dat", skip = 2)[, 1:2]

png("plot/n_w.png", width = 1000, height = 1000, res = 200)
par(mar = c(5, 5, 3, 2))
xlab = expression(lambda~~(nm))
ylab = expression(italic(n[w])~~(unitless))
plot(n_sw, ylim = c(1.26, 1.37), type = "l", xlab = xlab, ylab = ylab)
lines(n_w, lty = 2)
legend("topright", c("15ºC, 34 PSU", "15ºC, 00 PSU"), lty = 1:2, bty = "n")
dev.off()

#
# Compute surface reflectance for diffuse illumination:
#

# Uniform light:
fref_w_diff  <- rep(NA, nrow(n_w))
fref_sw_diff <- rep(NA, nrow(n_sw))

fun.w <- function(x) {
	sin(2 * x)
}

ifw <- integrate(fun.w, lower = 0, upper = pi / 2)$value

fun <- function(x, i) {
	refr_w <- refr.snell(teta.air = x, n_air = 1, n_water = n_w[i, 2])
	refl_w <- refl.fresnel(teta.i = x, teta.r = refr_w, n_air = 1, n_water = n_w[i, 2])
	refl_w <- refl_w * fun.w(x)
}

for(i in 1:nrow(n_w)) {
	fref_w_diff[i] <- integrate(fun, lower = 0, upper = pi / 2, i = i)$value / ifw
}

fun <- function(x, i) {
	refr_sw <- refr.snell(teta.air = x, n_air = 1, n_water = n_sw[i, 2])
	refl_sw <- refl.fresnel(teta.i = x, teta.r = refr_sw, n_air = 1, n_water = n_sw[i, 2])
	refl_sw <- refl_sw * fun.w(x)
}

for(i in 1:nrow(n_sw)) {
	fref_sw_diff[i] <- integrate(fun, lower = 0, upper = pi / 2, i = i)$value / ifw
}


# Cardioidal light:
fref_w_card  <- rep(NA, nrow(n_w))
fref_sw_card <- rep(NA, nrow(n_sw))

fun.w <- function(x) {
	(3 / (7 * pi)) * (1 + 2 * cos(x)) * sin(2 * x)
}

ifw <- integrate(fun.w, lower = 0, upper = pi / 2)$value

fun <- function(x, i) {
	refr_w <- refr.snell(teta.air = x, n_air = 1, n_water = n_w[i, 2])
	refl_w <- refl.fresnel(teta.i = x, teta.r = refr_w, n_air = 1, n_water = n_w[i, 2])
	refl_w <- refl_w * fun.w(x)
}

for(i in 1:nrow(n_w)) {
	fref_w_card[i] <- integrate(fun, lower = 0, upper = pi / 2, i = i)$value / ifw
}

fun <- function(x, i) {
	refr_sw <- refr.snell(teta.air = x, n_air = 1, n_water = n_sw[i, 2])
	refl_sw <- refl.fresnel(teta.i = x, teta.r = refr_sw, n_air = 1, n_water = n_sw[i, 2])
	refl_sw <- refl_sw * fun.w(x)
}

for(i in 1:nrow(n_sw)) {
	fref_sw_card[i] <- integrate(fun, lower = 0, upper = pi / 2, i = i)$value / ifw
}

png("plot/rho_surf_diff.png", width = 1000, height = 1000, res = 200)
par(mar = c(5, 5, 3, 2))
xlab = expression(lambda~~(nm))
ylab = expression(italic(rho[sky])~~(unitless))
plot(n_w[, 1], fref_sw_diff, type = "l", ylim = c(0.040, 0.07), xlab = xlab, ylab = ylab, 
     col = "grey")
lines(n_sw[, 1], fref_w_diff, lty = 2, col = "grey")
lines(n_sw[, 1], fref_sw_card, lty = 1)
lines(n_sw[, 1], fref_w_card, lty = 2)
legend("bottomleft", c("Cardioidal coastal", "Cardioidal lake", "Diffuse coastal", "Diffuse lake"), 
       lty = c(1, 2, 1, 2), col = c("black", "black", "grey", "grey"), bty = "n", cex = 0.8)
dev.off()

# Function: surfacelight
#
# Compute the (flat) surface radiance distribution for a given sky radiance 
# distribution. 
#
# Input:
# skydata - Sky radiance matrix as the one produced by the function skyrad.
# lambda  - Wavelength of choice. Must match exactly a wavelength in the skydata 
#           data. Only one lambda might be requested at a time.
# rhow    - water-leaving irradiance reflectance.
# eg      - global downwelling irradiance.
# simpar  - Simulation description data.
#
# Output:
# surrad - A long matrix format (same structure as skydata) where the first two 
# columns are the angle to the zenith and the relative azimuth.
#
# Details:
# Fresnel equation for unpolarized light is used to compute reflectance factors 
# for directions of the incoming radiance and those factors are applied before 
# flipping the matrix to represent specular directions. The factor only change 
# depending on the refractive indexes of air and water. In this version, the 
# refractive index of the water is fixed to a temperature of 15ºC and a salinity 
# of 34 PSU as computed from the WOPP - when the WOPP is implemented, n_water 
# can be computed for any condition of choice.
#

surfacelight <- function(skydata, lambda, rhow, eg, sim_par) {
	if(length(lambda) > 1)
		stop("The parameter lambda must be of length 1")
	id          <- which(lambda == sim_par$lambda) + 2
	lightdir    <- skydata[, 1:2]
	if(length(id) == 0)
		stop("Wavelength requested is not provided in the sky radiance data")
	sim_n_air   <- n_air(lambda, sim_par$p, sim_par$co2, sim_par$h, sim_par$t)
	sim_n_water <- approx(n_sw[, 1], n_sw[, 2], xout = lambda)$y
	sim_snell   <- refr.snell(teta.air = rad(0:90), teta.water = NA, 
			sim_n_air, sim_n_water)
	sim_fresnel <- ((sim_n_water - sim_n_air) / (sim_n_air + sim_n_water))^2
	sim_fresnel <- c(sim_fresnel, refl.fresnel(rad(1:90), sim_snell[-1], 
			sim_n_air, sim_n_water))
	surflight   <- (matrix(skydata[, id], ncol = 361) * sim_fresnel)[91:1, ]
	surflight   <- surflight + ((eg[id-2] * rhow[id-2]) / pi)
	surrad      <- cbind(skydata[, 1:2], as.vector(surflight))
	return(surrad)
}

