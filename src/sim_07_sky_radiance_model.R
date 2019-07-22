#
# Sky radiance model
#
# Contain functions to simulate the direct irradiance and the sky radiance 
# distribution according to an approximate solution to the RT equation.
#
# Version: 2018-05-31
#
# Alexandre Castagna
# alexandre.castagna@ugent.be
# Protistology and Aquatic Ecology
# Biology Department, Faculty of Sciences, Gent University
# Krijgslaan, 281, S8, 3rd floor, C-wing, Gent 9000, Belgium - BE
#
# Functions list:
# simpar      - L56  - Creates a data frame with the parameters for the simulation.
# read_simpar - L97  - Will read the parameters from a file.
# Lsky        - L139 - Sky radiance at a requested angle.
# sigma       - L183 - Auxiliary functions for Lsky.
# skyrad      - L242 - Wrapper to use Lsky for several angles.
# es_m        - L329 - Direct planar downwelling irradiance component at surface.
#

# Function: simpar
#
# Wrapper to create a data frame with the necessary components to control the 
# simulation. The simulation function and other post-processing functions use a 
# data frame that conveniently keep all the input parameters in a single 
# structure. The usefulness of this function is that it has default values that 
# could be used for many simulations.
#
# Input:
# lambda   - desired wavelengths (nm).
# c_w      - water vapor concentration (cm).
# c_o3     - ozone concentration (dobson).
# co2      - Carbon dioxide concentration (ppm).
# p        - surface pressure (KPa).
# h        - relative humidity (0 to 1).
# t        - temperature (ºC).
# tau_aer  - aerosol optical thickness at the desired wavelengths.
# w0_a     - single scattering albedo of the aerosol at the desired wavelengths.
# sza      - sun zenith angle of the simulation (º).
# rho      - surface reflectance for the desired wavelegths.
# f_0      - extraterrestrial solar irradiance for the desired wavelengths.
# zenith   - zeniths of the grid to be computed (º).
# razimuth - relative azimuths of the grid to be computed (º).
# p_a      - aerosol phase function. Should be a matrix with the first column 
#            the scattering angle and the scattering intensity in the second 
#            column.
#
# Output:
# A named list with the elements used as input.
#
# Details:
# The function will send an error if the wavelength dependent paramenters do not 
# have the same length. The only exception is the phase function; as of this 
# version, only one phase function can be supplied (spectrally independent).
#
# Examples:
# simpar()
# simpar(sza = 60)
#

simpar <- function(lambda   = 550,
                   c_w      = 2.5, 
                   c_o3     = 300, 
                   co2      = 450, 
                   p        = 101.25, 
                   h        = 0.95, 
                   t        = 15, 
                   tau_aer  = 0.15, 
                   w0_a     = 0.9, 
                   sza      = 30, 
                   rho      = 0.066, 
                   f_0      = 1, 
                   zenith   = 0:90, 
                   razimuth = 0:180, 
                   p_a      = NULL) {

	if(length(lambda) != length(tau_aer) | length(lambda) != length(w0_a) | 
	   length(lambda) != length(rho) | length(lambda) != length(f_0)) 
		stop("Wavelength dependent parameters must have the same ", 
			"length", call. = F)

	if((length(c_w) + length(c_o3) + length(co2) + length(p) + length(h) + 
	    length(t) + length(sza)) > 7)
 		stop("Except for zenith and razimuth, all wavelength ", 
			"independent parameters must be of length 1", call. = F)

	if(sza > 90) stop("Sun zenith angle must be between 0 and 90", call. = F) 

	if(is.null(p_a))
		p_a <- matrix(c(rad(0:180), p_a_HG(rad(0:180), a = 0.962, 
			g1 = 0.713, g2 = 0.759)), ncol = 2, dimnames = list(NULL, 
			c("angle", "phase_function")))

	simpar <- list(lambda, c_w, c_o3, co2, p, h, t, tau_aer, w0_a, sza, rho, 
			f_0, zenith, razimuth, p_a)
	names(simpar) <- c("lambda", "c_w", "c_o3", "co2", "p", "h", "t", 
			"tau_aer", "w0_a", "sza", "rho", "f_0", "zenith", 
			"razimuth", "p_a")
	return(simpar)
}

read_simpar <- function(fl) {
	siml <- lapply(strsplit(readLines(fl)[seq(2, 28, 2)], ","), as.numeric)
	simt <- read.table(fl, skip = 29)
	simp <- simpar(siml[[1]], siml[[2]], siml[[3]], siml[[4]], siml[[5]], siml[[6]], 
		siml[[7]], siml[[8]], siml[[9]], siml[[10]], siml[[11]], siml[[12]], 
		siml[[13]], siml[[14]], simt)
	return(simp)
}

# Function: Lsky
#
# Compute the directional distribution of sky radiance.
#
# Input:
# f_0     - extraterrestrial solar irradiance for the desired wavelengths.
# sigma   - atmospheric transmittance for the desired direction, as computed 
#           from the sigma function.
# theta_s - sun zenith angle (radians).
# tau_a   - total absorption optical thickness of the atmosphere.
#
# Output:
# Radiance from the requested direction.
#
# Details:
# This function is not intended to be called directly by the user, but from the 
# function skyrad. It is based on the model of Zimbordi and Voss 1989, but the 
# treatment of the airmass is updated with the equation of Kasten 1966. The 
# absorbers are treated as confined to the top layer of the atmosphere and the 
# atmosphere modeled as pure scattering. Is an approximation and should not be 
# used for large optical depth (> 0.2).
#
# Example:
# None.
#
# References:
# Kasten, F. 1966. A new table and approximate formula for relative optical air 
#    mass. Arch. Meteorol. Geophys. Bioklimatol. B14, 206:223.
# Zimbordi, G. and Voss, K. J. 1989. Geometrical and Spectral Distribution of 
#    Sky Radiance: Comparison between simulations and field measurements. Remote 
#    Sens. Environ. 27, 343-358.
#

Lsky <- function(f_0, sigma, theta_s, tau_a) {
  (f_0 / pi) * sigma * cos(theta_s) * exp(-tau_a * m_atm(theta_s))
 }

# Function: sigma
#
# Compute the transmittance of the atmosphere in the desired direction.
#
# Input:
# theta   - desired direction in the sky (radians).
# theta_s - sun zenith angle (radians).
# tau_d   - total scattering optical thickness of the atmosphere.
# rho     - irradiance reflectance of the surface (assumed to be lambertian).
# p1      - first term in the Legendre polynomial expansion of the global phase 
#           function.
# psi     - angle between the sun vector and the desired direction (radians).
# tau_r   - Rayleigh optical thickness.
# tau_as  - aerosol scattering optical thickness.
# p_a_fun - interpolation function for the aerosol phase function.
#
# Output:
# A vector with transmittance at the requested direction.
#
# Details:
# This function is not intended to be called directly by the user, but from the 
# function skyrad. It is based on the model of Sobolev 1975, but the treatment 
# of the airmass is updated with the equation of Kasten 1966.
#
# Example:
# None.
#
# References:
# Kasten, F. 1966. A new table and approximate formula for relative optical air 
#    mass. Arch. Meteorol. Geophys. Bioklimatol. B14, 206:223.
# Sobolev, V. V. 1975. Light scattering in planetary atmospheres. Pergamon 
#    Press, 256 pp.
# Zimbordi, G. and Voss, K. J. 1989. Geometrical and Spectral Distribution of 
#    Sky Radiance: Comparison between simulations and field measurements. Remote 
#    Sens. Environ. 27, 343-358.
#

# NOTE: In updating the airmass, is not clear if the denominator of function 
# sigma1 should also be changed. Tests I made show that it should stay as is.

sigma <- function(theta, theta_s, tau_d, rho, p1, psi, tau_r, tau_as, p_a_fun) {
 (((((1 - rho) * r(theta, tau_d)) + 2 * rho) * r(theta_s, tau_d)) / (4 + ((3 - p1) * (1 - rho) * tau_d))) - 
 (0.5 * (exp(-tau_d * m_atm(theta)) + exp(-tau_d * m_atm(theta_s)))) +
 (p_psi(psi, tau_r, tau_as, p_a_fun) - ((3 + p1) * cos(theta) * cos(theta_s))) * sigma_1(theta, theta_s, tau_d)
}

r <- function(theta, tau_d) {
 1 + ((3 / 2) * cos(theta)) + ((1 - ((3 / 2) * cos(theta))) * exp(-tau_d * m_atm(theta)))
}

sigma_1 <- function(theta, theta_s, tau_d) {
 0.25 * (exp(-tau_d * m_atm(theta)) - exp(-tau_d * m_atm(theta_s))) / (cos(theta) - cos(theta_s))
}
 
p_psi <- function(x, tau_r, tau_as, p_a_fun) { # Global phase function
 (tau_r * p_r(x) + tau_as * p_a_fun(x)) / (tau_r + tau_as)
}

p_psi_int <- function(x, tau_r, tau_as, p_a_fun) { # Wrapper of global phase function for integration
 p_psi(x, tau_r, tau_as, p_a_fun) * sin(x) * cos(x)
}

# Function: skyrad
#
# Compute the sky radiance distribution to all direction specified in a grid.
#
# Input:
# simpar - a named list with the control parameters for the simulation. See 
#          simpar function.
#
# Output:
# A three dimensional array with sky radiance distribution in the following 
# format: zenith, razimuth, wavelength.
#
# Details:
# This is a wrapper function that combine all the necessary ancillary 
# computations to calculate the sky radiance distribution. The fundamental 
# functions Lsky and sigma are based on Zimbordi and Voss 1989 and Sobolev 1975, 
# respectively. But the treatment of the airmass is updated with the equation 
# of Kasten 1966. In this model, the top layer of the atmosphere is purely 
# absorbing and confine all absorbing components while the rest of teh 
# atmosphere is treated as pure scattering. The surface is treated as a 
# lambertian reflector. Is an approximated solution of the radiative transfer 
# equation, and should not be used for large optical depths (> 0.2).
#
# Example:
# skyrad(simpar())
# skyrad(simpar(sza = 60))
#
# References:
# Kasten, F. 1966. A new table and approximate formula for relative optical air 
#    mass. Arch. Meteorol. Geophys. Bioklimatol. B14, 206:223.
# Sobolev, V. V. 1975. Light scattering in planetary atmospheres. Pergamon 
#    Press, 256 pp.
# Zimbordi, G. and Voss, K. J. 1989. Geometrical and Spectral Distribution of 
#    Sky Radiance: Comparison between simulations and field measurements. Remote 
#    Sens. Environ. 27, 343-358.
#

skyrad <- function(simpar) {
	theta_s     <- rad(simpar$sza)
	theta       <- rad(simpar$zenith)
	rphi        <- rad(simpar$razimuth)
	sim_tau_r   <- tau_r(simpar$lambda, simpar$p, simpar$co2, simpar$h, simpar$t)
	sim_tau_o3  <- tau_o3(simpar$lambda, simpar$c_o3)
	sim_tau_w   <- tau_w(simpar$lambda, simpar$c_w)
	sim_tau_g   <- tau_g(simpar$lambda)
	sim_tau_aer <- simpar$tau_aer
	sim_tau_as  <- sim_tau_aer * simpar$w0_a
	sim_tau_aa  <- sim_tau_aer * (1 - simpar$w0_a)
	sim_tau_d   <- sim_tau_r + sim_tau_as
	sim_tau_a   <- sim_tau_o3 + sim_tau_w + sim_tau_g + sim_tau_aa
	p_a_fun     <- approxfun(simpar$p_a[, 1], simpar$p_a[, 2])

	result <- array(NA, dim = c(length(simpar$zenith), length(simpar$razimuth), 
		length(simpar$lambda)))
	dimnames(result) <- list(simpar$zenith, simpar$razimuth, simpar$lambda)

	for(l in 1:length(simpar$lambda)) {
		cat(paste0("Computing at ", simpar$lambda[l], "..."))
		p1 <- integrate(p_psi_int, lower = 0, upper = pi, tau_r = sim_tau_r[l], 
                  tau_as = sim_tau_as[l], p_a_fun = p_a_fun, subdivisions = 1000L, 
                  stop.on.error = F)$value * (3 / 2)
		for(z in 1:length(simpar$zenith)) {
			for(f in 1:length(simpar$razimuth)) {
				sim_psi   <- psi(theta[z], theta_s, rphi[f])
				sim_sigma <- sigma(theta[z], theta_s, sim_tau_d[l], 
					simpar$rho[l], p1, sim_psi, sim_tau_r[l], 
					sim_tau_as[l], p_a_fun)
				result[z, f, l] <- Lsky(simpar$f_0[l], sim_sigma, 
					theta_s, sim_tau_a[l])
			}
  		}
		cat(" done!\n")
	}

	# Interpolate to fill NA if there is a theta == theta_s:
	id <- which(simpar$zenith == simpar$sza)
	if(length(id) > 0) {
		library(akima)
		xy <- expand.grid(simpar$zenith[id], simpar$razimuth)
		for(l in 1:length(simpar$lambda)) {
			result[id,,l] <- bicubic(x = simpar$zenith[-id], 
				y = simpar$razimuth, z = result[-id,,l], 
				x0 = xy[, 1], y0 = xy[, 2])$z
		}
	} 
	return(result)  
}

# Function: es_m
#
# Compute the direct component (Es) of the global irradiance (Eg).
#
# Input:
# lambda  - wavelength in nm (350 to 2500 nm).
# c_w     - water vapor concentration in cm.
# c_o3    - ozone concentration in dobson units.
# co2     - carbon dioxide concentration in ppm.
# p       - pressure in KPa.
# tau_aer - a function of lambda (see details).
# sza     - sun zenith angle in degrees (0 to 90).
# f0      - extraterrestrial irradiance.
#
# Output:
# A vector of the direct component of the global downweling irradiance.
#
# Details:
# Note that only one parameter should have length > 1, except for lambda and F0, 
# that can have the same length. tau_aer should be a function of lambda, 
# returning a vector.
#
# Examples:
# tf0   <- f0(lambda = seq(400, 2500, 100), doy = 1)
# taua <- approxfun(continental[[1]][, 1], continental[[1]][, 6])
# es_m(lambda = seq(400, 2500, 25), tau_aer = taua, sza = 30, f0 = tf0)
#
# Reference:
# Bird, R. E. and Riordan, C. 1986. Simple solar spectral model for direct and 
#    diffuse irradiance on horizontal and tilted planes at the earth’s surface 
#    for cloudless atmospheres. J. Climatol. Appl. Meteorol. 25, 87-97.
# Gregg, W. W. and Karder, K. L. 1990. A simple spectral solar irradiance model 
#    for cloudless maritime atmospheres. Limnology and Oceanography 35, 8, 
#    1657-1675.
#

es_m <- function(lambda, tau_aer, sza, f0, 
                 c_w  = 2.5,
                 c_o3 = 300, 
                 co2  = 450, 
                 p    = 101.325, 
                 h    = 0.8,  
                 t    = 15) {

	theta_s     <- rad(sza)
	sim_tau_r   <- tau_r(lambda, p, co2, h, t)
	sim_tau_o3  <- tau_o3(lambda, c_o3)
	sim_tau_w   <- tau_w(lambda, c_w)
	sim_tau_g   <- tau_g(lambda)
	sim_tau_aer <- tau_aer
	sim_tau     <- (sim_tau_r + sim_tau_w + sim_tau_g + sim_tau_aer) * 
			m_atm(theta_s) + sim_tau_o3 * m_atm(theta_s, oz = TRUE)
	es <- f0 * cos(theta_s) * exp(-sim_tau)
	return(es)
}


