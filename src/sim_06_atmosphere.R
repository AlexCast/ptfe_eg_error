#
# Atmospheric model
#
# Functions to setup the optical properties of the atmosphere.
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
# f0     - L66  - Gives the extraterrestrial irradiance corrected for Earth-Sun distance.
# p_ma   - L111 - Compute the density of density of moist air or its components.
# n_wv   - L156 - Compute the real part of the refractive index of water vapor.
# n_air  - L189 - Compute the real part of the refractive index of dry or moistured air.
# tau_r  - L245 - Compute the Rayleigh optical depth for 1 atmosphere.
# p_r    - L274 - Compute the Rayleigh phase function at the requested scattering angle.
# tau_o3 - L307 - Compute the ozone absorption optical thickness.
# tau_w  - L338 - Compute the water vapor absorption optical thickness.
# tau_g  - L364 - Absorption optical thickness of mixed gases.
# tau_aer_ang - L409 - Compute the aerosol optical thickness with Angstron's formulation.
# p_a_HG - L439 - Compute the phase function as described by a two term Henyey-Greenstein function.
# m_atm  - L472 - Compute the atmospheric path length for a given Sun zenith angle.
#
# NOTE: the phase functions described here follow the convention in atmospheric 
# optics where the normalization condition is:
# 1/2 * integral(pf(psi)*sin(psi)*dpsi) = 1
#

# Function: f0
#
# Gives the extraterrestrial irradiance corrected for Earth-Sun distance.
#
# Input:
# lambda  - wavelength (nm).
# doy     - day of the year (1 to 365).
# dataset - the data set of f0 at 1 AU. Options are thuillier.
#
# Output:
# A vector with f0 at the requested wavelengths.
#
# Details:
# As of this version only the data set of Thuillier is available. Extrapolation 
# outside the range is made repeating the most extreme value available, but will 
# send a warning message.
#
# Example:
# f0(lambda = seq(400, 2500, 100), doy = 1)
#
# Reference:
# Gordon, H. R., Brown, J. W., and Evans, R. H. 1988. Exact Rayleigh scattering 
#    calculations for use with the Nimbus-7 Coastal Zone Color Scanner. Appl. 
#    Opt. 27: 862-871.
# Gregg, W. W. and Karder, K. L. 1990. A simple spectral solar irradiance model 
#    for cloudless maritime atmospheres. Limnology and Oceanography 35, 8, 
#    1657-1675.
# Zibbordi, G. and Voss, K. J. 1989. Geometrical and Spectral Distribution of 
#    Sky Radiance: Comparison between Simulations and Field Measurements. Remote 
#    Sens. Environ. 27, 343-358.
#

f0_thuillier <- read.table("anc/Thuillier_f0.txt", skip = 1, header = F)

f0 <- function(lambda, doy, dataset = "thuillier") {
	if(doy < 1 | doy > 365) stop("jd must be in the range 1 to 365")
	if(dataset == "thuillier") f0 <- f0_thuillier
	if(min(lambda) < min(f0[, 1]) | max(lambda) > max(f0[, 1]))
		warning("Wavelength requested outside the limits of the data ", 
			"set. Using constant extrapolation.", call. = F)
	f0 <- approx(f0[, 1], f0[, 2], xout = lambda, rule = 2)$y
	f0 <- f0 * (1 + 1.67E-2 * cos(2 * pi * (doy - 3) / 365))^2
	return(f0)
}

# Function: p_ma
#
# Compute the density of density of moist air or its components.
#
# Input:
# p   - surface pressure (KPa).
# co2 - carbon dioxide concentration (ppm).
# h   - relative humidity (0 to 1).
# t   - temperature (ºC).
# cw  - molar fraction of water vapor in air.
# component - either NULL or "dry" or "moist".
#
# Output:
# A vector of densities.
#
# Details:
# The molar fraction of water vapor can be computed from t, p and h or can be 
# set directly, condition in which the h is ignored. The component option allow 
# to retrieve the density of the components of the moist air, either the dry air 
# component (dry) or the water vapor component (moist).
#
# Examples:
# p_ma(p = 101.325, co2 = 400, h = 0, t = 15)
# p_ma(p = 1.333, co2 = 400, h = 1, t = 20, cw = 1)
# p_ma(p = 101.325, co2 = 400, h = 0.8, t = 15, component = "dry")
# p_ma(p = 101.325, co2 = 400, h = 0.8, t = 15, component = "moist")
#
# Reference:
# Ciddor, P. E. 1996. Refractive index of air: new equations for the visible and 
#    near infrared. Applied Optics 35, 9, 1566-1573.
# Davis, R. S. 1992. Equation for the Determination of the Density of Moist Air 
#    (1981/91). Metrologia 29, 67-70.
#

p_ma <- function(p, co2, h, t, cw = NULL, component = NULL) {
	p   <- p * 1E+3
	f   <- 1.00062 + 3.14E-8 * p + 5.6E-7 * t^2
	t   <- t + 273.15
	ma  <- 1E-3 * (28.9635 + 12.011E-6 * (co2 - 400))
	svp <- exp(1.2378847E-5 * t^2 - 1.9121316E-2 * t + 33.93711047 - 
		6.3431645E+3 / t)
	if(is.null(cw)) cw  <- f * h * svp / p
	z   <- 1 - (p / t) * (1.58123E-6 - 2.9331E-8 * t + 1.1043E-10 * t^2 + 
		(5.707E-6 - 2.051E-8 * t) * cw + (1.9898E-4 - 2.376E-6 * t) * 
		cw^2) + (p / t)^2 * (1.83E-11 - 0.765E-8 * cw^2)
	p_ma <- (p * ma / (8.314510 * t * z)) * (1 - cw * (1 - 0.018015 / ma))
	if(!is.null(component)) {
		if(component == "dry") {
			p_ma_dry <- p * ma * (1 - cw)/ (8.314510 * t * z)
  			return(p_ma_dry)
		}
		if(component == "moist") {
			p_ma_moist <- p * 0.018015 * cw / (8.314510 * t * z)
			return(p_ma_moist)
		}
	}
	return(p_ma)
}

# Function: n_wv
#
# Compute the real part of the refractive index of water vapor.
#
# Input:
# lambda - wavelength (nm).
#
# Output:
# A vector refractive indexes.
#
# Example:
# n_wv(550)
#
# Reference:
# Ciddor, P. E. 1996. Refractive index of air: new equations for the visible and 
#    near infrared. Applied Optics 35, 9, 1566-1573.
# Owens, J. C. 1967. Optical Refractive Index of Air: Dependence on Pressure, 
#    Temperature and Composition. Applied Optics 6, 1, 51-59.
#

n_wv <- function(lambda) {
	lambda <- 1 / (lambda / 1E+3)
	1 + (1.022 * (295.235 + 2.6422 * lambda^2 - 0.032380 * lambda^4 + 
		0.004028 * lambda^6) * 1E-8)
}

# Function: n_air
#
# Compute the real part of the refractive index of dry or moistured air.
#
# Input:
# lambda  - wavelength (nm).
# p       - surface pressure (KPa).
# co2     - carbon dioxide concetration in (ppm).
# h       - relative humidity (0 to 1).
# t       - temperature (ºC).
#
# Output:
# A vector refractive indexes.
#
# Details:
# Ciddor data and formulation is valid for the range in the visible. Ideally, 
# refractive index of moist air in the NIR should be incorporated.
#
# Examples:
# n_air(lambda = c(400, 500, 600, 700), p = 101.325, co2 = 450, h = 0, t = 15)
# n_air(lambda = c(400, 500, 600, 700), p = 101.325, co2 = 450, h = 0.9, t = 15)
#
# Reference:
# Ciddor, P. E. 1996. Refractive index of air: new equations for the visible and 
#    near infrared. Applied Optics 35, 9, 1566-1573.
#

n_air <- function(lambda, p = 101.325, co2 = 450, h = 0, t = 15) {
	n_air <- 1 + (((5792105 / (238.0185 - (lambda / 1E+3)^-2)) + 
		(167917 / (57.362 - (lambda / 1E+3)^-2))) * 1E-8)
	if(co2 != 450)
		n_air <- 1 + ((n_air - 1) * (1 + 5.34E-7 * (co2 - 450)))
	if(h > 0) {
		n_wv  <- n_wv(lambda)
		la    <- (n_air^2 - 1) / (n_air^2 + 2)
		lw    <- (n_wv^2 - 1) / (n_wv^2 + 2)
		paxs  <- p_ma(p = 101.325, co2 = 400, h = 0, t = 15)
		pws   <- p_ma(p = 1.333, co2 = 400, h = 0, t = 20, cw = 1)
		pa    <- p_ma(p = p, co2 = co2, h = h, t = t, component = "dry")  
		pw    <- p_ma(p = p, co2 = co2, h = h, t = t, component = "moist")
		l     <- la * (pa / paxs) + lw * (pw / pws)
		n_air <- sqrt((1 + 2 * l) / (1 - l))
	}
	return(n_air)
}

# Function: tau_r 
#
# Compute the Rayleigh optical depth for 1 atmosphere.
#
# Input:
# lambda  - wavelength (nm).
# p       - surface pressure (KPa).
# co2     - carbon dioxide concentration in (ppm).
# h       - relative humidity (0 to 1).
# t       - temperature (ºC).
#
# Output:
# A vector of optical depths.
#
# Examples:
# tau_r(lambda = seq(400, 700, 100), p = 101.325, co2 = 450, h = 0, t = 15)
# tau_r(lambda = seq(400, 700, 100), p = 101.325, co2 = 450, h = 0.9, t = 15)
#
# References:
# Bodhaine, B. A., Wood, N. B., Dutton, E. G., Slusser, J. R. 1999. On Rayleigh 
#    Optical Depth Calculations. Journal of Atmospheric and Oceanic Technology 
#    16, 1854-1861.
#

f_n2 <- function(lambda) {
	1.034 + (3.17E-4 * (lambda / 1E+3)^-2)
}

f_o2 <- function(lambda) {
	1.096 + (1.385E-3 * (lambda / 1E+3)^-2) + (1.448E-4 * (lambda / 1E+3)^-4)
}

f_air <- function(lambda, co2) {
	((78.084 * f_n2(lambda)) + (20.946 * f_o2(lambda)) + (0.934) + 
		((co2 / 1E+4) * 1.150)) / (78.084 + 20.946 + 0.934 + (co2 / 1E+4))
} 

tau_r <- function(lambda, p = 101.325, co2 = 450, h = 0, t = 15) {
	ma    <- 15.0556 * (co2 / 10^6) + 28.9595
	Ns    <- 2.546899E+19
	n_air <- n_air(lambda, p, co2, h, t)
	cs_r  <- ((24 * pi^3 * (n_air^2 - 1)^2) / ((lambda * 10^-7)^4 * Ns^2 * 
		(n_air^2 + 2)^2)) * f_air(lambda, co2)
	tau_r <- cs_r * (p * 100 * (6.0221367E+23)) / (ma * 9.807)
	return(tau_r)
}

# Function: p_r
#
# Compute the Rayleigh phase function at the requested scattering angle.
#
# Input:
# psi - scattering angle (radians).
#
# Output:
# A vector of phase function at the requested directions.
#
# Example:
# p_r(rad(c(45, 135)))
#
# References:
# Zimbordi, G. and Voss, K. J. 1989. Geometrical and Spectral Distribution of 
#    Sky Radiance: Comparison between simulations and field measurements. Remote 
#    Sens. Environ. 27, 343-358.
#

p_r <- function(psi) {
	3 * (1 + cos(psi)^2) / 4
}

# Function: tau_o3
#
# Compute the ozone absorption optical thickness.
#
# Input:
# lambda - wavelength (nm).
# c_o3   - concentration of ozone (dobson).
# 
# Output:
# A vector with the optical thickness of ozone.
#
# Details:
# The ozone coefficients provided with SCIATRAN at 273 K are used in the 
# computations. The conversion used is 1 Dobson = 2.69E+16 molecules per cm^2.
#
# Example:
# tau_o3(350, c(200, 300))
#
# References:
# Iqbal, M. 1983. An Introduction to Solar Radiation, Academic, New York, pp. 
#    128-133
# Zimbordi, G. and Voss, K. J. 1989. Geometrical and Spectral Distribution of 
#    Sky Radiance: Comparison between simulations and field measurements. Remote 
#    Sens. Environ. 27, 343-358.
#

k_o3_sciamachy <- read.table("anc/SCIA_O3_Temp_cross-section_V4.1.DAT", 
	skip = 20, header = FALSE)

tau_o3 <- function(lambda, c_o3) {
	k_o3   <- approx(k_o3_sciamachy[, c(1, 5)], xout = lambda, 
			method = "linear", rule = 2)$y
	tau_o3 <- k_o3 * c_o3 * 2.69E+16
	return(tau_o3)
}

# Function: tau_w
#
# Compute the water vapor absorption optical thickness.
#
# Input:
# lambda - wavelength (nm).
# c_w    - preciptable water in the atmosphere (cm).
# 
# Output:
# A vector with the optical thickness of water vapor.
#
# Example:
# tau_w(750, 2.5)
#
# References:
# Iqbal, M. 1983. An Introduction to Solar Radiation, Academic, New York, pp. 
#    128-133.
# Zimbordi, G. and Voss, K. J. 1989. Geometrical and Spectral Distribution of 
#    Sky Radiance: Comparison between simulations and field measurements. Remote 
#    Sens. Environ. 27, 343-358.
#

k_w_leckner <- read.table("anc/Leckner_water_vapor.DAT", skip = 2, header = F)

tau_w <- function(lambda, c_w) {
	k_w <- approx(k_w_leckner, xout = lambda, rule = 2, method = "linear")$y
	tau_w <- (0.2385 * k_w * c_w) / (1 + 20.07 * k_w * c_w)^0.45
	return(tau_w)
}

# Function: tau_g
#
# Absorption optical thickness of mixed gases.
#
# Input:
# lambda - wavelength (nm).
#
# Output:
# A vector with the optical thickness of the mixed gases.
#
# Details:
# k_g, the effective absorption coefficient, is taken from Leckner (1978).
#
# Reference:
# Iqbal, M. 1983. An Introduction to Solar Radiation, Academic, New York, pp. 
#    128-133.
#

k_g_leckner <- read.table("anc/Leckner_mixed_gases.DAT", skip = 2, header = F)

tau_g <- function(lambda) { 
	k_g <- approx(k_g_leckner, xout = lambda, rule = 2, method = "linear")$y
	tau_g <- (1.41 * k_g) / (1 + 118.93 * k_g)^0.45
	return(tau_g)
}

# Function: tau_aer_ang
#
# Compute the aerosol optical thickness with Angstron's formulation.
#
# Imput:
# lambda - wavelength (nm).
# alpha  - spectral dependency coefficient.
# beta   - optical thickness at 1000 nm.
#
# Output:
# A vector with the aerosol optical thickness at the requested wavelengths.
#
# Example:
# tau_aer_ang(550, 0.7, beta(0.7, 28))
#
# Details:
# According to Iqbal 1983:
# Atmosphere	beta	alpha	Visibility (km)
# Clean 	0.0	1.30	340
# Clear 	0.1	1.30	28
# Turbid 	0.2	1.30	11
# Very  turbid  0.4	1.30	<5
#
# w0_a = 0.6 for urban environments
# w0_a = 0.9 for continental aerorol
# w0_a = 1.0 for marine
#
# Beta can be computed from visibility and an estimation of alphaaccording to 
# McClatchey and Selby.
#
# Reference:
# Iqbal, M. 1983. An Introduction to Solar Radiation, Academic, New York, pp. 
#    128-133.
#

beta <- function(alpha, vis) {
	(0.55)^alpha * (3.912 / vis - 0.01162) * (0.02472 * (vis - 5) + 1.132)
}
 
tau_aer_ang <- function(lambda, alpha, beta) {
	lambda <- lambda / 1000
	tau_aer <- beta * lambda^(-alpha)
	return(tau_aer)
}

# Function: p_a_HG
#
# Compute the phase function as described by a two term Henyey-Greenstein function.
#
# Input:
# a  - fraction of the first assymetry parameter.
# g1 - assymetry parameter 1.
# g2 - assymetry parameter 2. 
#
# Output:
# A vector with the phase function at the requested scattering angles.
#
# Details:
# None.
#
# Example:
# p_a_HG(rad(0:180), a = 0.962, g1 = 0.713, g2 = 0.759)
#
# Reference:
# Zimbordi, G. and Voss, K. J. 1989. Geometrical and Spectral Distribution of 
#    Sky Radiance: Comparison between simulations and field measurements. Remote
#    Sens. Environ. 27, 343-358.
#

p_a_HG <- function(psi, a = 0.962, g1 = 0.713, g2 = 0.759) {
	(((1 - g1^2) * a) / (1 + g1^2 - 2 * g1 * cos(psi))^1.5) + 
	(((1 - g2^2) * (1- a)) / (1 + g2^2 + 2 * g2 * cos(psi))^1.5)
}

# Function: m_atm
#
# Compute the atmospheric path length for a given Sun zenith angle. If the path 
# length is for Ozone an increased path length is required due to mainly 
# stratospheric distribution (Paltridge and Platt, as cited in Gregg and 
# Carder 1990).
# 
# Input:
# theta_s - a vector of Sun zenith angle in radians (0 to pi/2).
# oz      - logical. Should calculations be made for ozone?
#
# Output:
# A vector of mair mass.
#
# Example:
# m_atm(rad(80))
# m_atm(rad(80), oz = T)
#
# Reference:
# Gregg, W. W. and Karder, K. L. 1990. A simple spectral solar irradiance model 
#    for cloudless maritime atmospheres. Limnology and Oceanography 35, 8, 
#    1657-1675.
# Kasten, F. 1966. A new table and approximate formula for relative optical air 
#    mass. Arch. Meteorol. Geophys. Bioklimatol. B14, 206:223.
# Paltridge, G. W.; and Platt, C. M. R. 1976. Radiative processes in meteorology 
#    and climatology. Elsevier.
#

m_atm <- function(theta_s, oz = F) {
	if(oz) {
		1.0035 / sqrt(cos(theta_s)^2 + 0.007)
	} else {
		1 / (cos(theta_s) + 1.5 * (93.885 - deg(theta_s))^-1.253)
	}
}

