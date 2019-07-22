#
# Auxiliary geometrical functions.
#


#
# Scattering angle:
#

psi <- function(theta, theta_s, rphi) {
  acos((cos(theta) * cos(theta_s)) - (sin(theta) * sin(theta_s) * cos(rphi)))
}

#
# Convert degrees to radians:
#

rad <- function(degree) {
    rad <- degree*pi/180
    return(rad)
}

#
# Convert radians to degrees:
#

deg <- function(rad) {
    deg <- rad*180/pi
    return(deg)
}

#
# Simplified Snell's refraction formulation, as necessary for this simulations.
#
# Arguments:
# theta.air = Zenith angle of incidence in radians	0 <= theta.air <= pi/2
# n_air = Refractive index of air			
# n_water = Refractive index of seawater
#
# Result:
# The refraction angle in radians.
#
# Validated = TRUE

refr.snell <- function(theta.air = NA, n_air, n_water) {
	refr.ratio <- n_water / n_air
	Refr <- asin(sin(theta.air) / refr.ratio)
	return(Refr)
}

# Simplified Fresnel's reflection formulation, as necessary for this simulations. 
#
# Arguments:
# theta.i = Zenith angle of incidence in radians	0 <= theta.i <= pi/2
# theta.r = Nadir angle of refraction in radians	0 <= theta.r <= pi/2
# n_air = Refractive index of air			
# n_water = Refractive index of seawater
#
# Result:
# The unpolarized reflection.
#
# Validated = TRUE

refl.fresnel <- function(theta.i, theta.r, n_air, n_water) {
	Ref <- 1/2*((sin(theta.i-theta.r)/sin(theta.i+theta.r))^2)+
	1/2*((tan(theta.i-theta.r)/tan(theta.i+theta.r))^2)
	id <- which(theta.i == 0)
	if(length(id) > 0) Ref[id] <- ((n_water-n_air)/(n_air+n_water))^2
	return(Ref)
}


##############################################################################
# Update vector component cosines for tilt condition
#
# This function is based on Prahl et al. 1989. If mu.z is close to the axis 
# (i.e., theta ~0 or ~180º), than the new cosines are defined based only on 
# tilt angles because mu.x and mu.y of the incident bean will be zero.
#
# Arguments:
# psi = polar tilt angle in radians		0 <= psi <= pi
# phi = azimuthal tilt angle in radians		0 <= phi <= 2pi
# cos.matrix = matrix with directional cosines, on the following order: mu.x, 
# mu.y, mu.z.
#
# Default values:
# None.
#
# Result:
# A matrix with the new directional cosines in the order: mu.x, mu.y, 
# mu.z.
#
# Validated: TRUE
##############################################################################

update.cosine <- function(psi, phi, cos.m) {
	mu.s <- cos(psi)
	sin.theta <- sin(acos(cos.m[3]))
	sin.s <- sin(psi)
	matrix.b <- rbind(sin.s*cos(phi), sin.s*sin(phi), mu.s)
	if(abs(sin.theta) > 1e-12) {
		comp.matrix.a <- c(cos.m[1]*cos.m[3]/sin.theta, cos.m[2]*cos.m[3]/sin.theta, -sin.theta, 
        			-cos.m[2]/sin.theta, cos.m[1]/sin.theta, 0, cos.m[1], cos.m[2], cos.m[3])
		matrix.a <- matrix(comp.matrix.a, 3, 3)
		new.cosines <- apply(matrix.b, 2, FUN = function(x) { matrix.a %*% matrix(x, ncol = 1) })
	} else {
		matrix.b[3, ] <- sign(cos.m[3]) * matrix.b[3, ]
		new.cosines <- matrix.b
	}
	return(new.cosines)
}

##############################################################################
# Calculate the spherical coordinates from directional cosines
#
# If theta is parallel to z axis, than phi is set at random. Note that the 
# reference system is set to nadir angle, and therefore upwelling photon 
# packets will have theta > 90.
#
# Arguments:
# cos.matrix = matrix with directional cosines, on the following order: mu.x, 
# mu.y, mu.z.
#
# Default values:
# None.
#
# Result:
# A matrix with the spherical coordinates: theta, phi.
#
# Validated: TRUE
##############################################################################


dir.cos.to.sphere <- function(cos.m) {
	theta <- acos(cos.m[3, ])
	phi   <- rep(NA, length(theta))
	id1 <- which(1 - abs(cos.m[3, ]) < 1E-12)
	if(length(id1) > 0) phi[id1] <- 2 * pi * runif(length(id1))
	id1 <- setdiff(1:length(theta), id1)
	sin.t <- sin(theta[id1])
	id2 <- which(cos.m[2, id1] >= 0)
	if(length(id2) > 0) phi[id1][id2] <- acos(round(cos.m[1, id1, drop = F][, id2] / sin.t[id2], 
                                              digits = 12))
	id2 <- setdiff(1:length(id1), id2)
	if(length(id2) > 0) phi[id1][id2] <- (2 * pi) - acos(round(cos.m[1, id1, drop = F][, id2] / sin.t[id2], 
                                      digits = 12))
	return(cbind(theta, phi))
}

##############################################################################
# Calculate directional cosines from the spherical coordinates
#
# Arguments:
# polar = a matrix with polar coordinates on the order theta, phi.
#
# Default values:
# None.
#
# Result:
# A matrix with the directional cosines: mu.x, mu.y, mu.z.
#
# Validated: TRUE
##############################################################################

sphere.to.dir.cos <- function(polar) {
	mu.x <- sin(polar[1]) * cos(polar[2])
	mu.y <- sin(polar[1]) * sin(polar[2])
	mu.z <- cos(polar[1])
	cos.m <- rbind(mu.x, mu.y, mu.z)
	return(cos.m)
}


