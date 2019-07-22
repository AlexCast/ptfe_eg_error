#
# Plaque tilt
#
# Rotation of the radiance field to represent tilting of the plaque's normal 
# with reference to zenith.
#
# Version: 2018-07-22
#
# Alexandre Castagna
# alexandre.castagna@ugent.be
# Protistology and Aquatic Ecology
# Biology Department, Faculty of Sciences, Gent University
# Krijgslaan, 281, S8, 3rd floor, C-wing, Gent 9000, Belgium - BE
#


# Function: tilt
#
# Performs the rotation of the reference frame of the radiance field such to 
# represent the the relative change in angle between the PTFE plaque's 
# normal and the true zenith.
#
# Note that function does not intend to be generic and is dependent on the 
# specific spacial and angular references and data structure used in the 
# project.
#
# Arguments:
# ltab    - Sky radiance distribution in longtable format
# azitilt - Azimuthal angle of the normal vector to the local zenith.
# zentilt - Zenithal angle of the normal vector to the local zenith.
# sim_par - Simulation parameters
#

tilt <- function(ltab, azitilt, zentilt, sim_par) {
	azitilt   <- 180 - azitilt
	ltab      <- as.matrix(ltab[, -c(1:2), drop = F])
	dim(ltab) <- c(181, 361, ncol(ltab))
	lightfr   <- brick(ltab, xmn = -0.5, xmx = 360.5, ymn = -0.5, ymx = 180.5)

	tiltm     <- sphere.to.dir.cos(cbind(rad(zentilt), rad(0)))
	origdir   <- expand.grid(0:90, 0:360)
	origdir[, 2] <- (origdir[, 2] + azitilt) %% 360
	rotadir   <- deg(dir.cos.to.sphere(update.cosine(psi = rad(origdir[, 1]), 
                       phi = rad(origdir[, 2]), cos.m = tiltm)))
	rotadir[, 2] <- (rotadir[, 2] - azitilt) %% 360
	rotadir[, 1] <- 180 - rotadir[, 1]
	rotadir.sp   <- SpatialPoints(coords = rotadir[, 2:1])

	rtab  <- cbind(expand.grid(0:90, 0:360), extract(lightfr, rotadir.sp))
	colnames(rtab) <- c("tetha", "phi", paste0("l", sim_par$lambda))

#	write.table(format(rtab, digits = 4, scientific = T), file.path(path, fl), row.names = F, quote = F)
	return(rtab)
}

