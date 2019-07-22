# Eg estimation from PTFE plaques
#
# Analysis part 7 - Evaluate the average BRDFg in the absece of tilt
#
# The plaque measurement can be represented by:
#
# Lp = BRDFg * Eg
#
# Where BRDFg is the incident irradiance distribution averaged BRDF. For a good conversion
# the plugged in value for the BRDFg must be similar to the actual BRDFg. Values plugged 
# are the BRDF(Sun, Sensor) or rho(Sensor)/pi, where rho(Sensor) is the hemispherical 
# directional reflectance at the sensor view geometry. However, based on averaged sky 
# radiance distributions for clear skies, it is in principle possible to calculate the 
# average BRDFg. If BRDFg can then be modeled by an equation with dependency on lambda 
# and Sun zenith angle, then potentially errors can be reduced, to the extent that actual
# incident irradiance distributions approach those of calculated averages. This is no 
# more approximate than the cosine effect correction for deviation of irradiance sensor 
# from ideal cosine collectors. However, this possibility for now was not further pursued.
# To make this correction more generic for personal plaques, the BRDFg can be converted to 
# a factor that scales rho to achieve that BRDFg. This value is pi for a true lambertian 
# plaque.
#
# Version: 2018-12-03
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

library(akima)
library(abind)
library(sp)
library(raster)
library(rgeos)

#
# Calculate the weighted averaged BRDFs for a Sun zenith angle of 42 degrees:
#

res_brdfg  <- array(NA, dim = c(length(lambda), 10, 2, 2, length(type)), 
		dimnames = list(paste0(lambda, 'nm'), c('OC', 'ts00', 
		paste0('ts', seq(10, 80, 10))), c('tv00', 'tv40'), c('white', 'grey'), 
		type)) 

viewg  <- rbind(c(0, 270), c(40, 270))

for(s in 1:length(type)) {
	fls <- list.files(path = "sim/skyrad", pattern = type[s], full.names = T)
	par.fls <- c("sim/skyrad/overcast_simpar.txt", grep(pattern = "_simpar", fls, value = T))
	dif.fls <- c("sim/skyrad/overcast_diffuse.txt", grep(pattern = "_diffuse", fls, value = T))
	dir.fls <- c("sim/skyrad/overcast_direct.txt", grep(pattern = "_direct", fls, value = T))
	sky.fls <- c("sim/skyrad/overcast_sky.txt", grep(pattern = "_sky", fls, value = T))
	for(sza in 1:length(sky.fls)) {
		sim_par  <- read_simpar(par.fls[sza])
		dire     <- read.table(dir.fls[sza], header = T)
		dife     <- read.table(dif.fls[sza], header = T)
		skydata  <- read.table(sky.fls[sza], header = T)
		skydata  <- as.matrix(skydata[, -c(1:2), drop = F])
		dim(skydata)   <- c(length(sim_par$zenith), length(sim_par$razimuth)+180, length(sim_par$lambda))
		for(p in 1:2) {
			if(p == 1) {
				coef <- ptfe_white_all
			} else {
				coef <- ptfe_grey_532
			}
			for(l in 1:length(lambda)) {
				sdata <- skydata[,, l]
				for(g in 1:nrow(viewg)) {
					if(p == 1) {
						scale <- ptfe_99_rho[l, g]
						brdf  <- brdf_ptfe99[[g]][,, l]
					} else {
						scale <- ptfe_10_rho[l, g]
						brdf  <- brdf_ptfe10[[g]][,, l]
					}

					tdata  <- sdata * brdf * scale
					intp   <- bicubic(x = rad(0:90), y = rad(0:360), z = tdata, x0 = gxy[, 1], y0 = gxy[, 2])$z
					intp   <- intp * fx
					dim(intp) <- c(length(sx), length(sy))
					intp   <- intp[-length(sx), ] + intp[-1, ]
					intp   <- intp[, -length(sy)] + intp[, -1]
					tempdf <- sum(intp / 4)

					dbrdf   <- zernike_g(coef, sim_par$sza, 180, viewg[g, 1], viewg[g, 2], sim_par$lambda[l])
					if(p == 1) {
						dbrdf <- dbrdf * scale
					}

					tempdr  <- dire[l, 2] * dbrdf

					res_brdfg[l, sza, g, p, s] <- (tempdr + tempdf) / (dife[l, 2] + dire[l, 2])

					print(c(type[s], sim_par$sza, sim_par$lambda[l], g, round(res_brdfg[l, sza, g, p, s], 4)))
				}
			}
		}
	}
}

save(res_brdfg, file = "res/res_brdfg.Rdata")




