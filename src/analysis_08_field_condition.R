# Eg estimation from PTFE plaques
#
# Analysis part 8 - Evaluate the effect of cross-calibration with white plaque
#
# The grey plaque used in the field experiments of the UCONN cruise is a aged plaque, 
# used many times in the field and that has already been resurfaced at least once. BRDF
# and hemispherical-directional reflectance data were not available and plaque was cross
# calibrated with a pristine white (99%) reference plaque. The observed errors in the 
# field campaign are effectively lower than expected from the simulations for a grey 
# plaque using the lambertian assumption. This code investigates if the cross-calibration 
# in fact can improve performance of the grey plaque.
#
# The mechanism behind the process is that if the directional irradiance distribution is 
# fixed between the measurements of both plaques, than the ratios of the observed 
# radiances is equal to the ratios of the weighted average BRDFs - with the weighting 
# given by the incident directional irradiance distribution.
#
# Cross-calibration at nadir was made at 40 degrees SZA. Cross-calibration at 40 degrees
# view made at the same date and time failed and data was reconstructed using a cross 
# validation station. It was repeated at a 24 degrees SZA. The processing the field data 
# used the appropriate recalculated reflectance for each view angle.
#
# Version: 2018-12-03
#
# Alexandre Castagna
# alexandre.castagna@ugent.be
# Protistology and Aquatic Ecology
# Biology Department, Faculty of Sciences, Gent University
# Krijgslaan, 281, S8, 3rd floor, C-wing, Gent 9000, Belgium - BE

#
# Load libraries:
#

library(akima)
library(abind)
library(sp)
library(raster)
library(rgeos)

#
# Check field conditions:
#

s   <- 1:24
recalc <- cbind(
ptfe_99_rho[, 1] * apply(res_brdfg[, 6, 1, 2, s], 1, mean) / apply(res_brdfg[, 6, 1, 1, s], 1, mean), 
ptfe_99_rho[, 1] * apply(res_brdfg[, 4:5, 2, 2, s], 1, mean) / apply(res_brdfg[, 4:5, 2, 1, s], 1, mean)) # Because only 8/h was available for white plaque in the field, 
													  # although values at 40 degrees have only a negligible difference.
#
# Calculate new error without tilt:
#
# Simulated errors when cross-calibration values are used with the lambertian assumption 
# and for the sky conditions observed in the field (20 to 30 degrees).
#

viewg <- rbind(c(0, 270), c(40, 270))
shdwl <- list(shadow_1, shadow_2)

dimnms <- list(lambda, c('ts20', 'ts30'), c('tv00', 'tv40'), ptfes[2], type)

crosscal_shdw_tot_pi <- array(NA, dim = c(length(lambda), 2, nrow(viewg), 1, length(type)), dimnames = dimnms)

for(s in 1:length(type)) {
	fls <- list.files(path = "sim/skyrad", pattern = type[s], full.names = T)
	par.fls <- c("sim/skyrad/overcast_simpar.txt", grep(pattern = "_simpar", fls, value = T))
	dif.fls <- c("sim/skyrad/overcast_diffuse.txt", grep(pattern = "_diffuse", fls, value = T))
	dir.fls <- c("sim/skyrad/overcast_direct.txt", grep(pattern = "_direct", fls, value = T))
	sky.fls <- c("sim/skyrad/overcast_sky.txt", grep(pattern = "_sky", fls, value = T))
	for(sza in 4:5) {
		sim_par  <- read_simpar(par.fls[sza])
		dire     <- read.table(dir.fls[sza], header = T)
		dife     <- read.table(dif.fls[sza], header = T)
		skydata  <- read.table(sky.fls[sza], header = T)
		skydata  <- as.matrix(skydata[, -c(1:2), drop = F])
		dim(skydata)   <- c(length(sim_par$zenith), length(sim_par$razimuth)+180, length(sim_par$lambda))
		sun   <- SpatialPoints(cbind(0, 90 - sim_par$sza), proj4string = CRS(crs.wgs))
		for(p in 2) {
			coef <- ptfe_grey_532
			for(l in 1:length(lambda)) {
				sdata <- skydata[,,l]
				for(g in 1:nrow(viewg)) {
					scale <- ptfe_10_rho[l, g]
					brdf  <- brdf_ptfe10[[g]][,,l]
					crosscal <- recalc[l, g]
	
					shadow  <- shdwl[[g]]
					temp    <- mask(raster(sdata, xmn=-180.5, xmx=180.5, ymn=-0.5, ymx=90.5, crs = crs.wgs), shadow, inverse = T)
					tsdata  <- matrix(values(temp), ncol = 361, byrow = T)
					tsdata[is.na(tsdata)] <- 0

					tdata  <- tsdata * brdf * scale
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

					sunshdw <- gIntersects(shadow, sun)
					if(sunshdw) {
						tempdr <- 0
					} else {
						tempdr <- dire[l, 2] * dbrdf
					}

					meas_pi <- (c(tempdf, tempdr) * pi / crosscal)
					real    <- dife[l, 2] + dire[l, 2]

					crosscal_shdw_tot_pi[l, sza-3, g, 1, s] <- sum(meas_pi) / real
					print(c(type[s], sim_par$sza, sim_par$lambda[l], g, round(crosscal_shdw_tot_pi[l, sza-3, g, 1, s], 4)))
				}
			}
		}
	}
}

save(crosscal_shdw_tot_pi, file = "res/crosscal_shdw_tot_pi.Rdata")

#
# Calculate new error with tilt:
#

zets <- c(3, 6, 9, 12)
azts <- c(0, 45, 90, 135, 180, 225, 270, 315)

viewg <- rbind(c(0, 270), c(40, 270))
ptfes <- c("ptfe99", "ptfe10")

dimnms <- list(c('znt03', 'znt06', 'znt09', 'znt12'), c('azt000', 'azt045', 'azt090', 
		'azt135', 'azt180', 'azt225', 'azt270', 'azt315'), lambda, c('ts20', 
		'ts30'), c('tv00', 'tv40'), ptfes[2], type)

crosscal_all_tot_pi <- array(NA, dim = c(length(zets), length(azts), length(lambda), 2, 2, 1, length(type)), dimnames = dimnms)

for(s in 22:length(type)) {
	fls <- list.files(path = "sim/skyrad", pattern = type[s], full.names = T)
	par.fls <- c("sim/skyrad/overcast_simpar.txt", grep(pattern = "_simpar", fls, value = T))
	dif.fls <- c("sim/skyrad/overcast_diffuse.txt", grep(pattern = "_diffuse", fls, value = T))
	dir.fls <- c("sim/skyrad/overcast_direct.txt", grep(pattern = "_direct", fls, value = T))
	sky.fls <- c("sim/skyrad/overcast_sky.txt", grep(pattern = "_sky", fls, value = T))
	for(sza in 4:5) {
		sim_par  <- read_simpar(par.fls[sza])
		dire     <- read.table(dir.fls[sza], header = T)
		dife     <- read.table(dif.fls[sza], header = T)
		skydata  <- read.table(sky.fls[sza], header = T)

		rhow     <- approxfun(coastal[, 1], coastal[, 2], rule = 2)(sim_par$lambda)
		eg       <- dire[, 2] + dife[, 2]
		surfdata <- matrix(NA, ncol = length(sim_par$lambda), nrow = nrow(skydata))
		for(j in 1:length(sim_par$lambda)) {
			surfdata[, j] <- surfacelight(skydata, sim_par$lambda[j], rhow, eg, sim_par)[, 3]
		}
		skydata        <- as.matrix(skydata[, -c(1:2)])
		dim(skydata)   <- c(length(sim_par$zenith), length(sim_par$razimuth) + 180, length(sim_par$lambda))
		dim(surfdata)  <- c(length(sim_par$zenith), length(sim_par$razimuth) + 180, length(sim_par$lambda))
		lightdataS     <- abind(skydata, surfdata[2:91,,], along = 1)
		sun   <- SpatialPoints(cbind(0, 90 - sim_par$sza), proj4string = CRS(crs.wgs))

		for(geom in 1:2) {
			shadow <- get(paste0("shadow_", geom))
			lightraster <- brick(lightdataS, xmn = -180.5, xmx = 180.5, ymn = -90.5, ymx = 90.5)
			lightraster <- mask(lightraster, shadow, inverse = T, updatevalue = 0)
			lightdata   <- as.array(lightraster)
			dim(lightdata) <- c(181*361, length(sim_par$lambda))
			lightdata      <- cbind(expand.grid(0:180, 0:360), lightdata)

			for(zentilt in c(3, 6, 9, 12)) {
				if(geom == 1) tts <- cbind(rep(zentilt, 8), (c(0, 45, 90, 135, 180, 225, 270, 315)+0) %% 360)
				if(geom == 2) tts <- cbind(rep(zentilt, 8), (c(0, 45, 90, 135, 180, 225, 270, 315)+90) %% 360)
				viewB <- viewg[rep(geom, 8), ]
				for(j in 1:nrow(viewB)) {
					dirm <- sphere.to.dir.cos(cbind(rad(viewB[j, 1]), rad(viewB[j, 2])))
					viewB[j,] <- deg(dir.cos.to.sphere(update.cosine(psi = rad(tts[j, 1]), 
					phi = rad(tts[j, 2]), cos.m = dirm)))
				}
				idz  <- which(zets == zentilt)

				for(azitilt in c(0, 45, 90, 135, 180, 225, 270, 315)) {
					skydata <- tilt(lightdata, azitilt, zentilt, sim_par)
					skydata <- as.matrix(skydata[, -c(1:2), drop = F])
					dim(skydata) <- c(length(sim_par$zenith), length(sim_par$razimuth)+180, length(sim_par$lambda))
					ida        <- which(azts == azitilt)

					for(p in 2) {
						if(p == 1) {
							coef <- ptfe_white_all
						} else {
							coef <- ptfe_grey_532
						}
						dir <- paste0(paste0("sim/brdf/tilt/geometry_", geom, "/"), formatC(zentilt, width = 2, flag = 0))
						fl.brdf <- paste0(dir, "/brdf_", ptfes[p], "_", formatC(zentilt, width = 3, flag = 0), "_", formatC(azitilt, width = 3, flag = 0), ".Rda")
						load(fl.brdf)

						for(l in 1:length(lambda)) {
							tiltm    <- sphere.to.dir.cos(cbind(rad(sim_par$sza), rad(180)))
							sun_real <- deg(dir.cos.to.sphere(update.cosine(psi = rad(zentilt), phi = rad((azitilt + 180) %% 360), tiltm)))
							if(sza == 1 || sza == 2) sun_real[2] <- (sun_real[2] + 180) %% 360
							fact <- cos(rad(sun_real[1])) / cos(rad(sim_par$sza))
							if(sun_real[1] > 90) fact <- 0
		
							if(p == 1) {
								scale      <- ptfe_99_rho[l, geom]
								scale_real <- ptfe_99_rho_tilt[l, geom, idz, ida]
							} else {
								scale      <- ptfe_10_rho[l, geom]
								scale_real <- ptfe_10_rho_tilt[l, geom, idz, ida]
							}
							sdata  <- skydata[,, l]
						        tdata  <- sdata * brdf[,, l] * scale_real
							intp   <- bicubic(x = rad(0:90), y = rad(0:360), z = tdata, x0 = gxy[, 1], y0 = gxy[, 2])$z
							intp   <- intp * fx
							dim(intp) <- c(length(sx), length(sy))
							intp   <- intp[-length(sx), ] + intp[-1, ]
							intp   <- intp[, -length(sy)] + intp[, -1]
							tempdf <- sum(intp / 4)

							dbrdf      <- zernike_g(coef, sim_par$sza, 180, viewg[geom, 1], viewg[geom, 2], sim_par$lambda[l])
							dbrdf_real <- zernike_g(coef, sun_real[1], sun_real[2], viewB[ida, 1], viewB[ida, 2], sim_par$lambda[l])
							if(p == 1) {
								dbrdf      <- dbrdf * scale
								dbrdf_real <- dbrdf_real * scale_real
							}
							sunshdw <- gIntersects(shadow, sun)
							if(sunshdw) {
								tempdr <- 0
							} else {
								tempdr <- dire[l, 2] * dbrdf_real * fact
							}

							crosscal <- recalc[l, g]
							meas_pi <- (c(tempdf, tempdr) * pi / crosscal)
							real    <- dife[l, 2] + dire[l, 2]

							crosscal_all_tot_pi[idz, ida, l, sza-3, geom, 1, s] <- sum(meas_pi) / real
							print(c(type[s], round(c(sim_par$sza, sim_par$lambda[l], zentilt, azitilt, crosscal_all_tot_pi[idz, ida, l, sza-3, geom, 1, s]), 3)))
						}
					}
				}
			}
		}
	}
}

save(crosscal_all_tot_pi, file = "res/crosscal_all_tot_pi.Rdata")


#
# Extract average errors: 
#

library(abind)

temp     <- array(NA, dim = c(1, length(azts), length(lambda), 2, 2, 1, length(type)))
for(i in 1:length(azts)) temp[1, i, , , , , ] <- crosscal_shdw_tot_pi
crosscal_lamb_temp <- abind(temp, crosscal_all_tot_pi, along = 1)

(temp1 <- mean(abs(1 - crosscal_lamb_temp[1:5, , 6:36, , 1, 1, 1:24])))
sd(abs(1 - crosscal_lamb_temp[1:5, , 6:36, , 1, 1, 1:24]))
(temp2 <- mean(abs(1 - crosscal_lamb_temp[1:5, , 6:36, , 2, 1, 1:24])))
sd(abs(1 - crosscal_lamb_temp[1:5, , 6:36, , 2, 1, 1:24]))

#
# Extract tilt error for irradiance sensor:
#

temp <- array(1, dim = c(1, length(azts), length(lambda), 10, length(type)))		# Make array for the condition of no tilt. Under no Tilt, no BRDF and no Shadow, Eg estimated == Eg real
res_temp <- abind(temp, res_tilt_tot, along = 1)
irrade <- mean(abs(1 - res_temp[1:3, , 6:36, 4:5, , drop = F]))				# For tilts up to 5 degrees (actually 6 was modeled)

sqrt(temp1^2 + irrade^2)								# Quadrature sum
sqrt(temp2^2 + irrade^2)								# Quadrature sum


