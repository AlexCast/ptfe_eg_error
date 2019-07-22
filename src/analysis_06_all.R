# Eg estimation from PTFE plaques
#
# Analysis part 6 - Evaluate combined effects of BRDF, shadow and tilt
#
# Tilt effects are evaluated in combination with other effects. Since BRDF effects are 
# included, the BRDF and rho at the real view angle after tilt is used for calculations of 
# Lp, while the nominal BRDF and rho are used for the estimation of Eg, as would occur 
# under unintentional tilt. Because tilt also affect the shadow position since plaque tilt 
# is independent of structure/sensor, and because rotation of the shadow vector would 
# possibly produce invalid geometries, shadow mask is applied first before rotation of the 
# radiance field. If sensor and plaque were fixed and so experiencing the same tilt, 
# shadow mask would be applied after rotation - but this condition is not considered here.
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

library(akima)
library(abind)
library(sp)
library(raster)
library(rgeos)

#
# Perform calculations:
#

zets <- c(3, 6, 9, 12)
azts <- c(0, 45, 90, 135, 180, 225, 270, 315)

viewg <- rbind(c(0, 270), c(40, 270))
ptfes <- c("ptfe99", "ptfe10")

dimnms <- list(c('znt03', 'znt06', 'znt09', 'znt12'), c('azt000', 'azt045', 'azt090', 
		'azt135', 'azt180', 'azt225', 'azt270', 'azt315'), lambda, c('OC', 
		paste0('ts', formatC(seq(0, 80, 10), width = 2, flag = 0))), c('tv00', 
		'tv40'), ptfes, type)

res_all_dif <- array(NA, dim = c(length(zets), length(azts), length(lambda), 10, 2, 2, length(type)), dimnames = dimnms)
res_all_dir <- array(NA, dim = c(length(zets), length(azts), length(lambda), 10, 2, 2, length(type)), dimnames = dimnms)
res_all_tot <- array(NA, dim = c(length(zets), length(azts), length(lambda), 10, 2, 2, length(type)), dimnames = dimnms)
res_all_tot_pi <- array(NA, dim = c(length(zets), length(azts), length(lambda), 10, 2, 2, length(type)), dimnames = dimnms)

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

					for(p in 1:2) {
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

							meas    <- (c(tempdf, tempdr) / dbrdf)
							meas_pi <- (c(tempdf, tempdr) * pi / scale)
							real    <- dife[l, 2] + dire[l, 2]

							res_all_dif[idz, ida, l, sza, geom, p, s] <- meas[1] / dife[l, 2]
							res_all_dir[idz, ida, l, sza, geom, p, s] <- meas[2] / dire[l, 2]
							res_all_tot[idz, ida, l, sza, geom, p, s] <- sum(meas) / real
							res_all_tot_pi[idz, ida, l, sza, geom, p, s] <- sum(meas_pi) / real
							print(c(type[s], round(c(sim_par$sza, sim_par$lambda[l], zentilt, azitilt, res_all_tot[idz, ida, l, sza, geom, p, s], res_all_tot_pi[idz, ida, l, sza, geom, p, s]), 3)))
						}
					}
				}
			}
		}
	}
	save(res_all_dif, file = "res/res_all_dif.Rdata")
	save(res_all_dir, file = "res/res_all_dir.Rdata")
	save(res_all_tot, file = "res/res_all_tot.Rdata")
	save(res_all_tot_pi, file = "res/res_all_tot_pi.Rdata")
}


#
# Extract error statistics. 
# Get the mean absolute percentage error and the standard deviation of the absolute errors
# for the visible range under recommended skies and  overcast conditions. Note average 
# errors are computed from zero tilt to each step of tilt. Dimensions are: 
# maximum tilt, sky (SZA from 20 to 60 = 1, and OC = 2), geometry, plaque and inversion 
# (BRDF or pi)
#

temp     <- array(NA, dim = c(1, length(azts), length(lambda), 10, 2, 2, length(type)))	# Make array for the condition of no tilt. Values are those calculated for BRDF + shadow.
for(i in 1:length(azts)) temp[1, i, , , , , ] <- res_shdw_tot
res_brdf_temp <- abind(temp, res_all_tot, along = 1)
temp     <- array(NA, dim = c(1, length(azts), length(lambda), 10, 2, 2, length(type)))
for(i in 1:length(azts)) temp[1, i, , , , , ] <- res_shdw_tot_pi
res_lamb_temp <- abind(temp, res_all_tot_pi, along = 1)

dimnms <- list(c('znt03', 'znt06', 'znt09', 'znt12'), c('RSky', 'OC'), c('tv00', 'tv40'), 
		ptfes, c('Eq5A', 'Eq5B'))

error_all_avem <- array(NA, dim = c(4, 2, 2, 2, 2), dimnames = dimnms)
error_all_sdm  <- array(NA, dim = c(4, 2, 2, 2, 2), dimnames = dimnms)

for(g in 1:2) {
	for(p in 1:2) {
		for(i in 2:5) {
			temp <- abs(1 - res_brdf_temp[1:i, , 6:36, 4:8, g, p, ])
			error_all_avem[i-1, 1, g, p, 1] <- mean(temp)
			error_all_sdm[i-1, 1, g, p, 1]  <- sd(temp)

			temp <- abs(1 - res_brdf_temp[1:i, , 6:36, 1, g, p, ])
			error_all_avem[i-1, 2, g, p, 1] <- mean(temp)
			error_all_sdm[i-1, 2, g, p, 1]  <- sd(temp)

			temp <- abs(1 - res_lamb_temp[1:i, , 6:36, 4:8, g, p,])
			error_all_avem[i-1, 1, g, p, 2] <- mean(temp)
			error_all_sdm[i-1, 1, g, p, 2]  <- sd(temp)

			temp <- abs(1 - res_lamb_temp[1:i, , 6:36, 1, g, p,])
			error_all_avem[i-1, 2, g, p, 2] <- mean(temp)
			error_all_sdm[i-1, 2, g, p, 2]  <- sd(temp)

		}
	}
}

round(error_all_avem * 100, 1)
round(error_all_sdm * 100, 1)


# Graphical result:

mean_all_error    <- array(NA, c(58, 10, 2, 2))
mean_all_error_pi <- array(NA, c(58, 10, 2, 2))

for(p in 1:2) {
	for(geom in 1:2) {
		for(sza in 1:10) {
			for(lbd in 1:58) {
				xm <- res_all_tot[, , lbd, sza, geom, p, ]
				mean_all_error[lbd, sza, geom, p] <- mean(abs(1 - xm), na.rm = T)
				xm <- res_all_tot_pi[, , lbd, sza, geom, p, ]
				mean_all_error_pi[lbd, sza, geom, p] <- mean(abs(1 - xm), na.rm = T)
			}
		}
	}
}

plaque <- c("White", "Grey")

for(p in 1:2) {
	for(g in 1:2) {
		png(paste0("plot/plaque_error_brdfEshdwEtilt_effects_p", p, "_", viewg[g, 1], "_", viewg[g, 2], ".png"), height = 1200, width = 1200, res = 230, bg = "white")
		cols <- rev(rainbow(6, start = 0, end = 0.8))
		par(mar = c(5, 5, 3, 2))
		xlab = expression(lambda~(nm))
		ylab = expression("Average"~epsilon~("%"))
		plot(lambda[1:58],  mean_all_error[, 1, g, p]*100, type = "l", ylim = range(mean_all_error, na.rm = T)*100, col = cols[1], xlab = xlab, ylab = ylab, lwd = 1.5)
		lines(lambda[1:58], mean_all_error_pi[, 1, g, p]*100, col = cols[1], lwd = 1.5, lty = 2)
		lines(lambda[1:58], mean_all_error[, 2, g, p]*100, col = cols[2], lwd = 1.5)
		lines(lambda[1:58], mean_all_error_pi[, 2, g, p]*100, col = cols[2], lwd = 1.5, lty = 2)
		lines(lambda[1:58], mean_all_error[, 4, g, p]*100, col = cols[3], lwd = 1.5)
		lines(lambda[1:58], mean_all_error_pi[, 4, g, p]*100, col = cols[3], lwd = 1.5, lty = 2)
		lines(lambda[1:58], mean_all_error[, 6, g, p]*100, col = cols[4], lwd = 1.5)
		lines(lambda[1:58], mean_all_error_pi[, 6, g, p]*100, col = cols[4], lwd = 1.5, lty = 2)
		lines(lambda[1:58], mean_all_error[, 8, g, p]*100, col = cols[5], lwd = 1.5)
		lines(lambda[1:58], mean_all_error_pi[, 8, g, p]*100, col = cols[5], lwd = 1.5, lty = 2)
		lines(lambda[1:58], mean_all_error[, 10, g, p]*100, col = cols[6], lwd = 1.5)
		lines(lambda[1:58], mean_all_error_pi[, 10, g, p]*100, col = cols[6], lwd = 1.5, lty = 2)
		legend("right", c("OC", expression(theta[s]==0), expression(theta[s]==20), expression(theta[s]==40), expression(theta[s]==60), 
	  	     expression(theta[s]==80)), lty = 1, col = cols, bty = "n", cex = 0.8, lwd = 1.5)
#		legend("bottomright", c(expression(L[v]~to~E[g]*" convertion based on BRDF"*(theta[s]*", "*phi[s]*", "*theta[v]*", "*phi[v])^-1), 
#	         "Lambertian assumption"), lty = c(1, 3), col = c("black", "black"), lwd = 1.5, bty = "n", cex = 0.8)
		legend("topleft", c(paste(plaque[p], "PTFE"), eval(substitute(expression(theta[v]==a*"º,"~phi[v]==b*"º"), list(a = viewg[g, 1], b = viewg[g, 2]-180)))), bty = "n", cex = 0.9)
		dev.off()
	}
}


