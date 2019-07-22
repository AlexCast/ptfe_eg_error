# Eg estimation from PTFE plaques
#
# Analysis part 4 - Evaluate shadow effects with a leveled plaque, for view from nadir and 
#                   40º at 90º azimuth
#
# This script uses a fixed position (90º azimuth to the Sun) and two view angles 0º and 
# 40º to evaluate shadow effects alone and shadow effects interaction with BRDF effects, 
# both using the BRDF at the Sun-sensor geometry and the lambertian approximation. 
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
library(akima)
library(rgeos)

#
# Perform calculations shadow effect alone:
# 

viewg <- rbind(c(0, 270), c(40, 270))
shdwl <- list(shadow_1, shadow_2)

dimnms <- list(lambda, c('OC', paste0('ts', formatC(seq(0, 80, 10), width = 2, flag = 0))), 
		c('tv00', 'tv40'), type)

res_shdw_dif_nbrdf <- array(NA, dim = c(length(lambda), 10, nrow(viewg), length(type)), dimnames = dimnms)
res_shdw_dir_nbrdf <- array(NA, dim = c(length(lambda), 10, nrow(viewg), length(type)), dimnames = dimnms)
res_shdw_tot_nbrdf <- array(NA, dim = c(length(lambda), 10, nrow(viewg), length(type)), dimnames = dimnms)

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
		sun   <- SpatialPoints(cbind(0, 90 - sim_par$sza), proj4string = CRS(crs.wgs))
		for(l in 1:length(lambda)) {
			sdata <- skydata[,,l]
			for(geom in 1:nrow(viewg)) {

				shadow  <- shdwl[[geom]]
				temp    <- mask(raster(sdata, xmn=-180.5, xmx=180.5, ymn=-0.5, ymx=90.5, crs = crs.wgs), shadow, inverse = T)
				tsdata  <- matrix(values(temp), ncol = 361, byrow = T)
				tsdata[is.na(tsdata)] <- 0

				tdata  <- tsdata
				intp   <- bicubic(x = rad(0:90), y = rad(0:360), z = tdata, x0 = gxy[, 1], y0 = gxy[, 2])$z
				intp   <- intp * fx
				dim(intp) <- c(length(sx), length(sy))
				intp   <- intp[-length(sx), ] + intp[-1, ]
				intp   <- intp[, -length(sy)] + intp[, -1]
				tempdf <- sum(intp / 4)

				sunshdw <- gIntersects(shadow, sun)
				if(sunshdw) {
					tempdr <- 0
				} else {
					tempdr <- dire[l, 2]
				}

				meas    <- c(tempdf, tempdr)
				real    <- dife[l, 2] + dire[l, 2]

				res_shdw_dif_nbrdf[l, sza, geom, s] <- meas[1] / dife[l, 2]
				res_shdw_dir_nbrdf[l, sza, geom, s] <- meas[2] / dire[l, 2]
				res_shdw_tot_nbrdf[l, sza, geom, s] <- sum(meas) / real
				print(c(type[s], sim_par$sza, sim_par$lambda[l], geom, round(res_shdw_tot_nbrdf[l, sza, geom, s], 4)))
			}
		}
	}
	save(res_shdw_dif_nbrdf, file = "res/res_shdw_dif_nbrdf.Rdata")
	save(res_shdw_dir_nbrdf, file = "res/res_shdw_dir_nbrdf.Rdata")
	save(res_shdw_tot_nbrdf, file = "res/res_shdw_tot_nbrdf.Rdata")
}


#
# Perform calculations for BRDF + shadow interaction:
# 

ptfes <- c("ptfe99", "ptfe10")

dimnms <- list(lambda, c('OC', paste0('ts', formatC(seq(0, 80, 10), width = 2, flag = 0))), 
		c('tv00', 'tv40'), ptfes, type)

res_shdw_dif    <- array(NA, dim = c(length(lambda), 10, nrow(viewg), 2, length(type)), dimnames = dimnms)
res_shdw_dir    <- array(NA, dim = c(length(lambda), 10, nrow(viewg), 2, length(type)), dimnames = dimnms)
res_shdw_tot    <- array(NA, dim = c(length(lambda), 10, nrow(viewg), 2, length(type)), dimnames = dimnms)
res_shdw_tot_pi <- array(NA, dim = c(length(lambda), 10, nrow(viewg), 2, length(type)), dimnames = dimnms)

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
		sun   <- SpatialPoints(cbind(0, 90 - sim_par$sza), proj4string = CRS(crs.wgs))
		for(p in 1:2) {
			if(p == 1) {
				coef <- ptfe_white_all
			} else {
				coef <- ptfe_grey_532
			}
			for(l in 1:length(lambda)) {
				sdata <- skydata[,,l]
				for(geom in 1:nrow(viewg)) {
					if(p == 1) {
						scale <- ptfe_99_rho[l, geom]
						brdf  <- brdf_ptfe99[[geom]][,,l]
					} else {
						scale <- ptfe_10_rho[l, geom]
						brdf  <- brdf_ptfe10[[geom]][,,l]
					}
	
					shadow  <- shdwl[[geom]]
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

					dbrdf   <- zernike_g(coef, sim_par$sza, 180, viewg[geom, 1], viewg[geom, 2], sim_par$lambda[l])
					if(p == 1) {
						dbrdf <- dbrdf * scale
					}

					sunshdw <- gIntersects(shadow, sun)		
					if(sunshdw) {
						tempdr <- 0
					} else {
						tempdr <- dire[l, 2] * dbrdf
					}

					meas    <- (c(tempdf, tempdr) / dbrdf)
					meas_pi <- (c(tempdf, tempdr) * pi / scale)
					real    <- dife[l, 2] + dire[l, 2]

					res_shdw_dif[l, sza, geom, p, s] <- meas[1] / dife[l, 2]
					res_shdw_dir[l, sza, geom, p, s] <- meas[2] / dire[l, 2]
					res_shdw_tot[l, sza, geom, p, s] <- sum(meas) / real
					res_shdw_tot_pi[l, sza, geom, p, s] <- sum(meas_pi) / real
					print(c(type[s], sim_par$sza, sim_par$lambda[l], geom, round(res_shdw_tot[l, sza, geom, p, s], 4), round(res_shdw_tot_pi[l, sza, geom, p, s], 4)))
				}
			}
		}
	}
	save(res_shdw_dif, file = "res/res_shdw_dif.Rdata")
	save(res_shdw_dir, file = "res/res_shdw_dir.Rdata")
	save(res_shdw_tot, file = "res/res_shdw_tot.Rdata")
	save(res_shdw_tot_pi, file = "res/res_shdw_tot_pi.Rdata")
}

#
# Extract error statistics. 
# Get the mean absolute percentage error and the standard deviation of the absolute errors
# for the visible range under recommended skies and  overcast conditions. Dimensions are: 
# sky (SZA from 20 to 60 = 1, and OC = 2) and geometry for shadow effects alone and
# sky (SZA from 20 to 60 = 1, and OC = 2), geometry, plaque and inversion (BRDF or pi) for
# BRDF and shadow interaction.
#

dimnms <- list(c('RSky', 'OC'), c('tv00', 'tv40'))

error_shdw_avem <- array(NA, dim = c(2, 2), dimnames = dimnms)
error_shdw_sdm  <- array(NA, dim = c(2, 2), dimnames = dimnms)

for(geom in 1:2) {
	for(p in 1:2) {
		temp <- abs(1 - res_shdw_tot_nbrdf[6:36, 4:8, geom, ])
		error_shdw_avem[1, geom] <- mean(temp)
		error_shdw_sdm[1, geom]  <- sd(temp)

		temp <- abs(1 - res_shdw_tot_nbrdf[6:36, 1, geom, 1])
		error_shdw_avem[2, geom] <- mean(temp)
		error_shdw_sdm[2, geom]  <- sd(temp)
	}
}

round(error_shdw_avem * 100, 2)
round(error_shdw_sdm * 100, 2)

dimnms <- list(c('RSky', 'OC'), c('tv00', 'tv40'), ptfes, c('Eq5A', 'Eq5B'))

error_shdw_brdf_avem <- array(NA, dim = c(2, 2, 2, 2), dimnames = dimnms)
error_shdw_brdf_sdm  <- array(NA, dim = c(2, 2, 2, 2), dimnames = dimnms)

for(geom in 1:2) {
	for(p in 1:2) {
		temp <- abs(1 - res_shdw_tot[6:36, 4:8, geom, p, ])
		error_shdw_brdf_avem[1, geom, p, 1] <- mean(temp)
		error_shdw_brdf_sdm[1, geom, p, 1]  <- sd(temp)

		temp <- abs(1 - res_shdw_tot[6:36, 1, geom, p, 1])
		error_shdw_brdf_avem[2, geom, p, 1] <- mean(temp)
		error_shdw_brdf_sdm[2, geom, p, 1]  <- sd(temp)

		temp <- abs(1 - res_shdw_tot_pi[6:36, 4:8, geom, p, ])
		error_shdw_brdf_avem[1, geom, p, 2] <- mean(temp)
		error_shdw_brdf_sdm[1, geom, p, 2]  <- sd(temp)

		temp <- abs(1 - res_shdw_tot_pi[6:36, 1, geom, p, 1])
		error_shdw_brdf_avem[2, geom, p, 2] <- mean(temp)
		error_shdw_brdf_sdm[2, geom, p, 2]  <- sd(temp)
	}
}

round(error_shdw_brdf_avem * 100, 2)
round(error_shdw_brdf_sdm * 100, 2)

#
# Graphical result:
#

library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
	"\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
	"\\documentclass[18pt]{article}")

cols <- rev(rainbow(6, start = 0, end = 0.8))
cols[5] <- "#FFB84D"
plaque <- c("White", "Grey")
leg <- c("B", "A")

viewg <- rbind(c(0, 90), c(40, 90))

for(g in 1:2) {
	tikz(paste0("plot/plaque_shdw_effects_p0_", viewg[g, 1], "_", viewg[g, 2], ".tex"), width = 7, height = 7, standAlone = TRUE, pointsize = 18,
		packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
		"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
		"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
		par(mar = c(5, 5, 3, 2))
		xlab <- '$\\lambda$ (nm)'
		ylab <- '$E_{\\text{g}}^{\\text{ est.}}/E_{\\text{g}}^{\\text{{ real}}}$'
		plot(lambda, apply(res_shdw_tot_nbrdf[, 1, g, ], 1, mean), type = "l", ylim = c(0.935, 1), col = cols[1], xlab = '', ylab = '', lwd = 3, xaxt = "n", yaxt = "n")
		axis(2, at = seq(0.94, 1, 0.02), labels = T, cex.axis = 1.9)
		axis(1, at = seq(400, 1000, 100), labels = T, cex.axis = 1.9)
		mtext(xlab, side = 1, line = 3, cex = 2)
		mtext(ylab, side = 2, line = 3, cex = 2)
		lines(lambda, apply(res_shdw_tot_nbrdf[, 2, g, ], 1, mean), lty = 1, col = cols[2], lwd = 3)
		lines(lambda, apply(res_shdw_tot_nbrdf[, 4, g, ], 1, mean), lty = 1, col = cols[3], lwd = 3)
		lines(lambda, apply(res_shdw_tot_nbrdf[, 6, g, ], 1, mean), lty = 1, col = cols[4], lwd = 3)
		lines(lambda, apply(res_shdw_tot_nbrdf[, 8, g, ], 1, mean), lty = 1, col = cols[5], lwd = 3)
		lines(lambda, apply(res_shdw_tot_nbrdf[, 10, g, ], 1, mean), lty = 1, col = cols[6], lwd = 3)
#		legend("bottomleft", c('OC', '$\\theta_{\\text{s}}=0^{\\circ}$', '$\\theta_{\\text{s}}=20^{\\circ}$', '$\\theta_{\\text{s}}=40^{\\circ}$', '$\\theta_{\\text{s}}=60^{\\circ}$', 
#			'$\\theta_{\\text{s}}=80^{\\circ}$'), lty = 1, col = cols, bty = "n", cex = 1.2, lwd = 3, y.intersp = 0.5)
		legend("bottomright", eval(substitute(expression('$\\theta_{\\text{v}}$'==a*'$^{\\circ}$,'~'$\\phi_{\\text{v}}$'==b*'$^{\\circ}$'), list(a = viewg[g, 1], b = viewg[g, 2]))), bty = "n", cex = 1.5, y.intersp = 0.5)
		abline(h = 1, lty = 2, col = "grey")
#		legend("right", eval(substitute(expression(bold(letter)), list(letter = leg[g]))), bty = "n")
	dev.off()
	tools::texi2pdf(paste0("plot/plaque_shdw_effects_p0_", viewg[g, 1], "_", viewg[g, 2], ".tex"))
	file.rename(paste0("plaque_shdw_effects_p0_", viewg[g, 1], "_", viewg[g, 2], ".pdf"), paste0("plot/plaque_shdw_effects_p0_", viewg[g, 1], "_", viewg[g, 2], ".pdf"))
	file.remove(paste0("plaque_shdw_effects_p0_", viewg[g, 1], "_", viewg[g, 2], ".aux"), paste0("plaque_shdw_effects_p0_", viewg[g, 1], "_", viewg[g, 2], ".log"))
}

diff_f <- matrix(NA, ncol = length(dif.fls), nrow = length(lambda))
for(i in 1:length(dif.fls))
 {
  dire     <- read.table(dir.fls[i], header = T)
  dife     <- read.table(dif.fls[i], header = T)
  diff_f[, i] <- dife[, 2] / (dire[, 2] + dife[, 2])
 }

dvec <- c(1, 2, 4, 6, 8, 10)

for(g in 2)
 {
  png(paste0("plot/error_shdw_diffuse_fraction_", viewg[g, 1], "_", viewg[g, 2], ".png"), height = 1200, width = 1200, res = 230, bg = "transparent")
  xlab <- expression("Diffuse fraction "*E[d]/E[g])
  ylab <- expression(epsilon~("%"))
  par(mar = c(5, 5, 3, 2))
  plot(NA, xlim = range(diff_f[1:58, dvec]), ylim = range(abs(1 - res_shdw[, dvec, g] / res_brdf[, dvec, g])) * 100, xlab = xlab, ylab = ylab)
  for(i in 1:6) points(diff_f[1:58, dvec[i]], abs(1 - res_shdw[, dvec[i], g] / res_brdf[, dvec[i], g])  * 100, col = cols[i])
  legend("topleft", c("OC", expression(theta[s]==0*"º"), expression(theta[s]==20*"º"), expression(theta[s]==40*"º"), expression(theta[s]==60*"º"), 
       expression(theta[s]==80*"º")), pch = 1, col = cols, bty = "n", cex = 0.9)
  legend("bottomright", eval(substitute(expression(theta[v]==a*"º,"~phi[v]==b*"º"), list(a = viewg[g, 1], b = viewg[g, 2]))), bty = "n", cex = 0.9)
  dev.off()
 }


