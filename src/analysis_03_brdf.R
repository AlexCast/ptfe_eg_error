# Eg estimation from PTFE plaques
#
# Analysis part 3 - Evaluate BRDF effects for a leveled plaque, for view from nadir and 
#                   40º at 90º azimuth
#
# This script uses a fixed position (90º azimuth to the Sun) and two view angles 0º and 
# 40º and two plaques (white and grey) to evaluate BRDF effects alone, both using the BRDF 
# for the Sun-sensor geometry and the lambertian approximation.
#
# Even in the absence of shadowing or tilt, it is expected the the real BRDF divergence 
# from lambertian will impact the Eg estimation, since the conversion between Lp and Eg 
# can at best be performed with approximations. This effect is studied here using the BRDF 
# of a white and a grey spectralon plaque, published by Germer (2017, 2018). Since the sky 
# radiance data do not account for polarization, only the first Mueller matrix element of 
# the BRDF is used. Note that the sky radiance data is rotated such that the Sun is at 180 
# degrees, and as such the implementation of the functions also rotate for the convention 
# of Sun at 0 degrees.
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

#
# Perform calculations:
#

viewg <- rbind(c(0, 270), c(40, 270))
ptfes <- c("ptfe99", "ptfe10")

dimnms <- list(lambda, c('OC', paste0('ts', formatC(seq(0, 80, 10), width = 2, flag = 0))), 
		c('tv00', 'tv40'), ptfes, type)

res_brdf_tot    <- array(NA, dim = c(length(lambda), 10, nrow(viewg), 2, length(type)), dimnames = dimnms)
res_brdf_tot_pi <- array(NA, dim = c(length(lambda), 10, nrow(viewg), 2, length(type)), dimnames = dimnms)

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
				sdata <- skydata[,,l]
				for(geom in 1:nrow(viewg)) {
					if(p == 1) {
						scale <- ptfe_99_rho[l, geom]
						brdf  <- brdf_ptfe99[[geom]][,,l]
					} else {
						scale <- ptfe_10_rho[l, geom]
						brdf  <- brdf_ptfe10[[geom]][,,l]
					}

					tdata  <- sdata * brdf * scale
					intp   <- bicubic(x = rad(0:90), y = rad(0:360), z = tdata, x0 = gxy[, 1], y0 = gxy[, 2])$z
					intp   <- intp * fx
					dim(intp) <- c(length(sx), length(sy))
					intp   <- intp[-length(sx), ] + intp[-1, ]
					intp   <- intp[, -length(sy)] + intp[, -1]
					tempdf <- sum(intp / 4)

					dbrdf   <- zernike_g(coef, sim_par$sza, 180, viewg[geom, 1], viewg[geom, 2], sim_par$lambda[l])
					if(p == 1) {
						dbrdf <- dbrdf * scale			# The coefficients for the white plaque give a normalized BRDF
					}
					tempdr  <- dire[l, 2] * dbrdf
					meas    <- ((tempdf + tempdr) / dbrdf)
					meas_pi <- ((tempdf + tempdr) * pi / scale)
					real    <- dife[l, 2] + dire[l, 2]
					res_brdf_tot[l, sza, geom, p, s]    <- meas/real
					res_brdf_tot_pi[l, sza, geom, p, s] <- meas_pi/real
					print(c(type[s], sim_par$sza, sim_par$lambda[l], geom, round(res_brdf_tot[l, sza, geom, p, s], 4), round(res_brdf_tot_pi[l, sza, geom, p, s], 4)))
				}
			}
		}
	}
}

save(res_brdf_tot, file = "res/res_brdf_tot.Rdata")
save(res_brdf_tot_pi, file = "res/res_brdf_tot_pi.Rdata")

#
# Extract error statistics. 
# Get the mean absolute percentage error and the standard deviation of the absolute errors
# for the visible range under recommended skies and  overcast conditions. Dimensions are: 
# sky (SZA from 20 to 60 = 1, and OC = 2), geometry, plaque, inversion (BRDF or pi).
#

dimnms <- list(c('RSky', 'OC'), c('tv00', 'tv40'), ptfes, c('Eq5A', 'Eq5B'))

error_brdf_avem <- array(NA, dim = c(2, 2, 2, 2), dimnames = dimnms)
error_brdf_sdm  <- array(NA, dim = c(2, 2, 2, 2), dimnames = dimnms)

for(geom in 1:2) {
	for(p in 1:2) {
		temp <- abs(1 - res_brdf_tot[6:36, 4:8, geom, p, ])
		error_brdf_avem[1, geom, p, 1] <- mean(temp)
		error_brdf_sdm[1, geom, p, 1]  <- sd(temp)

		temp <- abs(1 - res_brdf_tot[6:36, 1, geom, p, 1])
		error_brdf_avem[2, geom, p, 1] <- mean(temp)
		error_brdf_sdm[2, geom, p, 1]  <- sd(temp)

		temp <- abs(1 - res_brdf_tot_pi[6:36, 4:8, geom, p, ])
		error_brdf_avem[1, geom, p, 2] <- mean(temp)
		error_brdf_sdm[1, geom, p, 2]  <- sd(temp)

		temp <- abs(1 - res_brdf_tot_pi[6:36, 1, geom, p, 1])
		error_brdf_avem[2, geom, p, 2] <- mean(temp)
		error_brdf_sdm[2, geom, p, 2]  <- sd(temp)
	}
}

round(error_brdf_avem * 100, 2)
round(error_brdf_sdm * 100, 2)

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

viewg <- rbind(c(0, 90), c(40, 90))

for(p in 1:2) {
	for(geom in 1:2) {
		tikz(paste0("plot/plaque_brdf_effects_p", p, "_", viewg[geom, 1], "_", viewg[geom, 2], ".tex"), width = 7, height = 7, standAlone = TRUE, pointsize = 18,
			packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
			"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
			"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
			par(mar = c(5, 5, 3, 2))
			xlab <- '$\\lambda$ (nm)'
			ylab <- '$E_{\\text{g}}^{\\text{ est.}}/E_{\\text{g}}^{\\text{{ real}}}$'
		 	plot(lambda, apply(res_brdf_tot[, 1, geom, p, ], 1, mean), type = "l", ylim = c(0.85, 1.15), col = cols[1], xlab = '', ylab = '', lwd = 3, xaxt = "n", yaxt = "n")
			axis(2, at = seq(0.85, 1.15, 0.05), labels = T, cex.axis = 1.9)
			axis(1, at = seq(400, 1000, 100), labels = T, cex.axis = 1.9)
			mtext(xlab, side = 1, line = 3, cex = 2)
			mtext(ylab, side = 2, line = 3, cex = 2)
		 	lines(lambda, apply(res_brdf_tot_pi[, 1, geom, p, ], 1, mean), lty = 2, col = cols[1], lwd = 3)
		 	lines(lambda, apply(res_brdf_tot[, 2, geom, p, ], 1, mean), lty = 1, col = cols[2], lwd = 3)
		 	lines(lambda, apply(res_brdf_tot_pi[, 2, geom, p, ], 1, mean), lty = 2, col = cols[2], lwd = 3)
		 	lines(lambda, apply(res_brdf_tot[, 4, geom, p, ], 1, mean), lty = 1, col = cols[3], lwd = 3)
		  	lines(lambda, apply(res_brdf_tot_pi[, 4, geom, p, ], 1, mean), lty = 2, col = cols[3], lwd = 3)
			lines(lambda, apply(res_brdf_tot[, 6, geom, p, ], 1, mean), lty = 1, col = cols[4], lwd = 3)
			lines(lambda, apply(res_brdf_tot_pi[, 6, geom, p, ], 1, mean), lty = 2, col = cols[4], lwd = 3)
			lines(lambda, apply(res_brdf_tot[, 8, geom, p, ], 1, mean), lty = 1, col = cols[5], lwd = 3)
			lines(lambda, apply(res_brdf_tot_pi[, 8, geom, p, ], 1, mean), lty = 2, col = cols[5], lwd = 3)
			lines(lambda, apply(res_brdf_tot[, 10, geom, p, ], 1, mean), lty = 1, col = cols[6], lwd = 3)
			lines(lambda, apply(res_brdf_tot_pi[, 10, geom, p, ], 1, mean), lty = 2, col = cols[6], lwd = 3)
			legend("bottomleft", c('OC', '$\\theta_{\\text{s}}=0^{\\circ}$', '$\\theta_{\\text{s}}=20^{\\circ}$', '$\\theta_{\\text{s}}=40^{\\circ}$', '$\\theta_{\\text{s}}=60^{\\circ}$', 
				'$\\theta_{\\text{s}}=80^{\\circ}$'), lty = 1, col = cols, bty = "n", cex = 1.2, lwd = 3, y.intersp = 0.5)
			legend("topright", c('$L_{\\text{p}}$ to $E_{\\text{g}}$ convertion based on $f_{\\text{r}}(\\theta_{\\text{s}}, \\phi_{\\text{s}}, \\theta_{\\text{v}}, \\phi_{\\text{v}})^{-1}$', 
		         "Lambertian assumption"), lty = c(1, 2), col = c("black", "black"), lwd = 3, bty = "n", cex = 1.2, y.intersp = 0.5)
			legend("bottomright", c(paste(plaque[p], "PTFE"), eval(substitute(expression('$\\theta_{\\text{v}}$'==a*'$^{\\circ}$,'~'$\\phi_{\\text{v}}$'==b*'$^{\\circ}$'), list(a = viewg[geom, 1], b = viewg[geom, 2])))), bty = "n", cex = 1.5, y.intersp = 0.5)
			abline(h = 1, lty = 2, col = "grey")
	  	dev.off()
		tools::texi2pdf(paste0("plot/plaque_brdf_effects_p", p, "_", viewg[geom, 1], "_", viewg[geom, 2], ".tex"))
		file.rename(paste0("plaque_brdf_effects_p", p, "_", viewg[geom, 1], "_", viewg[geom, 2], ".pdf"), paste0("plot/plaque_brdf_effects_p", p, "_", viewg[geom, 1], "_", viewg[geom, 2], ".pdf"))
		file.remove(paste0("plaque_brdf_effects_p", p, "_", viewg[geom, 1], "_", viewg[geom, 2], ".aux"), paste0("plaque_brdf_effects_p", p, "_", viewg[geom, 1], "_", viewg[geom, 2], ".log"))
	}
}


diff_f <- matrix(NA, ncol = length(dif.fls), nrow = length(lambda))
for(i in 1:length(dif.fls))
 {
  dire     <- read.table(dir.fls[i], header = T)
  dife     <- read.table(dif.fls[i], header = T)
  diff_f[, i] <- dife[, 2] / (dire[, 2] + dife[, 2])
 }

dvec <- c(1, 2, 4, 6, 8, 10)

for(geom in 2)
 {
  png(paste0("plot/error_brdf_diffuse_fraction_", viewg[geom, 1], "_", viewg[geom, 2], ".png"), height = 1200, width = 1200, res = 230, bg = "transparent")
  xlab <- expression("Diffuse fraction "*E[d]/E[g])
  ylab <- expression(epsilon~("%"))
  par(mar = c(5, 5, 3, 2))
  plot(NA, xlim = range(diff_f[1:58, dvec]), ylim = range(abs(1 - res_brdf[, dvec, geom])) * 100, xlab = xlab, ylab = ylab)
  for(i in 1:6) points(diff_f[1:58, dvec[i]], abs(1 - res_brdf[, dvec[i], geom]) * 100, col = cols[i])
  legend("topleft", c("OC", expression(theta[s]==0*"º"), expression(theta[s]==20*"º"), expression(theta[s]==40*"º"), expression(theta[s]==60*"º"), 
       expression(theta[s]==80*"º")), pch = 1, col = cols, bty = "n", cex = 0.9)
  legend("bottomright", eval(substitute(expression(theta[v]==a*"º,"~phi[v]==b*"º"), list(a = viewg[geom, 1], b = viewg[geom, 2]))), bty = "n", cex = 0.9)
  dev.off()
 }



