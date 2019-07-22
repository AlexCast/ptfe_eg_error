# Eg estimation from PTFE plaques
#
# Analysis part 5: Evaluate Tilt effect in absence of other error sources.
#
# Tilt effects are evaluated in isolation from other error sources considered in this
# study: BRDF effects and shadowing effects. That condition can be replicated by 
# calculating plaque emergent radiance at the view angles for a perfectly lambertian 
# surface, in the absence of disturbances to the light field. For efficiency, the code 
# actually computes the irradiance estimation by direct integration of the incidence 
# radiance distribution since for a lambertian plaque there is a perfect relation between 
# plaque indirect estimates and direct observation, regardless of plaque view angles.
#
# Tilt effects are evaluated by rotating the radiance distribution (including upwelling) 
# so to simulate the incidence light field over the plaque surface under a tilted 
# condition. Tilts are defined as angular departures from the plaque normal to the local 
# zenith, in the directions specified as azimuth angles from Sun direction in a counter-
# clockwise rotation. The tilts were evaluated as the combination of zenith tilts at 3, 6, 
# 9 and 12 degrees and azimuth tilts at 0, 45, 90, 135, 180, 225, 270 and 315. But since 
# in the absence of BRDF and shadow effects tilt errors depend only on relative azimuth
# (symmetry is preserved), calculations are only made for 0, 45, 90, 135, 180 and mirrored
# for 45 -> 315 (-45), 90 -> 270 (-90) and 135 -> 225 (-135).
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
library(raster)
library(sp)

#
# Perform calculations:
#

dimnms <- list(c('znt03', 'znt06', 'znt09', 'znt12'), c('azt000', 'azt045', 'azt090', 
		'azt135', 'azt180', 'azt225', 'azt270', 'azt315'), lambda, c('OC', 
		paste0('ts', formatC(seq(0, 80, 10), width = 2, flag = 0))), type)

zets <- c(3, 6, 9, 12)
azts <- c(0, 45, 90, 135, 180, 225, 270, 315)

res_tilt_dif <- array(NA, dim = c(length(zets), length(azts), length(lambda), 10, length(type)), dimnames = dimnms)
res_tilt_dir <- array(NA, dim = c(length(zets), length(azts), length(lambda), 10, length(type)), dimnames = dimnms)
res_tilt_tot <- array(NA, dim = c(length(zets), length(azts), length(lambda), 10, length(type)), dimnames = dimnms)

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
		lightdata      <- abind(skydata, surfdata[2:91,,], along = 1)
		dim(lightdata) <- c(181*361, length(sim_par$lambda))							
		lightdata      <- cbind(expand.grid(0:180, 0:360), lightdata)						

		for(zentilt in c(3, 6, 9, 12)) {
			for(azitilt in c(0, 45, 90, 135, 180)) {							
				skydata <- tilt(lightdata, azitilt, zentilt, sim_par)					
				skydata <- as.matrix(skydata[, -c(1:2), drop = F])					
				dim(skydata) <- c(length(sim_par$zenith), length(sim_par$razimuth)+180, length(sim_par$lambda))

				for(l in 1:length(lambda)) {
					sdata <- skydata[,, l]
				        tdata <- sdata
	
					tiltm    <- sphere.to.dir.cos(cbind(rad(sim_par$sza), rad(180)))		
					sun_real <- deg(dir.cos.to.sphere(update.cosine(psi = rad(zentilt), phi = rad((azitilt + 180) %% 360), tiltm))) 
					if(sza == 1 || sza == 2) sun_real[2] <- (sun_real[2] + 180) %% 360		
					fact <- cos(rad(sun_real[1])) / cos(rad(sim_par$sza))				
					if(sun_real[1] > 90) fact <- 0							
		
					ida        <- which(azts == azitilt)

					intp   <- bicubic(x = rad(0:90), y = rad(0:360), z = tdata, x0 = gxy[, 1], y0 = gxy[, 2])$z	
					intp   <- intp * fx
					dim(intp) <- c(length(sx), length(sy))
					intp   <- intp[-length(sx), ] + intp[-1, ]
					intp   <- intp[, -length(sy)] + intp[, -1]
					tempdf <- sum(intp / 4)										
					tempdr <- dire[l, 2] * fact									

					idz  <- which(zets == zentilt)
					meas <- c(tempdf, tempdr)
					real <- dife[l, 2] + dire[l, 2]
		
					res_tilt_dif[idz, ida, l, sza, s] <- meas[1] / dife[l, 2]
					res_tilt_dir[idz, ida, l, sza, s] <- meas[2] / dire[l, 2]
					res_tilt_tot[idz, ida, l, sza, s] <- sum(meas) / real
					print(c(type[s], round(c(sim_par$sza, sim_par$lambda[l], zentilt, azitilt, res_tilt_tot[idz, ida, l, sza, s]), 3)))
				}
			}
		}
	}
	save(res_tilt_dif, file = "res/res_tilt_dif.Rdata")
	save(res_tilt_dir, file = "res/res_tilt_dir.Rdata")
	save(res_tilt_tot, file = "res/res_tilt_tot.Rdata")
}

#
# Mirroring of the results to avoid redundant calculations:
# Without shadow or BRDF, error depends only on relative azimuth.
#

res_tilt_dif[, c(6, 7, 8), , , ] <- res_tilt_dif[, c(4, 3, 2), , , ]
res_tilt_dir[, c(6, 7, 8), , , ] <- res_tilt_dir[, c(4, 3, 2), , , ]
res_tilt_tot[, c(6, 7, 8), , , ] <- res_tilt_tot[, c(4, 3, 2), , , ]

save(res_tilt_dif, file = "res/res_tilt_dif.Rdata")
save(res_tilt_dir, file = "res/res_tilt_dir.Rdata")
save(res_tilt_tot, file = "res/res_tilt_tot.Rdata")

#
# Extract error statistics. 
# Get the mean absolute percentage error and the standard deviation of the absolute errors
# for the visible range under recommended skies and  overcast conditions. Note average 
# errors are computed from zero tilt to each step of tilt. Dimensions are: 
# maximum tilt, sky (SZA from 20 to 60 = 1, and OC = 2)
#

dimnms <- list(c('znt03', 'znt06', 'znt09', 'znt12'), c('RSky', 'OC'))

error_tilt_avem <- array(NA, dim = c(length(zets), 2), dimnames = dimnms)
error_tilt_sdm  <- array(NA, dim = c(length(zets), 2), dimnames = dimnms)

temp <- array(1, dim = c(1, length(azts), length(lambda), 10, length(type)))		# Make array for the condition of no tilt. Under no Tilt, no BRDF and no Shadow, Eg estimated == Eg real
res_temp <- abind(temp, res_tilt_tot, along = 1)

for(i in 2:5) {
	temp <- abs(1 - res_temp[1:i, , 6:36, 4:8, , drop = F])
	error_tilt_avem[i-1, 1] <- mean(temp)
	temp <- abs(1 - res_temp[1:i, , 6:36, 1, , drop = F])
	error_tilt_avem[i-1, 2] <- mean(temp)
	temp <- abs(1 - res_temp[1:i, , 6:36, 4:8, , drop = F])
	error_tilt_sdm[i-1, 1] <- sd(temp)
	temp <- abs(1 - res_temp[1:i, , 6:36, 1, , drop = F])
	error_tilt_sdm[i-1, 2] <- sd(temp)
}

round(error_tilt_avem * 100, 1)
round(error_tilt_sdm * 100, 1)

#
# Graphical result:
#

library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
	"\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
	"\\documentclass[18pt]{article}")

col <- colorRampPalette(c("blue", "steelblue", "white", "orange", "red"))(256)

# At 550 nm for all tilts/azimuths:
for(sza in 2:10) {
	xmi <- res_tilt_tot[, 1:5, 21, sza, ]
	xm <- matrix(NA, ncol = 5, nrow = 4)
	for(i in 1:4) {
		for(j in 1:5) {
			xm[i, j] <- mean(xmi[i, j, ])
		}
	}
	fl <- paste0("plot/tilt_sza_", formatC((sza-2)*10, width = 2, flag = 0), ".tex")
	if(sza == 1) fl <- paste0("plot/tilt_sza_OC.pdf")
	tikz(fl, width = 7*1.2, height = 7, standAlone = TRUE, pointsize = 18,
		packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
		"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
		"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
		par(mar = c(5, 5, 2, 2))
		xlab = '$\\Delta \\phi_{\\text{N}}\\,(^{\\circ})$'
		ylab = '$\\theta_{\\text{N}}\\,(^{\\circ})$'
		plot(NA, xlim = c(1, 6), ylim = c(1, 5), xaxs = 'i', yaxs = 'i', xaxt = 'n', 
		yaxt = 'n', xlab = '', ylab = '')
		mtext(xlab, side = 1, line = 3, cex = 2)
		mtext(ylab, side = 2, line = 3, cex = 2)
		xrst <- matrix(.map2color(xm, col, limits = c(0.5, 1.5)), ncol = ncol(xm))[nrow(xm):1, ]
		rasterImage(xrst, xleft = 1, ybottom = 1, xright = 6, ytop = 5, 
		interpolate = F)
		axis(1, at = seq(1.5, 5.5, 1), labels = azts[1:5], cex.axis = 1.9)
		axis(2, at = seq(1.5, 4.5, 1), labels = zets, cex.axis = 1.9)
		box()
		if(sza == 1) {
			mtext("OC", side = 3, line = 0.5, cex = 2.5)
		} else {
			mtext(eval(substitute(expression('$\\theta_{\\text{s}}$'==ang*'$^{\\circ}$'), list(ang = (sza-2)*10))), side = 3, line = 0.5, cex = 2.5)
		}
	dev.off()
	tools::texi2pdf(fl)
	file.rename(paste0("tilt_sza_", formatC((sza-2)*10, width = 2, flag = 0), ".pdf"), paste0("plot/tilt_sza_", formatC((sza-2)*10, width = 2, flag = 0), ".pdf"))
	file.remove(paste0("tilt_sza_", formatC((sza-2)*10, width = 2, flag = 0), ".aux"), paste0("tilt_sza_", formatC((sza-2)*10, width = 2, flag = 0), ".log"))
}


tikz('plot/tilt_legend.tex', width = 7, height = 7*0.35, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	par(mar = c(5, 1, 1, 1))
	xlab <- '$E_{\\text{g}}^{\\text{ est.}}/E_{\\text{g}}^{\\text{{ real}}}$'
	plot(NA, xlim = c(0.5, 1.5), ylim = c(0, 1), xaxs = 'i', yaxs = 'i',
	yaxt = 'n', xlab = '', ylab = '', xaxt = 'n')
	mtext(xlab, side = 1, line = 3, cex = 2)
	axis(1, at = seq(0.6, 1.4, 0.2), labels = TRUE, cex.axis = 1.9)
	xrst <- t(as.raster(.map2color(seq(0.5, 1.5, length.out = length(col)), col)))
	rasterImage(xrst, xleft = 0.5, ybottom = 0, xright = 1.5, ytop = 1, interpolate = T)
dev.off()
tools::texi2pdf('plot/tilt_legend.tex')
file.rename('tilt_legend.pdf', 'plot/tilt_legend.pdf')
file.remove('tilt_legend.aux', 'tilt_legend.log')

#
# Combined plot with legend:
#

tikz('plot/error_tilt_sza50_legend.tex', width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	layout(matrix(c(1,2),nrow=2), heights = c(3,1))
	sza = 7
	xmi <- res_tilt_tot[-5, 1:5, 21, sza, ]
	xm <- matrix(NA, ncol = 5, nrow = 4)
	for(i in 1:4) {
		for(j in 1:5) {
			xm[i, j] <- mean(xmi[i, j, ])
		}
	}
	par(mar = c(5, 5, 4, 2.3))
	xlab = '$\\Delta \\phi_{\\text{N}}\\,(^{\\circ})$'
	ylab = '$\\theta_{\\text{N}}\\,(^{\\circ})$'
	plot(NA, xlim = c(1, 6), ylim = c(1, 5), xaxs = 'i', yaxs = 'i', xaxt = 'n', 
		yaxt = 'n', xlab = '', ylab = '')
	mtext(xlab, side = 1, line = 3, cex = 2)
	mtext(ylab, side = 2, line = 3, cex = 2)
	xrst <- matrix(.map2color(xm, col, limits = c(0.5, 1.5)), ncol = ncol(xm))[nrow(xm):1, ]
	rasterImage(xrst, xleft = 1, ybottom = 1, xright = 6, ytop = 5, interpolate = F)
	axis(1, at = seq(1.5, 5.5, 1), labels = azts[1:5], cex.axis = 1.9)
	axis(2, at = seq(1.5, 4.5, 1), labels = zets, cex.axis = 1.9)
	box()
	mtext(eval(substitute(expression('$\\theta_{\\text{s}}$'==ang*'$^{\\circ}$'), list(ang = (sza-2)*10))), side = 3, line = 0.5, cex = 2.5)
	par(mar = c(5, 5, 0.1, 2.3))
	xlab <- '$E_{\\text{g}}^{\\text{ est.}}/E_{\\text{g}}^{\\text{{ real}}}$'
	plot(NA, xlim = c(0.7, 1.3), ylim = c(0, 1), xaxs = 'i', yaxs = 'i',
	yaxt = 'n', xlab = '', ylab = "", xaxt = 'n')
	mtext(xlab, side = 1, line = 3, cex = 2)
	axis(1, at = seq(0.7, 1.3, 0.1), labels = TRUE, cex.axis = 1.9)
	xrst <- t(as.raster(.map2color(seq(0.5, 1.5, length.out = length(col)), col)))
	rasterImage(xrst, xleft = 0.5, ybottom = 0, xright = 1.5, ytop = 1, interpolate = T)
	box()
dev.off()
tools::texi2pdf('plot/error_tilt_sza50_legend.tex')
file.rename('error_tilt_sza50_legend.pdf', 'plot/error_tilt_sza50_legend.pdf')
file.remove('error_tilt_sza50_legend.aux', 'error_tilt_sza50_legend.log')

#
# Mean tilt error per SZA, for all tilts and azimuths considered:
#

mean_tilt_error <- matrix(NA, ncol = 58, nrow = 10)
for(sza in 1:10) {
	for(lbd in 1:58) {
		xm <- res_tilt_tot[, , lbd, sza, ]
		mean_tilt_error[sza, lbd] <- mean(abs(1 - xm))
	}
}

cols <- rev(rainbow(6, start = 0, end = 0.8))
cols[5] <- "#FFB84D"

tikz('plot/error_tilt.tex', width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	cols <- rev(rainbow(6, start = 0, end = 0.8))
	par(mar = c(5, 5, 3, 2))
	xlab <- '$\\lambda$ (nm)'
	ylab = 'Average $\\epsilon$ (\\%)'
	plot(lambda[1:58], mean_tilt_error[1, ]*100, type = "l", ylim = range(mean_tilt_error)*100, col = cols[1], xlab = '', ylab = '', lwd = 3, yaxt = 'n', xaxt = 'n')
	axis(2, at = seq(0, 40, 10), labels = T, cex.axis = 1.9)
	axis(1, at = seq(400, 1000, 100), labels = T, cex.axis = 1.9)
	mtext(xlab, side = 1, line = 3, cex = 2)
	mtext(ylab, side = 2, line = 3, cex = 2)
	lines(lambda[1:58], mean_tilt_error[2, ]*100, col = cols[2], lwd = 3)
	lines(lambda[1:58], mean_tilt_error[4, ]*100, col = cols[3], lwd = 3)
	lines(lambda[1:58], mean_tilt_error[6, ]*100, col = cols[4], lwd = 3)
	lines(lambda[1:58], mean_tilt_error[8, ]*100, col = cols[5], lwd = 3)
	lines(lambda[1:58], mean_tilt_error[10, ]*100, col = cols[6], lwd = 3)
	legend(720, 30, c('OC', '$\\theta_{\\text{s}}=0^{\\circ}$', '$\\theta_{\\text{s}}=20^{\\circ}$', '$\\theta_{\\text{s}}=40^{\\circ}$', '$\\theta_{\\text{s}}=60^{\\circ}$', 
		'$\\theta_{\\text{s}}=80^{\\circ}$'), lty = 1, col = cols, bty = "n", cex = 1.2, lwd = 3, y.intersp = 0.5)
dev.off()
tools::texi2pdf('plot/error_tilt.tex')
file.rename('error_tilt.pdf', 'plot/error_tilt.pdf')
file.remove('error_tilt.aux', 'error_tilt.log')

