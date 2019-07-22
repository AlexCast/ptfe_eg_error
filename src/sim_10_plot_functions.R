

#
# Visualize BRDF
#

plot_brdf <- function(x, viewg = NULL, ste = F, ort = T, leg = T) {
	require("raster")
	if(is.null(viewg)) viewg <- c(NA, NA)
	rbrdf  <- raster(x, xmn = -180.5, xmx = 180.5, ymn = -0.5, ymx = 90.5, 
		crs = crs.wgs)
	base   <- raster(nrows = 90, ncols = 360, xmn = -180, xmx = 180, ymn = 0, 
		ymx = 90, crs = crs.wgs)
	base.c <- crop(base, extent(-180, 180, 0, 90))
	rbrdf  <- resample(rbrdf, base)
	if(ste) {
		rbrdf <- projectRaster(crop(rbrdf, extent(-180, 180, 0, 90)), 
			crs = crs.ste)
		plot(rbrdf, bty = "n", box = FALSE, zlim = c(0.20, 0.45), 
			axes = F, col = heat.colors(255), axis.args = list(at = 
			c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45), labels = TRUE, 
			cex.axis = 0.8), legend.args = list(text = 
			eval(substitute(expression("BRDF"*(theta[i]*", "*
			phi[i]*", "*a*", "*b)), list(a = viewg[1], b = viewg[2]))), 
			side = 4, line = 2.9, cex = 0.8), legend = leg)
		plot(zlin_ste, add = TRUE, col = "grey", lwd = 2)
		lines(c(0, 0), c(-12350000, 12350000), lwd = 2, col = "grey")
		lines(c(-12350000, 12350000), c(0, 0), lwd = 2, col = "grey")
	} else if(ort) {
		rbrdf <- projectRaster(crop(rbrdf, extent(-180, 180, 0, 90)), 
			crs = crs.ort)
		plot(rbrdf, bty = "n", box = FALSE, zlim = c(0.20, 0.45), 
			axes = F, col = heat.colors(255), axis.args = list(at = 
			c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45), labels = TRUE, 
			cex.axis = 0.8), legend.args = list(text = 
			eval(substitute(expression("BRDF"*(theta[i]*", "*
			phi[i]*", "*a*", "*b)), list(a = viewg[1], b = viewg[2]))), 
			side = 4, line = 2.9, cex = 0.8), legend = leg)
		plot(zlin_ort, add = TRUE, col = "grey45", lwd = 2)
		lines(c(0, 0), c(-6.35E6, 6.35E6), lwd = 2, col = "grey45")
		lines(c(-6.35E6, 6.35E6), c(0,0), lwd = 2, col = "grey45")
		text(x = 3.2E6, y = 3.2E6, labels = '40°', col = 'grey45')
		text(x = 1.45E6, y = 1.45E6, labels = '15°', col = 'grey45')
		text(x = 0, y = -6.7E6, labels = '0°', cex = 1, xpd = T, col = 'grey45')
		text(x = -7.0E6, y = 0, labels = '270°', cex = 1, xpd = T, col = 'grey45')
		text(x = 0, y = 6.7E6, labels = '180°', cex = 1, xpd = T, col = 'grey45')
		text(x = 6.9E6, y = 0, labels = '90°', cex = 1, xpd = T, col = 'grey45')
	} else {
		plot(rbrdf, col = heat.colors(255))
	}
 }

plot_brdf_latex <- function(x, viewg = NULL, leg = T) {
	require("raster")
	if(is.null(viewg)) viewg <- c(NA, NA)
	rbrdf  <- raster(x, xmn = -180.5, xmx = 180.5, ymn = -0.5, ymx = 90.5, 
		crs = crs.wgs)
	base   <- raster(nrows = 90, ncols = 360, xmn = -180, xmx = 180, ymn = 0, 
		ymx = 90, crs = crs.wgs)
	base.c <- crop(base, extent(-180, 180, 0, 90))
	rbrdf  <- resample(rbrdf, base)
	rbrdf  <- projectRaster(crop(rbrdf, extent(-180, 180, 0, 90)), crs = crs.ort)
	plot(rbrdf, bty = "n", box = FALSE, zlim = c(0.20, 0.45), 
		axes = F, col = heat.colors(255), axis.args = list(at = 
		c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45), labels = 
		c('0.20', '0.25', '0.30', '0.35', '0.40', '$>0.45$'), 
		cex.axis = 1), legend.args = list(text = 
		eval(substitute(expression('$f_{\\text{r}}(\\theta_{\\text{i}}, \\phi_{\\text{i}}, $'*a*'$^{\\circ}, $'*b*'$^{\\circ}, 550)$'~'sr$^{-1}$'), list(a = viewg[1], b = viewg[2]))), 
		side = 4, line = 2.5, cex = 1.3), legend = leg)
	plot(zlin_ort, add = TRUE, col = "grey45", lwd = 4)
	lines(c(0, 0), c(-6.35E6, 6.35E6), lwd = 4, col = "grey45")
	lines(c(-6.35E6, 6.35E6), c(0,0), lwd = 4, col = "grey45")
	text(x = 4.92E6, y = 4.92E6, labels = '$90^{\\circ}$', col = 'grey45', cex = 1.5)
	text(x = 3.27E6, y = 3.27E6, labels = '$40^{\\circ}$', col = 'grey45', cex = 1.5)
	text(x = 1.52E6, y = 1.52E6, labels = '$15^{\\circ}$', col = 'grey45', cex = 1.5)
	text(x = 0.18E6, y = -6.9E6, labels = '$0^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
	text(x = -7.2E6, y = 0, labels = '$270^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
	text(x = 0.18E6, y = 6.7E6, labels = '$180^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
	text(x = 7.1E6, y = 0, labels = '$90^{\\circ}$', cex = 1.5, xpd = T, col = 'grey45')
 }

plot_brdf_latex(brdf_ptfe99[[geom]][ , , 21], viewg = c(0, 90))

#
# Visualize sky radiance
#

plot_skyr <- function(x, ste = T, cont = T, leg = T, main = '') {
	require("raster")
	radaf  <- raster(x, xmn = -180.5, xmx = 180.5, ymn = -0.5, ymx = 90.5, 
			crs = crs.wgs)
	base   <- raster(nrows = 90, ncols = 360, xmn = -180, xmx = 180, ymn = 0, 
			ymx = 90, crs = crs.wgs)
	base.c <- crop(base, extent(-180, 180, 0, 90))
	radaf  <- resample(radaf, base)
	cols <- colorRampPalette(c("steelblue4", "steelblue", "white", "orange", "red"))(256)
	if(ste) {
		radaf <- projectRaster(crop(radaf, extent(-180, 180, 0, 90)), 
			crs = crs.ste)
		radaf[radaf > 1000] <- 1000
		radaf[radaf < 10] <- 10
		plot(log(radaf, 10), bty = "n", box = FALSE, zlim = c(log(10, 
			10), log(1000, 10)), axes = F, axis.args = list(at = 
			log(c(10, 30, 100, 300, 1000), 10), labels = c("<10", 
			"30", "100", "300", ">1000"), cex.axis = 0.8), 
			legend.args = list(text=expression(L[550]~(mW~m^2~nm^-1~sr^-1)), 
			side = 4, line = 2.9, cex = 0.8), legend = leg, col = cols, main = main)
		if(cont) contour(log(radaf, 10), add = T, levels = log(c(50, 
			100, 300, 600, 800, 1000), 10), labels = c(50, 100, 300, 
			600, 800, 1000))
		plot(zlin_ste, add = TRUE, col = "grey", lwd = 2)
		lines(c(0, 0), c(-12350000, 12350000), lwd = 2, col = "grey")
		lines(c(-12350000, 12350000), c(0, 0), lwd = 2, col = "grey")
	} else if(ort) {
		rbrdf <- projectRaster(crop(rbrdf, extent(-180, 180, 0, 90)), 
			crs = crs.ort)
		radaf[radaf > 1000] <- 1000
		radaf[radaf < 10] <- 10
		plot(log(radaf, 10), bty = "n", box = FALSE, zlim = c(log(10, 
			10), log(1000, 10)), axes = F, axis.args = list(at = 
			log(c(10, 30, 100, 300, 1000), 10), labels = c("<10", 
			"30", "100", "300", ">1000"), cex.axis = 0.8), 
			legend.args = list(text=expression(L[550]~(mW~m^2~nm^-1~sr^-1)), 
			side = 4, line = 2.9, cex = 0.8), legend = leg, col = cols, main = main)
		if(cont) contour(log(radaf, 10), add = T, levels = log(c(50, 
			100, 300, 600, 800, 1000), 10), labels = c(50, 100, 300, 
			600, 800, 1000))
		plot(zlin_ort, add = TRUE, col = "grey45", lwd = 2)
		lines(c(0, 0), c(-6.35E6, 6.35E6), lwd = 2, col = "grey45")
		lines(c(-6.35E6, 6.35E6), c(0,0), lwd = 2, col = "grey45")
		text(x = 3.2E6, y = 3.2E6, labels = '40°', col = 'grey45')
		text(x = 1.45E6, y = 1.45E6, labels = '15°', col = 'grey45')
		text(x = 0, y = -6.7E6, labels = '0°', cex = 1, xpd = T, col = 'grey45')
		text(x = -7.0E6, y = 0, labels = '270°', cex = 1, xpd = T, col = 'grey45')
		text(x = 0, y = 6.7E6, labels = '180°', cex = 1, xpd = T, col = 'grey45')
		text(x = 6.9E6, y = 0, labels = '90°', cex = 1, xpd = T, col = 'grey45')
	} else {
		plot(radaf, col = cols, main = main)
	}
}

plot_skyr_latex <- function(x, cont = T, leg = T, main = '') {
	require("raster")
	radaf  <- raster(x, xmn = -180.5, xmx = 180.5, ymn = -0.5, ymx = 90.5, 
			crs = crs.wgs)
	base   <- raster(nrows = 90, ncols = 360, xmn = -180, xmx = 180, ymn = 0, 
			ymx = 90, crs = crs.wgs)
	base.c <- crop(base, extent(-180, 180, 0, 90))
	radaf  <- resample(radaf, base)
	cols <- colorRampPalette(c("steelblue4", "steelblue", "white", "orange", "red"))(256)
	radaf <- projectRaster(crop(radaf, extent(-180, 180, 0, 90)), crs = crs.ort)
	radaf[radaf > 1000] <- 1000
	radaf[radaf < 10] <- 10
	plot(log(radaf, 10), bty = "n", box = FALSE, zlim = c(log(10, 
		10), log(1000, 10)), axes = F, axis.args = list(at = 
		log(c(10, 30, 100, 300, 1000), 10), labels = c("$<10$", 
		"30", "100", "300", "$>1000$"), cex.axis = 1), 
		legend.args = list(text = '$L_{\\text{d}}(\\theta_{\\text{i}}, \\phi_{\\text{i}}, 550)$ (mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$)', 
		side = 4, line = 2.5, cex = 1.3), legend = leg, col = cols, main = main)
	plot(zlin_ort, add = TRUE, col = "grey45", lwd = 4)
	lines(c(0, 0), c(-6.35E6, 6.35E6), lwd = 4, col = "grey45")
	lines(c(-6.35E6, 6.35E6), c(0,0), lwd = 4, col = "grey45")
	text(x = 4.85E6, y = 4.85E6, labels = '$90^{\\circ}$', col = 'grey45')
	text(x = 3.2E6, y = 3.2E6, labels = '$40^{\\circ}$', col = 'grey45')
	text(x = 1.45E6, y = 1.45E6, labels = '$15^{\\circ}$', col = 'grey45')
	text(x = 0, y = -6.7E6, labels = '$0^{\\circ}$', cex = 1, xpd = T, col = 'grey45')
	text(x = -7.0E6, y = 0, labels = '$270^{\\circ}$', cex = 1, xpd = T, col = 'grey45')
	text(x = 0, y = 6.7E6, labels = '$180^{\\circ}$', cex = 1, xpd = T, col = 'grey45')
	text(x = 6.9E6, y = 0, labels = '$90^{\\circ}$', cex = 1, xpd = T, col = 'grey45')
	if(cont) contour(log(radaf, 10), add = T, levels = log(c(50, 
		100, 300, 600, 800, 1000), 10), labels = c(50, 100, 300, 
		600, 800, 1000), labcex = 1.3, lwd = 2)
}

