
################################################################################
# Aerosol models:
#
# Aerosol models were retrieved from OPAC by using 4 different aerosol models 
# (maritime clean, maritime polluted, continental clean and continental 
# polluted), with 3 humidity levels (50, 80 and 95%) and two AOT_550 levels 
# (0.1 and 0.2). The input files used with the OPAC code are stored in the anc 
# folder. To model different aerosol loads and make them match to the specified 
# levels (0.1 and 0.2 at 550 nm), the scale height was changed from the default 
# values for the mixed layer. The mixed layer height was kept at 2 km and all 
# other vertical distribution profiles kept at default values.
#

aer.fls <-  paste0("anc/opac_out/", rep(c("maritime", "continental"), 
	each = 12), "_", rep(rep(c("clean", "polluted"), each = 6), 2), "_", 
	rep(rep(c(50, 80, 95), each = 2), 2), "_", rep(c("01", "02"), 12), 
	".out")

aer.mod <-  paste0(rep(c("maritime", "continental"), each = 12), "_", 
	rep(rep(c("c", "p"), each = 6), 2), rep(rep(c(50, 80, 95), 
	each = 2), 2), "_", rep(c("01", "02"), 12))

rh <- rep(rep(c(50, 80, 95), each = 2), 4)

for(aer in 1:length(aer.fls)) {
	assign(aer.mod[aer], read_opac(aer.fls[aer], rh = rh[aer]))
}

# Check AOT at 550:
taus <- numeric(length(aer.mod))
for(i in 1:length(aer.mod)) {
	taus[i] <- get(aer.mod[i])[[1]][5, 6]
}
data.frame(aer.mod, round(taus, 1)) 

# Normalize the volume scattering function to the phase function:
# phase function in atmospheric optics is multiplied by 4pi, so that it has no 
# units. So the normalization condition takes that into consideration and is:
# 1/4pi * integral_phi integral_theta beta_tilde * sin(theta) dtheta dphi
# = 1/2 * integral_theta beta_tilde * sin(theta) dtheta
# Therefore, to normalize the phase function, the integral is also divided by 2.

fun <- function(x, data)
 {
  temp <- approx(x = data[, 1], y = data[, 2], xout = x)$y
  return(temp * sin(x))
 }

for(aer in 1:length(aer.mod)) {
	aer.int <- get(aer.mod[aer])
	for(lbd in 1:ncol(aer.int[[3]])) {
		aer.int[[3]][, lbd] <- aer.int[[3]][, lbd] / (integrate(f = fun, 
			lower = 0, upper = pi, data = cbind(rad(aer.int[[2]]), 
			aer.int[[3]][, lbd]), subdivisions = 1000)$value / 2)
	}
	assign(aer.mod[aer], aer.int)
}

################################################################################
# RT model run:
# 
# sim_par controls the simulation and has a fixed condition for all simulations 
# except aerosol model. RT is run for all angles in 1º spacing and wavelengths 
# from 350 to 1000 in 10 nm steps and Sun zenith angles from 0 to 80º in 10º 
# steps.
# 

library(abind)

type <- sub(".out", "", sub("anc/opac_out/", "", aer.fls))

for(aer in 1:length(aer.mod)) {
	lambda  <- seq(350, 1000, 10)
	aer_m   <- get(aer.mod[aer])

	sim_par <- simpar(
		lambda  = lambda,
		h       = rh[aer],
		rho     = approxfun(n_w[, 1], fref_sw_diff)(lambda) + approxfun(coastal[, 1], 
				coastal[, 2], rule = 2)(lambda),
		sza     = 0,
		tau_aer = approx(aer_m$opp[, 1], aer_m$opp[, 6], xout = lambda)$y,
		w0_a    = approx(aer_m$opp[, 1], aer_m$opp[, 5], xout = lambda)$y,
		p_a     = cbind(rad(aer_m[[2]]), aer_m[[3]][, 4]),
		f_0     = f0(lambda, doy = 1, dataset = "thuillier")
	)

	for(sza in seq(0, 80, 10)) {
		sim_par$sza <- sza
		skylight <- skyrad(sim_par)
		skylight <- abind(skylight, skylight[, 180:1, ], along = 2)
		diffuse  <- numeric(length(lambda))

		cat("\nIntegrating sky radiance to diffuse irradiance...")
		for(lbd in 1:length(lambda)) {
				intp   <- bicubic(x = rad(0:90), y = rad(0:360), z = skylight[,,lbd], x0 = gxy[, 1], y0 = gxy[, 2])$z
				intp   <- intp * fx
				dim(intp) <- c(length(sx), length(sy))
				intp   <- intp[-length(sx), ] + intp[-1, ]
				intp   <- intp[, -length(sy)] + intp[, -1]
				diffuse[lbd] <- sum(intp / 4)
		}
		cat(" done!\n")

		cat("Computing direct irradiance...")
		Es <- es_m(sim_par$lambda, sim_par$tau_aer, sim_par$sza, sim_par$f_0, 
			sim_par$c_w, sim_par$c_o3, sim_par$co2, sim_par$p, sim_par$h, 
			sim_par$t)
		cat(" done!\n")

		cat("Writing to disk...")
		simpar.fl <- file(paste0("sim/skyrad/", type[aer], "_sza_", formatC(sim_par$sza, 
				width = 2, flag = 0), "_simpar.txt"), "w")
		for(i in 1:14) {
			writeLines(c(names(sim_par)[i], paste0(sim_par[[i]], 
				collapse = ",")), con = simpar.fl)
		}
		writeLines(names(sim_par)[15], con = simpar.fl)
		write.table(sim_par$p_a, file = simpar.fl, append = T, row.names = F, 
			col.names = F)
		close(simpar.fl)

		Ed <- cbind(lambda, diffuse)
		colnames(Ed) <- c("wavelength", "diffuse_irrad")
		write.table(format(Ed, digits = 4, scientific = T), 
		paste0("sim/skyrad/", type[aer], "_sza_", formatC(sim_par$sza, width = 2, flag = 0), 
			"_diffuse.txt"), row.names = F, quote = F) 
		Es <- cbind(lambda, Es)
		colnames(Es) <- c("wavelength", "direct_irrad")
		write.table(format(Es, digits = 4, scientific = T), paste0("sim/skyrad/", type[aer], 
			"_sza_", formatC(sim_par$sza, width = 2, flag = 0), 
			"_direct.txt"), row.names = F, quote = F)

		longtab <- skylight
		dim(longtab) <- c(length(sim_par$zenith) * (length(sim_par$razimuth) + 180), 
				  length(sim_par$lambda))
		longtab <- cbind(expand.grid(sim_par$zenith, c(sim_par$razimuth, sim_par$razimuth[180:1])), 
				longtab)
		colnames(longtab) <- c("tetha", "phi", paste0("l", sim_par$lambda))
		write.table(format(longtab, digits = 4, scientific = T), paste0("sim/skyrad/", 
			type[aer], "_sza_", formatC(sim_par$sza, width = 2, flag = 0), 
			"_sky.txt"), row.names = F, quote = F)
		cat(" done!\n\nFinished!\n\n")
	}
}

# Cardioidal distribution, from Gordon 1985:
cardioidal <- function(Eg, thetai) {
	(3/(7*pi)) * Eg * (1 + 2 * cos(thetai))
}

plot(0:90, cardioidal(1, rad(0:90)), type = "l")

# Direct irradiance is zero and diffuse irradiance is 100:

skylight <- matrix(cardioidal(100, rad(0:90)), ncol = 1)
skylight <- skylight[, rep(1, 361)]
dim(skylight) <- c(91*361, 1)
skylight <- cbind(expand.grid(0:90, 0:360), skylight)
skylight <- skylight[, c(1:2, rep(3, length(lambda)))]
colnames(skylight) <- c("tetha", "phi", paste0("l", lambda))
write.table(format(skylight, digits = 4, scientific = T), "sim/skyrad/overcast_sky.txt", row.names = F, quote = F)

Ed <- cbind(lambda, 100)
colnames(Ed) <- c("wavelength", "diffuse_irrad")
write.table(format(Ed, digits = 4, scientific = T), 
"sim/skyrad/overcast_diffuse.txt", row.names = F, quote = F) 
Es <- cbind(lambda, 0)
colnames(Es) <- c("wavelength", "direct_irrad")
write.table(format(Es, digits = 4, scientific = T), 
"sim/skyrad/overcast_direct.txt", row.names = F, quote = F)

file.copy("sim/skyrad/maritime_clean_50_01_sza_00_simpar.txt", "sim/skyrad/overcast_simpar.txt")

#
# Make plots:
#

library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
	"\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
	"\\documentclass[18pt]{article}")

fls   <- paste0("sim/skyrad/", type[seq(2, 24, 2)], rep("_sza_40_sky.txt", 12))
ylabs <- rep(c("Maritime clean", "Maritime polluted", "Continental clean", "Continental polluted"), each = 3)
xlabs <- rep(c("RH: 50 \\%", "RH: 80 \\%", "RH: 95 \\%"), 4)
for(i in 1:length(fls)) {
	tikz(paste0("plot/skyrad_", type[seq(2, 24, 2)][i], ".tex"), width = 7, height = 7, standAlone = TRUE, pointsize = 18,
		packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
		"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
		"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
		temp <- as.matrix(read.table(fls[i], header = T)[, -c(1, 2)])
		dim(temp) <- c(91, 361, 66)
		plot_skyr_latex(temp[, , 21])
		if(i == 1 | i == 4 | i == 7 | i == 10) mtext(ylabs[i], side = 2, line = 1, cex = 3.3)
		if(i == 1 | i == 2 | i == 3) mtext(xlabs[i], side = 3, line = 1, cex = 3.0)
	dev.off()
	tools::texi2pdf(paste0("plot/skyrad_", type[seq(2, 24, 2)][i], ".tex"))
	file.rename(paste0("skyrad_", type[seq(2, 24, 2)][i], ".pdf"), paste0("plot/skyrad_", type[seq(2, 24, 2)][i], ".pdf"))
	file.remove(paste0("skyrad_", type[seq(2, 24, 2)][i], ".aux"), paste0("skyrad_", type[seq(2, 24, 2)][i], ".log"))
}

.map2color <- function(x, pal, limits){
    if(missing(limits)) limits <- range(x, na.rm = T)
    pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal)+1), all.inside = TRUE)]
}

cols <- colorRampPalette(c("steelblue4", "steelblue", "white", "orange", "red"))(256)

tikz("plot/skyrad_legend.tex", width = 7, height = 2.3, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	par(mar = c(5, 2, 2, 2))
	plot(NA, xlim = log(c(10, 1000), 10), ylim = c(0, 1), xaxs = 'i', yaxs = 'i', xaxt = 'n', 
	yaxt = 'n', xlab = "", ylab = "")
	axis(1, at = log(c(10, 30, 100, 300, 1000), 10), labels = c("$<$10", "30", "100", "300", "$>$1000"), cex.axis = 2)
	xrst <- t(as.raster(.map2color(seq(1, 3, length.out = length(cols)), cols)))
	rasterImage(xrst, xleft = 1, ybottom = 0, xright = 3, ytop = 1, interpolate = T)
	mtext('$L_{\\text{d}}(\\theta_{\\text{i}}, \\phi_{\\text{i}}, 550)$ (mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$)', side = 1, line = 3, cex = 2.2)
	box()
dev.off()
tools::texi2pdf("plot/skyrad_legend.tex")
file.rename("skyrad_legend.pdf", "plot/skyrad_legend.pdf")
file.remove("skyrad_legend.aux", "skyrad_legend.log")


tikz("plot/aer_beta.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	xlab <- '$\\psi\\,(^{\\circ})$'
	ylab <- '$\\tilde{\\beta}_{\\text{a}}(\\psi, 550)$ (unitless)'
	par(mar = c(5, 5, 3, 2))
	aer.t <- get(aer.mod[1])
	plot(aer.t[[2]], aer.t[[3]][, 5], type = "l", log = "xy", col = "blue", 
		lty = 1, lwd = 4, ylim = c(0.05, 600), xlab = '', ylab = '', 
		xaxt = "n", yaxt = "n")
	axis(2, at = c(0.05, 0.5, 5, 50, 500), labels = c("0.05", "0.5", "5", "50", "500"), cex.axis = 1.9)
	axis(1, at = c(0.1, 1, 10, 45, 180), labels = c("0.1", "1", "10", "45", "180"), cex.axis = 1.9)
	mtext(xlab, side = 1, line = 3, cex = 2)
	mtext(ylab, side = 2, line = 3, cex = 2)

	aer.t <- get(aer.mod[3])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "blue", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[5])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "blue", lty = 3, lwd = 4)

	aer.t <- get(aer.mod[7])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "orange", lty = 1, lwd = 4)
	aer.t <- get(aer.mod[9])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "orange", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[11])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "orange", lty = 3, lwd = 4)

	aer.t <- get(aer.mod[13])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "deepskyblue", lty = 1, lwd = 4)
	aer.t <- get(aer.mod[15])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "deepskyblue", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[17])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "deepskyblue", lty = 3, lwd = 4)

	aer.t <- get(aer.mod[19])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "red", lty = 1, lwd = 4)
	aer.t <- get(aer.mod[21])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "red", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[23])
	lines(aer.t[[2]], aer.t[[3]][, 5], col = "red", lty = 3, lwd = 4)
	legend("topright", c("RH:", "50 \\%", "80 \\%", "95 \\%"), lty = c(NA, 1, 2, 3), 
		bty = "n", lwd = 4, cex = 1.5, y.intersp = 0.5, x.intersp = 0.5)
	legend("bottomleft", c("Model:", "Maritime clean", "Maritime polluted", 
		"Continental clean", "Continental polluted"), lty = c(NA, 1, 1, 1, 1), 
		col = c(NA, "blue", "orange", "deepskyblue", "red"), bty = "n", lwd = 4, cex = 1.5, y.intersp = 0.5, x.intersp = 0.5)
dev.off()
tools::texi2pdf("plot/aer_beta.tex")
file.rename("aer_beta.pdf", "plot/aer_beta.pdf")
file.remove("aer_beta.aux", "aer_beta.log")

tikz("plot/aer_w0.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	xlab <- '$\\lambda$ (nm)'
	ylab <- '$\\omega_0$ (unitless)'
	par(mar = c(5, 5, 3, 2))
	aer.t <- get(aer.mod[1])[[1]][1:12, ]
	plot(aer.t[, 1], aer.t[, 5], col = "blue", lty = 1, lwd = 4, type = "l", 
		xlim = c(350, 1000), ylim = c(0.75, 1), xlab = '', ylab = '', 
		xaxt = "n", yaxt = "n")
	axis(2, at = seq(0.75, 1, 0.05), labels = T, cex.axis = 1.9)
	axis(1, at = seq(400, 1000, 100), labels = T, cex.axis = 1.9)
	mtext(xlab, side = 1, line = 3, cex = 2)
	mtext(ylab, side = 2, line = 3, cex = 2)

	aer.t <- get(aer.mod[3])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "blue", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[5])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "blue", lty = 3, lwd = 4)

	aer.t <- get(aer.mod[7])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "orange", lty = 1, lwd = 4)
	aer.t <- get(aer.mod[9])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "orange", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[11])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "orange", lty = 3, lwd = 4)

	aer.t <- get(aer.mod[13])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "deepskyblue", lty = 1, lwd = 4)
	aer.t <- get(aer.mod[15])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "deepskyblue", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[17])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "deepskyblue", lty = 3, lwd = 4)

	aer.t <- get(aer.mod[19])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "red", lty = 1, lwd = 4)
	aer.t <- get(aer.mod[21])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "red", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[23])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 5], col = "red", lty = 3, lwd = 4)
dev.off()
tools::texi2pdf("plot/aer_w0.tex")
file.rename("aer_w0.pdf", "plot/aer_w0.pdf")
file.remove("aer_w0.aux", "aer_w0.log")

tikz("plot/aer_tau.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	xlab <- '$\\lambda$ (nm)'
	ylab <- '$\\tau_{\\text{a}}(\\lambda) / \\tau_{\\text{a}}(550)$'
	par(mar = c(5, 5, 3, 2))
	aer.t <- get(aer.mod[1])[[1]][1:12, ]
	plot(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "blue", lty = 1, lwd = 4, 
		type = "l", xlim = c(350, 1000), ylim = c(0.4, 1.7), xlab = '', ylab = '',
	xaxt = "n", yaxt = "n")
	axis(2, at = seq(0.4, 1.6, 0.2), labels = T, cex.axis = 1.9)
	axis(1, at = seq(400, 1000, 100), labels = T, cex.axis = 1.9)
	mtext(xlab, side = 1, line = 3, cex = 2)
	mtext(ylab, side = 2, line = 3, cex = 2)

	aer.t <- get(aer.mod[3])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "blue", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[5])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "blue", lty = 3, lwd = 4)
	
	aer.t <- get(aer.mod[7])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "orange", lty = 1, lwd = 4)
	aer.t <- get(aer.mod[9])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "orange", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[11])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "orange", lty = 3, lwd = 4)
	
	aer.t <- get(aer.mod[13])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "deepskyblue", lty = 1, lwd = 4)
	aer.t <- get(aer.mod[15])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "deepskyblue", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[17])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "deepskyblue", lty = 3, lwd = 4)
	
	aer.t <- get(aer.mod[19])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "red", lty = 1, lwd = 4)
	aer.t <- get(aer.mod[21])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "red", lty = 2, lwd = 4)
	aer.t <- get(aer.mod[23])[[1]][1:12, ]
	lines(aer.t[, 1], aer.t[, 6] / aer.t[5, 6], col = "red", lty = 3, lwd = 4)
dev.off()
tools::texi2pdf("plot/aer_tau.tex")
file.rename("aer_tau.pdf", "plot/aer_tau.pdf")
file.remove("aer_tau.aux", "aer_tau.log")

fls.dir <- paste0("sim/skyrad/", type[seq(2, 24, 2)], rep("_sza_40_direct.txt", 12))
fls.dif <- paste0("sim/skyrad/", type[seq(2, 24, 2)], rep("_sza_40_diffuse.txt", 12))

tikz("plot/aer_diff.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	par(mar = c(5, 5, 3, 2))
	xlab <- '$\\lambda$ (nm)'
	ylab <- '$E_{\\text{d}}/E_{\\text{g}}$'
	cols <- rep(c("blue", "orange", "deepskyblue", "red"), each = 3)
	ltys <- rep(c(1, 2, 3), 4)
	plot(NA, xlim = range(lambda), ylim = c(0.05, 0.5), xlab = '', ylab = '',
	xaxt = "n", yaxt = "n")
	axis(2, at = seq(0.1, 0.5, 0.1), labels = T, cex.axis = 1.9)
	axis(1, at = seq(400, 1000, 100), labels = T, cex.axis = 1.9)
	mtext(xlab, side = 1, line = 3, cex = 2)
	mtext(ylab, side = 2, line = 3, cex = 2)
	for(i in 1:12) {
		dir <- read.table(fls.dir[i], header = T)[, 2]
		dif <- read.table(fls.dif[i], header = T)[, 2]
		lines(lambda, dif/(dir+dif), col = cols[i], lty = ltys[i], lwd = 4)
	}
	legend("topright", c('$\\theta_{\\text{s}} = 40^{\\circ}$', '$\\tau_{\\text{a}}(550) = 0.2$'), bty = "n", cex = 2, y.intersp = 0.5)
dev.off()
tools::texi2pdf("plot/aer_diff.tex")
file.rename("aer_diff.pdf", "plot/aer_diff.pdf")
file.remove("aer_diff.aux", "aer_diff.log")



