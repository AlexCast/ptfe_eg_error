# Eg estimation from PTFE plaques
#
# Analysis part 2 - Pre-compute the BRDF for each plaque, view geometry and tilt condition
#
# The normalized BRDF for each plaque and view geometry is calculated. Additionally the 
# normalized BRDF for each plaque and view geometry under tilt (tilt as modeled here 
# result in change of the view geometry) is also calculated. The hemispherical-directional
# reflectance is also calculated both for the plaque level and tilted for use as a scaling
# coefficient of the normalized BRDFs. The normalized BRDFs have an integral of 1.
#
# Note that this source code requires a compiled C program to calculate the BRDF given the
# Zernike coefficients. The C source code is provided together with a Linux 64 
# compilation. Windows and Mac users have to compile as necessary.
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
# Calculate the nominal hemispherical-directional reflectance, rho(theta_v, lambda), for 
# each plaque at each view geometry and wavelength. The BRDF are interpolated for all 
# relevant directions in a 1x1 degree grid for the measurement wavelengths. This BRDF is 
# then integrated to rho and later interpolated to the wavelength grid used in this study.
# rho values are used as a scaling factor for the normalized BRDFs (integral = 1).
#

dir.create("sim/brdf")

coeffs      <- c(ptfe_white_351, ptfe_white_532, ptfe_white_633, ptfe_white_1064)	# original Zernike coefficients of the white PFTE
lbds        <- c(351, 532, 633, 1064)							# wavelength of BRDF measurements
ptfe_99_rho <- matrix(NA, ncol = 2, nrow = length(coeffs))				# Initial mstructure to receive white plaque results
colnames(ptfe_99_rho) <- c("v00", "v40")
rownames(ptfe_99_rho) <- lbds

for(lbd in 1:length(coeffs)) {
	temp <- zernike(coeffs[lbd], 0, 270, lbds[lbd])
	dim(temp) <- c(91, 361)
	intp   <- bicubic(x = rad(0:90), y = rad(0:360), z = temp, x0 = gxy[, 1], 
			y0 = gxy[, 2])$z
	intp   <- intp * fx
	dim(intp) <- c(length(sx), length(sy))
	intp   <- intp[-length(sx), ] + intp[-1, ]
	intp   <- intp[, -length(sy)] + intp[, -1]
	ptfe_99_rho[lbd, 1] <- sum(intp / 4)

	temp <- zernike(coeffs[lbd], 40, 270, lbds[lbd])
	dim(temp) <- c(91, 361)
	intp   <- bicubic(x = rad(0:90), y = rad(0:360), z = temp, x0 = gxy[, 1], 
			y0 = gxy[, 2])$z
	intp   <- intp * fx
	dim(intp) <- c(length(sx), length(sy))
	intp   <- intp[-length(sx), ] + intp[-1, ]
	intp   <- intp[, -length(sy)] + intp[, -1]
	ptfe_99_rho[lbd, 2] <- sum(intp / 4)
}

ptfe_99_rho <- cbind(approx(x = lbds, y = ptfe_99_rho[, 1], xout = lambda, rule = 2)$y,
                     approx(x = lbds, y = ptfe_99_rho[, 2], xout = lambda, rule = 2)$y)	# Interpolation for all wavelengths
colnames(ptfe_99_rho) <- c("v00", "v40")
rownames(ptfe_99_rho) <- lambda

coeffs      <- c(ptfe_grey_532)								# original Zernike coefficients of the grey PFTE
lbds        <- c(532)									# wavelength of BRDF measurements
ptfe_10_rho <- matrix(NA, ncol = 2, nrow = length(coeffs))				# Initial mstructure to receive grey plaque results
colnames(ptfe_10_rho) <- c("v00", "v40")
rownames(ptfe_10_rho) <- lbds

for(lbd in 1:length(coeffs)) {
	temp <- zernike(coeffs[lbd], 0, 270, lbds[lbd])
	dim(temp) <- c(91, 361)
	intp   <- bicubic(x = rad(0:90), y = rad(0:360), z = temp, x0 = gxy[, 1], 
			y0 = gxy[, 2])$z
	intp   <- intp * fx
	dim(intp) <- c(length(sx), length(sy))
	intp   <- intp[-length(sx), ] + intp[-1, ]
	intp   <- intp[, -length(sy)] + intp[, -1]
	ptfe_10_rho[lbd, 1] <- sum(intp / 4)

	temp <- zernike(coeffs[lbd], 40, 270, lbds[lbd])
	dim(temp) <- c(91, 361)
	intp   <- bicubic(x = rad(0:90), y = rad(0:360), z = temp, x0 = gxy[, 1], 
			y0 = gxy[, 2])$z
	intp   <- intp * fx
	dim(intp) <- c(length(sx), length(sy))
	intp   <- intp[-length(sx), ] + intp[-1, ]
	intp   <- intp[, -length(sy)] + intp[, -1]
	ptfe_10_rho[lbd, 2] <- sum(intp / 4)
}

#ptfe_10_rho <- cbind(approx(x = lbds, y = ptfe_10_rho[, 1], xout = lambda, rule = 2)$y,
#                     approx(x = lbds, y = ptfe_10_rho[, 2], xout = lambda, rule = 2)$y)
ptfe_10_rho <- ptfe_10_rho[rep(1, length(lambda)),]					# Since for now measurements are only available at 532 nm,  
colnames(ptfe_10_rho) <- c("v00", "v40")						# no interpolation is possible and data is just replicated.
rownames(ptfe_10_rho) <- lambda

#
# Compute the spectral BRDF for view geometries with the plaque level. This section 
# calculates the normalized BRDF for each wavelength and view geometry. The normalized 
# Zernike coefficients are used for the white plaque and the BRDF for the grey plaque is
# normalized on the fly.
# BRDF data is saved as three-dimensional arrays (zenith, azimuth, wavelength) and stored 
# in lists, with each geometry as a list element.
#

viewg <- rbind(c(0, 270), c(40, 270))

brdf_ptfe99 <- list()
for(g in 1:2) {
	brdf_ptfe99[[g]] <- list()
	brdf_ptfe99[[g]] <- rep(NA, 91*361)
	brdf_ptfe99[[g]] <- zernike_sp(ptfe_white_all, viewg[g, 1], viewg[g, 2]) 
	dim(brdf_ptfe99[[g]]) <- c(91, 361, length(lambda))
}
save(brdf_ptfe99, file = "sim/brdf/brdf_ptfe99.Rda")

brdf_ptfe10 <- list()
for(g in 1:2) {
	brdf_ptfe10[[g]] <- list()
	brdf_ptfe10[[g]] <- rep(NA, 91*361)
	brdf_ptfe10[[g]] <- zernike(ptfe_grey_532, viewg[g, 1], viewg[g, 2], 532) / ptfe_10_rho[1, g]
 	brdf_ptfe10[[g]] <- brdf_ptfe10[[g]][, rep(1, length(lambda))]
	dim(brdf_ptfe10[[g]]) <- c(91, 361, length(lambda))
}
save(brdf_ptfe10, file = "sim/brdf/brdf_ptfe10.Rda")

#
# Calculate the real BRDF under tilt, caused by changes in the view angle. The real 
# hemispherical-directional reflectance (rho) for each tilt/plaque/geometry is also 
# calculated for a correct scaling.
#

dir.create("sim/brdf/tilt", showWarnings = FALSE)
dir.create("sim/brdf/tilt/geometry_1", showWarnings = FALSE)
dir.create("sim/brdf/tilt/geometry_2", showWarnings = FALSE)

viewg <- rad(rbind(c(0, 270), c(40, 270)))
ptfes <- c("ptfe99", "ptfe10")

dimnms <- list(lambda, c('tv00', 'tv40'), c('znt03', 'znt06', 'znt09', 'znt12'), 
		c('azt000', 'azt045', 'azt090', 'azt135', 'azt180', 'azt225', 'azt270', 'azt315'))

ptfe_10_rho_tilt <- array(NA, dim = c(length(lambda), 2, 4, 8), dimnames = dimnms)	# rho under tilt for the grey plaque 
ptfe_99_rho_tilt <- array(NA, dim = c(length(lambda), 2, 4, 8), dimnames = dimnms)	# rho under tilt for the white plaque

coeffs <- c(ptfe_white_351, ptfe_white_532, ptfe_white_633, ptfe_white_1064)
cflbd  <- c(351, 532, 633, 1064)
ptfe_99_rho_tilt_l <- array(NA, dim = c(length(cflbd), 2, 4, 8))

for(pq in 1:length(ptfes)) {
	for(geom in 1:nrow(viewg)) {
		for(zentilt in c(3, 6, 9, 12)) {
			if(geom == 1) tts <- cbind(rep(zentilt, 8), (c(0, 45, 90, 135, 180, 225, 270, 315)+0) %% 360)
			if(geom == 2) tts <- cbind(rep(zentilt, 8), (c(0, 45, 90, 135, 180, 225, 270, 315)+90) %% 360)
			dir <- paste0(paste0("sim/brdf/tilt/geometry_", geom, "/"), formatC(zentilt, width = 2, flag = 0))
			dir.create(path = dir, showWarnings = F)
			viewB <- deg(viewg[rep(geom, 8), ])
			for(j in 1:nrow(viewB)) {
				dirm <- sphere.to.dir.cos(cbind(rad(viewB[j, 1]), rad(viewB[j, 2])))
				viewB[j,] <- deg(dir.cos.to.sphere(update.cosine(psi = rad(tts[j, 1]), 
						phi = rad(tts[j, 2]), cos.m = dirm)))
			}
			for(g in 1:nrow(viewB)) {
				if(pq == 1) {
					brdf <- zernike_sp(ptfe_white_all, viewB[g, 1], viewB[g, 2])
					dim(brdf) <- c(91, 361, length(lambda))

					for(lbd in 1:length(coeffs)) {
						brdf_l <- zernike(coeffs[lbd], viewB[g, 1], viewB[g, 2], cflbd[lbd])

						dim(brdf_l) <- c(91, 361)
						intp <- bicubic(x = rad(0:90), y = rad(0:360), z = brdf_l, x0 = gxy[, 1], y0 = gxy[, 2])$z
						intp <- intp * fx
						dim(intp) <- c(length(sx), length(sy))
						intp <- intp[-length(sx), ] + intp[-1, ]
						intp <- intp[, -length(sy)] + intp[, -1]
						idz  <- which(zets == zentilt)
						ptfe_99_rho_tilt_l[lbd, geom, idz, g] <- sum(intp / 4)
					}
				} else {
					brdf  <- zernike(ptfe_grey_532, viewB[g, 1], viewB[g, 2], 532)

					dim(brdf) <- c(91, 361)
					intp <- bicubic(x = rad(0:90), y = rad(0:360), z = brdf, x0 = gxy[, 1], y0 = gxy[, 2])$z
					intp <- intp * fx
					dim(intp) <- c(length(sx), length(sy))
					intp <- intp[-length(sx), ] + intp[-1, ]
					intp <- intp[, -length(sy)] + intp[, -1]
					idz  <- which(zets == zentilt)
					ptfe_10_rho_tilt[, geom, idz, g] <- rep(sum(intp / 4), length(lambda))
					brdf <- brdf / ptfe_10_rho_tilt[1, geom, idz, g]
					dim(brdf) <- c(91*361, 1)
					brdf <- brdf[, rep(1, length(lambda))]
					dim(brdf) <- c(91, 361, length(lambda))
				}
				if(geom == 1) azi = tts[g, 2]
				if(geom == 2) azi = c(0, 45, 90, 135, 180, 225, 270, 315)[g]
				fl <- paste0(dir, "/brdf_", ptfes[pq], "_", formatC(zentilt, width = 3, flag = 0), "_", formatC(azi, width = 3, flag = 0), ".Rda")
			 	save(brdf, file = fl)
			}
		}
	}
}

for(i in 1:4) {
	for(j in 1:8) {
		ptfe_99_rho_tilt[,,i,j] <- 
			c(approx(x = cflbd, y = ptfe_99_rho_tilt_l[, 1, i, j], xout = lambda, rule = 2)$y,
			  approx(x = cflbd, y = ptfe_99_rho_tilt_l[, 2, i, j], xout = lambda, rule = 2)$y)
	}
}

#
# Plots:
#

library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
	"\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
	"\\documentclass[18pt]{article}")

g = 1
l = 36
tikz("plot/brdf_pfte99_00_90.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
plot_brdf_latex(brdf_ptfe99[[g]][,,l] / ptfe_99_rho[l, g], viewg = c('0', '90'))
mtext('White PTFE (99 \\%)', side = 3, line = 1, cex = 2, col = 'black')
mtext('$\\theta_{\\text{v}} = 0^{\\circ}, \\phi_{\\text{v}} = 90^{\\circ}$', side = 2, line = 2, cex = 2, col = 'black')
dev.off()
tools::texi2pdf("plot/brdf_pfte99_00_90.tex")
file.rename("brdf_pfte99_00_90.pdf", "plot/brdf_pfte99_00_90.pdf")
file.remove("brdf_pfte99_00_90.aux", "brdf_pfte99_00_90.log")

tikz("plot/brdf_pfte10_00_90.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
plot_brdf_latex(brdf_ptfe10[[g]][,,l] / ptfe_99_rho[l, g], viewg = c('0', '90'))
mtext('Grey PTFE (10 \\%)', side = 3, line = 1, cex = 2, col = 'black')
#mtext('$\\theta_{\\text{v}} = 0^{\\circ}, \\phi_{\\text{v}} = 90^{\\circ}$', side = 2, line = 2, cex = 2, col = 'black')
dev.off()
tools::texi2pdf("plot/brdf_pfte10_00_90.tex")
file.rename("brdf_pfte10_00_90.pdf", "plot/brdf_pfte10_00_90.pdf")
file.remove("brdf_pfte10_00_90.aux", "brdf_pfte10_00_90.log")

g = 2

tikz("plot/brdf_pfte99_40_90.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
plot_brdf_latex(brdf_ptfe99[[g]][,,l] / ptfe_99_rho[l, g], viewg = c('40', '90'))
#mtext('White PTFE (99 \\%)', side = 3, line = 1, cex = 2, col = 'black')
mtext('$\\theta_{\\text{v}} = 40^{\\circ}, \\phi_{\\text{v}} = 90^{\\circ}$', side = 2, line = 2, cex = 2, col = 'black')
dev.off()
tools::texi2pdf("plot/brdf_pfte99_40_90.tex")
file.rename("brdf_pfte99_40_90.pdf", "plot/brdf_pfte99_40_90.pdf")
file.remove("brdf_pfte99_40_90.aux", "brdf_pfte99_40_90.log")

tikz("plot/brdf_pfte10_40_90.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
plot_brdf_latex(brdf_ptfe10[[g]][,,l] / ptfe_99_rho[l, g], viewg = c('40', '90'))
#mtext('Grey PTFE (10 \\%)', side = 3, line = 1, cex = 2, col = 'black')
#mtext('$\\theta_{\\text{v}} = 0^{\\circ}, \\phi_{\\text{v}} = 90^{\\circ}$', side = 2, line = 2, cex = 2, col = 'black')
dev.off()
tools::texi2pdf("plot/brdf_pfte10_40_90.tex")
file.rename("brdf_pfte10_40_90.pdf", "plot/brdf_pfte10_40_90.pdf")
file.remove("brdf_pfte10_40_90.aux", "brdf_pfte10_40_90.log")

cols <- heat.colors(255)
tikz("plot/brdf_legend.tex", width = 7, height = 2.3, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	par(mar = c(5, 2, 2, 2))
	plot(NA, xlim = c(0.20, 0.45), ylim = c(0, 1), xaxs = 'i', yaxs = 'i', xaxt = 'n', 
	yaxt = 'n', xlab = "", ylab = "")
	axis(1, at = c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45), labels = c('0.20', '0.25', '0.30', '0.35', '0.40', '$>0.45$'), cex.axis = 2)
	xrst <- t(as.raster(.map2color(seq(0.2, 0.45, length.out = length(cols)), cols)))
	rasterImage(xrst, xleft = 0.20, ybottom = 0, xright = 0.45, ytop = 1, interpolate = T)
	mtext('$f_{\\text{r}}(\\theta_{\\text{i}}, \\phi_{\\text{i}}, \\theta_{\\text{v}}, \\phi_{\\text{v}}, 550) / \\rho(\\theta_{\\text{v}}, 550)$', side = 1, line = 3, cex = 2.2)
	box()
dev.off()
tools::texi2pdf("plot/brdf_legend.tex")
file.rename("brdf_legend.pdf", "plot/brdf_legend.pdf")
file.remove("brdf_legend.aux", "brdf_legend.log")

