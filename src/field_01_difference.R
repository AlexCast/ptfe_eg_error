#
# Process data from field inter-comparison experiments
#
# Irradiance data is provided as W/m2/nm.
#


#
# Import data from field inter-comparison:
#

asdzpocr_m <- list(
	hdasd_mat = cbind(read.csv('field/nadir_plaque_Eg_uconn_cruise.csv', row.names = 1), read.csv('field/side_plaque_Eg_uconn_cruise.csv', row.names = 1)),
	zpocr_mat = cbind(read.csv('field/nadir_ocr_1_Eg_uconn_cruise.csv', row.names = 1), read.csv('field/side_ocr_1_Eg_uconn_cruise.csv', row.names = 1))
)

asdhdocr_m <- list(
	hdasd_mat = cbind(read.csv('field/nadir_plaque_Eg_uconn_cruise.csv', row.names = 1), read.csv('field/side_plaque_Eg_uconn_cruise.csv', row.names = 1)),
	hdocr_mat = cbind(read.csv('field/nadir_ocr_2_Eg_uconn_cruise.csv', row.names = 1), read.csv('field/side_ocr_2_Eg_uconn_cruise.csv', row.names = 1))
)


#
# Calculate statistics for mean and standard deviation of absolute difference:
#

library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
	"\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
	"\\documentclass[18pt]{article}")

# Nadir:

mae <- apply(cbind(abs(asdzpocr_m[[1]][, 1:10] - asdzpocr_m[[2]][, 1:10]) / asdzpocr_m[[2]][, 1:10],
                   abs(asdhdocr_m[[1]][, 1:10] - asdhdocr_m[[2]][, 1:10]) / asdhdocr_m[[2]][, 1:10]), 1, mean, na.rm = T)
sae <- apply(cbind(abs(asdzpocr_m[[1]][, 1:10] - asdzpocr_m[[2]][, 1:10]) / asdzpocr_m[[2]][, 1:10],
                   abs(asdhdocr_m[[1]][, 1:10] - asdhdocr_m[[2]][, 1:10]) / asdhdocr_m[[2]][, 1:10]), 1, sd, na.rm = T)

round(mean(mae[51:351])*100, 2)
round(mean(sae[51:351])*100, 2)

tikz("plot/10_mean_relative_error_OCR_nadir.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	par(mar = c(5, 5, 3, 2))
	xlab <- '$\\lambda$ (nm)'
	ylab <- 'Mean Absolute Percentage Difference (\\%)'
	plot(350:800, mae * 100, xlim = c(355, 800), ylim = c(0, 20), type = "l", xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', lwd = 3)
	axis(2, at = seq(0, 20, 5), labels = TRUE, cex.axis = 1.9)
	axis(1, at = seq(400, 1000, 100), labels = TRUE, cex.axis = 1.9)
	mtext(xlab, side = 1, line = 3, cex = 2)
	mtext(ylab, side = 2, line = 3, cex = 1.7)
	polygon(x = c(355:795, 795:355), y = c(mae[-c(1:5, 447:451)]-sae[-c(1:5, 447:451)], rev(mae[-c(1:5, 447:451)]+sae[-c(1:5, 447:451)])) * 100, col = "grey90", border = "grey90")
	lines(350:800, mae * 100)
	lines(x = c(400, 700), y = c(mean(mae[51:351]), mean(mae[51:351])) * 100, col = "red", lwd = 3)
	abline(v = c(400, 700), col = "grey", lty = 2)
	text(x = 550, y = 15, labels = eval(substitute(expression('$\\overline{\\text{MAPD}}_{400:700}$'== x*'$\\,\\%$'~(N == y)), list(x = round(mean(mae[51:351])*100, 2), y = sum(!is.na(c(asdzpocr_m[[2]][51, 1:10], asdhdocr_m[[2]][51, 1:10])))))), col = "red", cex = 1.5)
	text(x = 550, y = 10, labels = '$\\text{MAPD} = \\frac{100}{\\text{N}}\\sum\\limits_{\\text{Sta.}}{\\frac{|\\text{ASD}-\\text{HyperOCR}|}{\\text{HyperOCR}}}$', cex = 1.5)
	text(x = 555, y = 20, labels = '$\\theta_{\\text{v}}=0^{\\circ}, \\phi_{\\text{v}} = 90^{\\circ}$', cex = 2)
dev.off()
tools::texi2pdf("plot/10_mean_relative_error_OCR_nadir.tex")
file.rename("10_mean_relative_error_OCR_nadir.pdf", "plot/10_mean_relative_error_OCR_nadir.pdf")
file.remove("10_mean_relative_error_OCR_nadir.aux", "10_mean_relative_error_OCR_nadir.log")


# Side:

mae <- apply(cbind(abs(asdzpocr_m[[1]][, 11:20] - asdzpocr_m[[2]][, 11:20]) / asdzpocr_m[[2]][, 11:20],
                   abs(asdhdocr_m[[1]][, 11:20] - asdhdocr_m[[2]][, 11:20]) / asdhdocr_m[[2]][, 11:20]), 1, mean, na.rm = T)
sae <- apply(cbind(abs(asdzpocr_m[[1]][, 11:20] - asdzpocr_m[[2]][, 11:20]) / asdzpocr_m[[2]][, 11:20],
                   abs(asdhdocr_m[[1]][, 11:20] - asdhdocr_m[[2]][, 11:20]) / asdhdocr_m[[2]][, 11:20]), 1, sd, na.rm = T)

round(mean(mae[51:351])*100, 2)
round(mean(sae[51:351])*100, 2)

tikz("plot/10_mean_relative_error_OCR_side.tex", width = 7, height = 7, standAlone = TRUE, pointsize = 18,
	packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
	"\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
	"\\usepackage{amssymb}", "\\usepackage{amsmath}"))
	par(mar = c(5, 5, 3, 2))
	xlab <- '$\\lambda$ (nm)'
	ylab <- 'Mean Absolute Percentage Difference (\\%)'
	plot(350:800, mae * 100, xlim = c(355, 800), ylim = c(0, 20), type = "l", xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', lwd = 3)
	axis(2, at = seq(0, 20, 5), labels = TRUE, cex.axis = 1.9)
	axis(1, at = seq(400, 1000, 100), labels = TRUE, cex.axis = 1.9)
	mtext(xlab, side = 1, line = 3, cex = 2)
	mtext(ylab, side = 2, line = 3, cex = 1.7)
	polygon(x = c(355:795, 795:355), y = c(mae[-c(1:5, 447:451)]-sae[-c(1:5, 447:451)], rev(mae[-c(1:5, 447:451)]+sae[-c(1:5, 447:451)])) * 100, col = "grey90", border = "grey90")
	lines(350:800, mae * 100)
	lines(x = c(400, 700), y = c(mean(mae[51:351]), mean(mae[51:351])) * 100, col = "red", lwd = 3)
	abline(v = c(400, 700), col = "grey", lty = 2)
	text(x = 550, y = 15, labels = eval(substitute(expression('$\\overline{\\text{MAPD}}_{400:700}$'== x*'$\\,\\%$'~(N == y)), list(x = round(mean(mae[51:351])*100, 2), y = sum(!is.na(c(asdzpocr_m[[2]][51, 11:20], asdhdocr_m[[2]][51, 11:20])))))), col = "red", cex = 1.5)
#	text(x = 550, y = 10, labels = '$\\text{MAPD} = \\frac{100}{\\text{N}}\\sum\\limits_{\\text{Sta.}}{\\frac{|\\text{ASD}-\\text{HyperOCR}|}{\\text{HyperOCR}}}$', cex = 1.5)
	text(x = 555, y = 20, labels = '$\\theta_{\\text{v}}=40^{\\circ}, \\phi_{\\text{v}} = 90^{\\circ}$', cex = 2)
dev.off()
tools::texi2pdf("plot/10_mean_relative_error_OCR_side.tex")
file.rename("10_mean_relative_error_OCR_side.pdf", "plot/10_mean_relative_error_OCR_side.pdf")
file.remove("10_mean_relative_error_OCR_side.aux", "10_mean_relative_error_OCR_side.log")


png("plot/10_mean_relative_error_OCR_side.png", height = 1200, width = 1200, res = 230, bg = "transparent")
par(mar = c(5, 5, 3, 2))
plot(350:800, mae*100, xlim = c(355, 800), ylim = c(0, 20), type = "l", xlab = expression(lambda~(nm)), ylab = "Mean Absolute Percentage Error")
polygon(x = c(355:795, 795:355), y = c(mae[-c(1:5, 447:451)]-sae[-c(1:5, 447:451)], rev(mae[-c(1:5, 447:451)]+sae[-c(1:5, 447:451)])) * 100, col = "grey90", border = "grey90")
lines(350:800, mae * 100)
lines(x = c(400, 700), y = c(mean(mae[51:351]), mean(mae[51:351])) * 100, col = "red")
abline(v = c(400, 700), col = "grey", lty = 2)
mtext(text = eval(substitute(expression(bar(MAPE)[400:700] == x*"%"~(N == y)), list(x = round(mean(mae[51:351])*100, 2), y = sum(!is.na(c(asdzpocr_m[[2]][51, 11:20], asdhdocr_m[[2]][51, 11:20])))))), side = 1, line = -10.5, col = "red")
#mtext(text = expression(MAPE == over(sum(over(abs(ASD-OCR), OCR)),N) * 100), side = 3, line = -3, cex = 0.8)
legend("top", expression(theta[v]==40*"º,"~phi[v]==90*"º"), bty = "n", cex = 0.9)
dev.off()

