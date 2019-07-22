#
# sim_08_plaque_brdf.R
#
# Uses a C code extracted from the the SCATMECH library functions.
#

# Coefficients paths fro C program:

ptfe_white_all  <- "anc/ptfe_brdf_coef/PTFE.coeff.csv"
ptfe_white_351  <- "anc/ptfe_brdf_coef/PTFE.351.coeff.csv"
ptfe_white_532  <- "anc/ptfe_brdf_coef/PTFE.532.coeff.csv"
ptfe_white_633  <- "anc/ptfe_brdf_coef/PTFE.633.coeff.csv"
ptfe_white_1064 <- "anc/ptfe_brdf_coef/PTFE.1064.coeff.csv"
ptfe_grey_532   <- "anc/ptfe_brdf_coef/PTFE.grey10.532.coef.csv"

# Function:

zernike <- function(coefp, thetav, phiv, lambda) {
	system(paste("brdf_c/brdf.exe", coefp, rad(thetav), rad(phiv), lambda/1000, "> brdf_c/temp"))
	as.matrix(read.table("brdf_c/temp", header = F))
}

zernike_sp <- function(coefp, thetav, phiv) {
	system(paste("brdf_c/brdf_spec.exe", coefp, rad(thetav), rad(phiv), "> brdf_c/temp"))
	as.matrix(read.table("brdf_c/temp", header = F))
}

zernike_g <- function(coefp, thetai, phii, thetav, phiv, lambda) {
	as.numeric(system(paste("brdf_c/brdf_geom.exe", coefp, rad(thetai), rad(phii), rad(thetav), rad(phiv), lambda/1000), intern = T))
}


