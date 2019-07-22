#
# Example code.
#
# The examples below are provided to help the extraction of statistics over different data
# ranges than those used in the study.
#

# IMPORTANT - It is necessary first to restore the R workspace:
load('restore.RData')

# Example 1
#
# Extract error statistics for the BRDF effect alone:
#
# First the arrays to receive the results (error_brdf_avem and error_brdf_sdm) are created
# with names in each dimension to facilitate interpretation. Than data is extracted for 
# two sets of conditions: (1) RSky: Clear skies with Sun zenith angle between 20 and 60 
# degrees (corresponding data indexes of 4 to 8); and (2) OC: Overcast condition 
# (corresponding data index of 1). Statistics are retrieved in the visible range only 
# (wavelengths from 400 to 700 nm, with corresponding data indexes of 6 to 36). Statistics 
# are calculated both for retrieval with Eq. 5a and with Eq. 5b, and looped through both 
# plaques and geometries. Statistics are the mean absolute percentage error and the 
# standard deviation of the absolute percentage error. The statistics are calculated over 
# all wavelengths and sun zeniths included and for all 24 sky models. 
#

dimnms <- list(c('RSky', 'OC'), c('tv00', 'tv40'), ptfes, c('Eq5A', 'Eq5B'))

error_brdf_avem <- array(NA, dim = c(2, 2, 2, 2), dimnames = dimnms)
error_brdf_sdm  <- array(NA, dim = c(2, 2, 2, 2), dimnames = dimnms)

for(geom in 1:2) {
	for(p in 1:2) {
		temp <- abs(1 - res_brdf_tot[6:36, 4:8, geom, p, ])			# Convert from ratio to absolute fraction difference
		error_brdf_avem[1, geom, p, 1] <- mean(temp)
		error_brdf_sdm[1, geom, p, 1]  <- sd(temp)

		temp <- abs(1 - res_brdf_tot[6:36, 1, geom, p, 1])			# Just for the overcast condition. Note that only one sky model is used because under OC there is no difference.
		error_brdf_avem[2, geom, p, 1] <- mean(temp)
		error_brdf_sdm[2, geom, p, 1]  <- sd(temp)

		temp <- abs(1 - res_brdf_tot_pi[6:36, 4:8, geom, p, ])			# The same as the two above, but for equation 5b
		error_brdf_avem[1, geom, p, 2] <- mean(temp)
		error_brdf_sdm[1, geom, p, 2]  <- sd(temp)

		temp <- abs(1 - res_brdf_tot_pi[6:36, 1, geom, p, 1])
		error_brdf_avem[2, geom, p, 2] <- mean(temp)
		error_brdf_sdm[2, geom, p, 2]  <- sd(temp)
	}
}

round(error_brdf_avem * 100, 2)
round(error_brdf_sdm * 100, 2)

# Example 2
#
# Extract error statistics for the tilt effect alone:
#
# This can be of interest since is also applicable to irradiance sensors. Statistic are 
# calculated for the interval range of tilt from 0 to each tilt step (3 to 12, in 3 degree
# step). So first is necessary to include the error when no tilt is present - since for
# the effect of tilt alone BRDF and shadow are not considered, when under tilt of 0, there
# is no error. Note that in addition to the previous example, statistics are also 
# calculated over all azimuths. Since BRDF and shadow are not included, errors matrix is
# simpler.
#

dimnms <- list(c('znt03', 'znt06', 'znt09', 'znt12'), c('RSky', 'OC'))

error_tilt_avem <- array(NA, dim = c(length(zets), 2), dimnames = dimnms)
error_tilt_sdm  <- array(NA, dim = c(length(zets), 2), dimnames = dimnms)

temp <- array(1, dim = c(1, length(azts), length(lambda), 10, length(type)))		# Make array for the condition of no tilt. Under no Tilt, no BRDF and no Shadow, Eg estimated == Eg real
library(abind)
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

# Example 3
#
# Visualize sky radiance distribution.
#

library(sp)
library(raster)

sim_par <- read_simpar('sim/skyrad/continental_clean_95_01_sza_40_simpar.txt')
skym1   <- read.table('sim/skyrad/continental_clean_95_01_sza_40_sky.txt', header = T)
skym1   <- as.matrix(skym1[, -c(1:2), drop = F])
dim(skym1) <- c(length(sim_par$zenith), length(sim_par$razimuth)+180, length(sim_par$lambda))

plot_skyr(skym1[ , , 21]) 								# Raster legend fixed to 550 nm, corresponding to index 21.


# Example 4
#
# Visualize BRDF.
#

library(sp)
library(raster)

geom <- 1 
plot_brdf(brdf_ptfe99[[geom]][ , , 21], viewg = c(0, 90))
plot_brdf(brdf_ptfe10[[geom]][ , , 21], viewg = c(0, 90))

geom <- 2 
plot_brdf(brdf_ptfe99[[geom]][ , , 21], viewg = c(40, 90))
plot_brdf(brdf_ptfe10[[geom]][ , , 21], viewg = c(40, 90))

