#
# Read OPAC output
#
# This function imports standard OPAC outputs into R. Is not meant to be generic.
#
# fl - OPAC out file path
# rh - desired relative humidity
#

read_opac <- function(fl, rh) {
	id <- which(c(0, 50, 70, 80, 90, 95, 98, 99) == rh)
	if(length(id) == 0) stop("Relative humidity must be one of: 0, 50, 70, 80, 90, 95, 98 or 99")
	data  <- readLines(fl)
	id.scat <- grep("112 scattering angles", data)
	scang <- as.numeric(t(matrix(unlist(strsplit(data[(id.scat+1):(id.scat+14)], "\\s+")), ncol = 9, byrow = T)[, -1]))
	id.pf <- grep("phase function", data)[-1]
	rh <- c("rh_00", "rh_50", "rh_70", "rh_80", "rh_90", "rh_95", "rh_98", "rh_99")

	opp <- array(as.numeric(unlist(strsplit(trimws(data[id.pf-1]), "\\s+"))), dim = c(6, 17, 8), 
               dimnames = list(c("lambda", "c",  "b",  "a",  "w0",  "tau"), NULL, rh))
	opp <- aperm(opp, perm = c(2, 1, 3))
	opp[,1,] <- opp[,1,] * 1000

	temp <- trimws(data[as.vector(apply(cbind(id.pf+1, id.pf+6), 1, function(x) x[1]:x[2]))])
	vsf   <- array(as.numeric(unlist(strsplit(temp, "\\s+"))), dim = c(112, 17, 8), dimnames = list(
                NULL, paste0("l", opp[,1,1]), rh))
	return(list(opp = opp[,,id], scang, vsf[,,id]))
}

