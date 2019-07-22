This folder contains data reported in 

Thomas A. Germer, "Full four-dimensional and reciprocal Mueller matrix 
bidirectional reflectance distribution function of sintered 
polytetrafluoroethylene," Applied Optics, vol. xx, pp. xxxx-xxxx (2017).

The files are:

	README.txt - This file

	PTFE.1064.coeff.csv - Coefficients obtained using only regularly-spaced 1064 nm data
	PTFE.633.coeff.csv - Coefficients obtained using only regularly-spaced 633 nm data
	PTFE.532.coeff.csv - Coefficients obtained using only regularly-spaced 532 nm data
	PTFE.351.coeff.csv - Coefficients obtained using only regularly-spaced 351 nm data
	PTFE.coeff.csv - Coefficients obtained using all of the regularly-spaced data

	PTFE.data.csv - The regularly-spaced BRDF data
	PTFE.data.random.csv - The randomly-spaced BRDF data

	zernikeexpansion.cpp - C++ source code file for ZernikeExpansion_BRDF_Model, compatible with SCATMECH
	zernikeexpansion.h - C++ header file for ZernikeExpansion_BRDF_Model, compatible with SCATMECH

The columns for the coefficient files are:

	i,j,m,n,k,l,p,a

where 

	i,j - Mueller matrix index
	m,n - Zernike radial orders
	k,l - Zernike azimuthal orders
	p   - power of lambda (where lambda is in micrometers)
	a   - coefficient $a_{ij,mn}^{klp}$


The columns for the BRDF data files are

	lambda,thetai,phii,thetar,phir,fr11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44

where

	lambda - wavelength in micrometers
	thetai - incident polar angle in degrees
	phii   - incident azimuthal angle in degrees
	thetar - reflected polar angle in degrees
	phir   - reflected azimuthal angle in degrees
	fr11   - unpolarized BRDF in inverse steradians
	m12    - normalized Mueller matrix element 12
	m13    - normalized Mueller matrix element 13
	m14    - normalized Mueller matrix element 14
	m21    - normalized Mueller matrix element 21
	m22    - normalized Mueller matrix element 22
	m23    - normalized Mueller matrix element 23
	m24    - normalized Mueller matrix element 24
	m31    - normalized Mueller matrix element 31
	m32    - normalized Mueller matrix element 32
	m33    - normalized Mueller matrix element 33
	m34    - normalized Mueller matrix element 34
	m41    - normalized Mueller matrix element 41
	m42    - normalized Mueller matrix element 42
	m43    - normalized Mueller matrix element 43
	m44    - normalized Mueller matrix element 44

