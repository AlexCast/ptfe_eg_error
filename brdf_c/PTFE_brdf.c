// Retrieved and adapted from the SCATMECH library
// PTFE_brdf.c
// gcc PTFE_brdf.c -o brdf.exe -lm


#include <stdio.h>		// Functions: fprintf, printf,
#include <stdlib.h>		// Functions: calloc
#include <math.h>

double factorial_list[171] =
    {   1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880.,3.6288e6,3.99168e7,4.790016e8,6.2270208e9,
        8.71782912e10,1.307674368e12,2.0922789888e13,3.55687428096e14,6.402373705728e15,
        1.21645100408832e17,
        2.43290200817664e18,5.109094217170944e19,1.1240007277776077e21,2.585201673888498e22,
        6.204484017332394e23,1.5511210043330986e25,4.0329146112660565e26,1.0888869450418352e28,
        3.0488834461171387e29,8.841761993739702e30,2.6525285981219107e32,8.222838654177922e33,
        2.631308369336935e35,8.683317618811886e36,2.9523279903960416e38,1.0333147966386145e40,
        3.7199332678990125e41,1.3763753091226346e43,5.230226174666011e44,2.0397882081197444e46,
        8.159152832478977e47,3.345252661316381e49,1.40500611775288e51,6.041526306337383e52,
        2.658271574788449e54,1.1962222086548019e56,5.502622159812089e57,2.5862324151116818e59,
        1.2413915592536073e61,6.082818640342675e62,3.0414093201713376e64,1.5511187532873822e66,
        8.065817517094388e67,4.2748832840600255e69,2.308436973392414e71,1.2696403353658276e73,
        7.109985878048635e74,4.0526919504877214e76,2.3505613312828785e78,1.3868311854568984e80,
        8.32098711274139e81,5.075802138772248e83,3.146997326038794e85,1.98260831540444e87,
        1.2688693218588417e89,8.247650592082472e90,5.443449390774431e92,3.647111091818868e94,
        2.4800355424368305e96,1.711224524281413e98,1.1978571669969892e100,8.504785885678623e101,
        6.1234458376886085e103,4.4701154615126844e105,3.307885441519386e107,2.48091408113954e109,
        1.8854947016660504e111,1.4518309202828587e113,1.1324281178206297e115,8.946182130782976e116,
        7.156945704626381e118,5.797126020747368e120,4.753643337012842e122,3.945523969720659e124,
        3.314240134565353e126,2.81710411438055e128,2.4227095383672734e130,2.107757298379528e132,
        1.8548264225739844e134,1.650795516090846e136,1.4857159644817615e138,1.352001527678403e140,
        1.2438414054641308e142,1.1567725070816416e144,1.087366156656743e146,1.032997848823906e148,
        9.916779348709496e149,9.619275968248212e151,9.426890448883248e153,9.332621544394415e155,
        9.332621544394415e157,9.42594775983836e159,9.614466715035127e161,9.90290071648618e163,
        1.0299016745145628e166,1.081396758240291e168,1.1462805637347084e170,1.226520203196138e172,
        1.324641819451829e174,1.4438595832024937e176,1.588245541522743e178,1.7629525510902446e180,
        1.974506857221074e182,2.2311927486598138e184,2.5435597334721877e186,2.925093693493016e188,
        3.393108684451898e190,3.969937160808721e192,4.684525849754291e194,5.574585761207606e196,
        6.689502913449127e198,8.094298525273444e200,9.875044200833601e202,1.214630436702533e205,
        1.506141741511141e207,1.882677176888926e209,2.372173242880047e211,3.0126600184576594e213,
        3.856204823625804e215,4.974504222477287e217,6.466855489220474e219,8.47158069087882e221,
        1.1182486511960043e224,1.4872707060906857e226,1.9929427461615188e228,2.6904727073180504e230,
        3.659042881952549e232,5.012888748274992e234,6.917786472619489e236,9.615723196941089e238,
        1.3462012475717526e241,1.898143759076171e243,2.695364137888163e245,3.854370717180073e247,
        5.5502938327393044e249,8.047926057471992e251,1.1749972043909107e254,1.727245890454639e256,
        2.5563239178728654e258,3.80892263763057e260,5.713383956445855e262,8.62720977423324e264,
        1.3113358856834524e267,2.0063439050956823e269,3.0897696138473508e271,4.789142901463394e273,
        7.471062926282894e275,1.1729568794264145e278,1.853271869493735e280,2.9467022724950384e282,
        4.7147236359920616e284,7.590705053947219e286,1.2296942187394494e289,2.0044015765453026e291,
        3.287218585534296e293,5.423910666131589e295,9.003691705778438e297,1.503616514864999e300,
        2.5260757449731984e302,4.269068009004705e304,7.257415615307999e306
    };

double mpow(int m) {
	int result;
	if (m<0) m=-m;
	if (m & 0x1) result = -1;
	else result = 1;
	return result;
}

double Fact(int j) {
	if (j <= 170 && j>=0) return factorial_list[j];
	else exit(-1);
}

double Radial_Zernike_Polynomial(int n, int l, double rho) {
	int m = abs(l);
	if ((n-m)%2 == 1) return 0.;
	double sum = 0;
	for (int s = 0; s <= (n-m) / 2; ++s) {
		sum += mpow(s)*Fact(n-s)/Fact(s)/Fact((n+m)/2-s)/Fact((n-m)/2-s)*pow(rho,(double)(n-2*s));		
	}
	return sum;
}


double ZernikeExpansion_BRDF_Model(     double  thetai, 
					double  thetas, 
					double  phii, 
					double  phis, 
					double  lambda, 
					int     lines,
					int*    iv, 
					int*    jv, 
					int*    nv, 
					int*    mv, 
					int*    kv, 
					int*    lv, 
					int*    pv, 
					double* av) {
	double rotation = phii;
	double brdf = 0;

	double rhoi = sqrt(2.)*sin(thetai/2.);
	double rhor = sqrt(2.)*sin(thetas/2.);
	phii = phii - rotation;
	double phir = phis - rotation;

	double lambdapower[5];
	lambdapower[0] = 1.;
	lambdapower[1] = lambda;
	lambdapower[2] = lambdapower[1] * lambda;
	lambdapower[3] = lambdapower[2] * lambda;
	lambdapower[4] = lambdapower[3] * lambda;

	for(int q = 0; q < lines; q++) {
		if(iv[q] == 1 && jv[q] == 1) {
			int i = iv[q];
			int j = jv[q];
			int n = nv[q];
			int m = mv[q];
			int k = kv[q];
			int l = lv[q];
			int p = pv[q];
			double a = av[q];

			double power = lambdapower[p];

			double prefact = 0.5/M_PI*sqrt((n+1.)*(m+1.)/((n==0||(n==m&&l==0))?4.:((n==m||l==0)?2.:1.)))*power;
			brdf += a * prefact * (Radial_Zernike_Polynomial(n,l,rhoi)*Radial_Zernike_Polynomial(m,l,rhor) +
				 Radial_Zernike_Polynomial(m,l,rhoi)*Radial_Zernike_Polynomial(n,l,rhor)) * cos(l*(phir-phii));
		}
	}

	return brdf;
}

int main(int argc, char *argv[]) {

	FILE*   fpi;						// Pointer for input file
	int     lines = 0;					// Number of lines in coefficient file
	int*    i, *j, *n, *m, *k, *l, *p;			// Pointers to integers
	double* a;						// Pointers to doubles
	int c;

	// Read coefficient file:

	if (argc != 5) {
		printf ("Usage: %s inputfile thetas phis lambda\n", argv[0]);
		exit (-1);
	}
 
	fpi = fopen(argv[1], "r");
	double thetas = atof(argv[2]);
	double phis = atof(argv[3]);
	double lambda = atof(argv[4]);


//printf("%.4f %.4f %.4f \n", thetas, phis, lambda);

	if(fpi == NULL) {
		printf("Failed to open input file %s\n", argv[1]);
		exit(-1);
	}

	char ch;
	while(!feof(fpi)) {
		ch = fgetc(fpi);
		if(ch == '\n') {
			lines++;
		}
	}
	fseek(fpi, 0, SEEK_SET);

	i = (int*) calloc(lines, sizeof(int));
	j = (int*) calloc(lines, sizeof(int));
	n = (int*) calloc(lines, sizeof(int));
	m = (int*) calloc(lines, sizeof(int));
	k = (int*) calloc(lines, sizeof(int));
	l = (int*) calloc(lines, sizeof(int));
	p = (int*) calloc(lines, sizeof(int));
	a = (double*) calloc(lines, sizeof(double));

	for(c = 0; c < lines; c++) {
		fscanf(fpi, "%d, %d, %d, %d, %d, %d, %d, %lf", &i[c], &j[c], 
			&n[c], &m[c], &k[c], &l[c], &p[c], &a[c]);
//		printf("%d, %d, %d, %d, %d, %d, %d, %f \n", i[c], j[c], 
//			n[c], m[c], k[c], l[c], p[c], a[c]);
//		printf("\n\n");
	}

	// Create grid:

	double thetai_v[32851] = {0};
	double phii_v[32851] = {0};
	double valt = 0;
	double valp = 0;

	for(c = 0; c < 32851; c++) {
		thetai_v[c] = valt * M_PI / 180.;
		phii_v[c] = valp * M_PI / 180.;
		valt += 1.;
		if(valt > 90.) { 
			valt = 0.;
			valp += 1.;
		}
	}

	double brdf[32851] = {0};
	for(c = 0; c < 32851; c++) {
		brdf[c] = ZernikeExpansion_BRDF_Model(
						thetai_v[c], 
						thetas, 
						phii_v[c], 
						phis,
						lambda, 
						lines, i, j, n, m, k, l, p, a);	
		printf("%.4f\n", brdf[c]);
	}

	fclose(fpi);

	return 0;
}
