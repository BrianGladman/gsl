#ifndef COULOMB_H_
#define COULOMB_H_


/* Normalized hydrogenic bound states, radial dependence. */
double gsl_sf_hydrogenicR_1(double Z, double r);
double gsl_sf_hydrogenicR_2(int l, double Z, double r);
double gsl_sf_hydrogenicR(int n, int l, double Z, double r);


/* Coulomb wave function normalization constant.
   See Abramowitz+Stegun 14.1.8 14.1.9.
   coulomb_CL() calculates recursively in the same
   manner as coulomb_CL_list()
 */
double gsl_sf_coulomb_CL(double lam, double eta);
void   gsl_sf_coulomb_CL_list(double l_min, int count, double eta, double *cl);


/* Coulomb wave functions F,G for general lambda > -1
   conventions of Abramowitz and Stegun
 */
void gsl_sf_coulomb_wave_FGp(double x, double eta,
 		      	     double lam_min, double lam_max,
 		      	     double * fc, double * gc,
 		      	     double * fc_prime, double * gc_prime
	      	      	     );
void gsl_sf_coulomb_wave_FG(double x, double eta,
 		      	    double lam_min, double lam_max,
 		      	    double * fc, double * gc
	      	      	    );
void gsl_sf_coulomb_wave_F(double x, double eta,
	      	      	   double lam_min, double lam_max,
	              	   double * fc
	      	      	   );

/* Coulomb wave function divided by the argument,
   F(xi, eta)/xi. This is the function which reduces to
   spherical Bessel functions in the limit eta->0.
   The only issue is handling x near zero, which is trapped
   and handled using the explicit result for the behaviour
   of F(xi, eta) at the origin. If the x=~0 trap is sprung,
   the function returns 0, else it returns 1.
 */
int gsl_sf_coulomb_wave_sphF(double x, double eta,
		      	     double lam_min, double lam_max,
	              	     double * fc
	      	      	     );


/* value of overflow exponent for last Coulomb wave calculation
   If this integer, n, is nonzero, then the actual values of
   F,G,Fp,Gp are given by
            F(actual) = 10^(-n)  F(coulomb_wave)
           Fp(actual) = 10^(-n) Fp(coulomb_wave)
            G(actual) = 10^(n)   G(coulomb_wave)
           Gp(actual) = 10^(n)  Gp(coulomb_wave)
 */
int coul_wave_overflow_exp(void);


#endif /* !COULOMB_H_ */
