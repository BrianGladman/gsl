/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_COUPLING_H_
#define GSL_SF_COUPLING_H_


/* Wigner coefficients: (ja jb ma mb | ja jb jc mc) 
 */
int gsl_sf_coupling_3j_impl(int two_ja, int two_jb, int two_jc,
                            int two_ma, int two_mb, int two_mc,
			    double * result
			    );
int gsl_sf_coupling_3j_e(int two_ja, int two_jb, int two_jc,
                         int two_ma, int two_mb, int two_mc,
		         double * result
			 );
double gsl_sf_coupling_3j(int two_ja, int two_jb, int two_jc,
                          int two_ma, int two_mb, int two_mc
			  );


/* 6j Coefficients: (ja jb jc | jd je jf)
 */
int gsl_sf_coupling_6j_impl(int two_ja, int two_jb, int two_jc,
                            int two_jd, int two_je, int two_jf,
			    double * result
			    );
int gsl_sf_coupling_6j_e(int two_ja, int two_jb, int two_jc,
                         int two_jd, int two_je, int two_jf,
			 double * result
			 );
double gsl_sf_coupling_6j(int two_ja, int two_jb, int two_jc,
                          int two_jd, int two_je, int two_jf
			  );


/* 9j Coefficients: ( )
 */
int gsl_sf_coupling_9j_impl(int two_ja, int two_jb, int two_jc,
                            int two_jd, int two_je, int two_jf,
			    int two_jg, int two_jh, int two_ji,
			    double * result
			    );
int gsl_sf_coupling_9j_e(int two_ja, int two_jb, int two_jc,
                         int two_jd, int two_je, int two_jf,
			 int two_jg, int two_jh, int two_ji,
			 double * result
			 );
double gsl_sf_coupling_9j(int two_ja, int two_jb, int two_jc,
                          int two_jd, int two_je, int two_jf,
			  int two_jg, int two_jh, int two_ji
			  );


#endif  /* !GSL_SF_COUPLING_H_ */
