/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_COUPLING_H_
#define GSL_SF_COUPLING_H_

#include <gsl/gsl_sf_result.h>


/* 3j Symbols:  / ja jb jc \
 *              \ ma mb mc /
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_coupling_3j_impl(int two_ja, int two_jb, int two_jc,
                            int two_ma, int two_mb, int two_mc,
			    gsl_sf_result * result
			    );
int gsl_sf_coupling_3j_e(int two_ja, int two_jb, int two_jc,
                         int two_ma, int two_mb, int two_mc,
		         gsl_sf_result * result
			 );


/* 6j Symbols:  / ja jb jc \
 *              \ jd je jf /
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_coupling_6j_impl(int two_ja, int two_jb, int two_jc,
                            int two_jd, int two_je, int two_jf,
			    gsl_sf_result * result
			    );
int gsl_sf_coupling_6j_e(int two_ja, int two_jb, int two_jc,
                         int two_jd, int two_je, int two_jf,
			 gsl_sf_result * result
			 );


/* 9j Symbols:  / ja jb jc \
 *              | jd je jf |
 *              \ jg jh ji /
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_coupling_9j_impl(int two_ja, int two_jb, int two_jc,
                            int two_jd, int two_je, int two_jf,
			    int two_jg, int two_jh, int two_ji,
			    gsl_sf_result * result
			    );
int gsl_sf_coupling_9j_e(int two_ja, int two_jb, int two_jc,
                         int two_jd, int two_je, int two_jf,
			 int two_jg, int two_jh, int two_ji,
			 gsl_sf_result * result
			 );


#endif  /* !GSL_SF_COUPLING_H_ */
