/* Author: G. Jungman
 * RCS:    $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_psi.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* Chebyshev fit for f(y) = Re(Psi(1+Iy)) + M_EULER - y^2/(1+y^2) - y^2/(2(4+y^2))
 * 1 < y < 10
 *   ==>
 * y(x) = (9x + 11)/2,  -1 < x < 1
 * x(y) = (2y - 11)/9
 *
 * g(x) := f(y(x))
 */
static double r1py_data[] = {
   1.59888328244976954803168395603,
   0.67905625353213463845115658455,
  -0.068485802980122530009506482524,
  -0.0057881841830958667920088311823,
   0.0085112581671086159804198556475,
  -0.0040426561346996934343345564091,
   0.00135232840615940260177846295622,
  -0.000311646563930660566674525382102,
   0.0000185075637852491354372191392113,
   0.0000283487054275298502964921455903,
  -0.0000194875360145745355675419596539,
   8.0709788710834469408621587335e-6,
  -2.29835643213405180370603465611e-6,
   3.05066295996047498438559626587e-7,
   1.30422386324183646107742848462e-7,
  -1.23086571810489505894646902083e-7,
   5.7710855710682427240667414345e-8,
  -1.82755593424509639660926363536e-8,
   3.10204713006265894207595189301e-9,
   6.8989327480593812470039430640e-10,
  -8.7182290258923059852334818997e-10,
   4.4069147710243611798213548777e-10,
  -1.47273110991985359634672002769e-10,
   2.75896825232626447488258442482e-11,
   4.1871826756975856411554363568e-12,
  -6.5673460487260087541400767340e-12,
   3.4487900886723214020103638000e-12,
  -1.18072514174486906079737940779e-12,
   2.37983143439695892587093155740e-13,
   2.16636304108188318242594658208e-15
};
static gsl_sf_cheb_series r1py_cs = {
  r1py_data,
  29,
  -1,1,
  (double *)0,
  (double *)0,
  18
};


/* Chebyshev fits from SLATEC code for psi(x)

 Series for PSI        on the interval  0.	   to  1.00000D+00
				       with weighted error   2.03E-17
					log weighted error  16.69
			      significant figures required  16.39
				   decimal places required  17.37

 Series for APSI       on the interval  0.	   to  2.50000D-01
				       with weighted error   5.54E-17
					log weighted error  16.26
			      significant figures required  14.42
				   decimal places required  16.86

*/

static double psics_data[23] = {
  -.038057080835217922,
   .491415393029387130, 
  -.056815747821244730,
   .008357821225914313,
  -.001333232857994342,
   .000220313287069308,
  -.000037040238178456,
   .000006283793654854,
  -.000001071263908506,
   .000000183128394654,
  -.000000031353509361,
   .000000005372808776,
  -.000000000921168141,
   .000000000157981265,
  -.000000000027098646,
   .000000000004648722,
  -.000000000000797527,
   .000000000000136827,
  -.000000000000023475,
   .000000000000004027,
  -.000000000000000691,
   .000000000000000118,
  -.000000000000000020
};
static gsl_sf_cheb_series psi_cs = {
  psics_data,
  22,
  -1, 1,
  (double *)0,
  (double *)0,
  17
};

static double apsics_data[16] = {    
  -.0204749044678185,
  -.0101801271534859,
   .0000559718725387,
  -.0000012917176570,
   .0000000572858606,
  -.0000000038213539,
   .0000000003397434,
  -.0000000000374838,
   .0000000000048990,
  -.0000000000007344,
   .0000000000001233,
  -.0000000000000228,
   .0000000000000045,
  -.0000000000000009,
   .0000000000000002,
  -.0000000000000000 
};    
static gsl_sf_cheb_series apsi_cs = {
  apsics_data,
  15,
  -1, 1,
  (double *)0,
  (double *)0,
  9
};

#define PSI_TABLE_NMAX 100
static double psi_table[PSI_TABLE_NMAX+1] = {
  0.0,  /* Infinity */              /* psi(0) */
 -M_EULER,                          /* psi(1) */
  0.42278433509846713939348790992,  /* ...    */
  0.92278433509846713939348790992,
  1.25611766843180047272682124325,
  1.50611766843180047272682124325,
  1.70611766843180047272682124325,
  1.87278433509846713939348790992,
  2.01564147795560999653634505277,
  2.14064147795560999653634505277,
  2.25175258906672110764745616389,
  2.35175258906672110764745616389,
  2.44266167997581201673836525479,
  2.52599501330914535007169858813,
  2.60291809023222227314862166505,
  2.67434666166079370172005023648,
  2.74101332832746036838671690315,
  2.80351332832746036838671690315,
  2.86233685773922507426906984432,
  2.91789241329478062982462539988,
  2.97052399224214905087725697883,
  3.02052399224214905087725697883,
  3.06814303986119666992487602645,
  3.11359758531574212447033057190,
  3.15707584618530734186163491973,
  3.1987425128519740085283015864,
  3.2387425128519740085283015864,
  3.2772040513135124700667631249,
  3.3142410883505495071038001619,
  3.3499553740648352213895144476,
  3.3844381326855248765619282407,
  3.4177714660188582098952615740,
  3.4500295305349872421533260902,
  3.4812795305349872421533260902,
  3.5115825608380175451836291205,
  3.5409943255438998981248055911,
  3.5695657541153284695533770196,
  3.5973435318931062473311547974,
  3.6243705589201332743581818244,
  3.6506863483938174848844976139,
  3.6763273740348431259101386396,
  3.7013273740348431259101386396,
  3.7257176179372821503003825420,
  3.7495271417468059598241920658,
  3.7727829557002943319172153216,
  3.7955102284275670591899425943,
  3.8177324506497892814121648166,
  3.8394715810845718901078169905,
  3.8607481768292527411716467777,
  3.8815815101625860745049801110,
  3.9019896734278921969539597029,
  3.9219896734278921969539597029,
  3.9415975165651470989147440166,
  3.9608282857959163296839747858,
  3.9796962103242182164764276160,
  3.9982147288427367349949461345,
  4.0163965470245549168131279527,
  4.0342536898816977739559850956,
  4.0517975495308205809735289552,
  4.0690389288411654085597358518,
  4.0859880813835382899156680552,
  4.1026547480502049565823347218,
  4.1190481906731557762544658694,
  4.1351772229312202923834981274,
  4.1510502388042361653993711433,
  4.1666752388042361653993711433,
  4.1820598541888515500147557587,
  4.1972113693403667015299072739,
  4.2121367424746950597388624977,
  4.2268426248276362362094507330,
  4.2413353784508246420065521823,
  4.2556210927365389277208378966,
  4.2697055997787924488475984600,
  4.2835944886676813377364873489,
  4.2972931188046676391063503626,
  4.3108066323181811526198638761,
  4.3241399656515144859531972094,
  4.3372978603883565912163551041,
  4.3502848733753695782293421171,
  4.3631053861958823987421626300,
  4.3757636140439836645649474401,
  4.3882636140439836645649474401,
  4.4006092930563293435772931191,
  4.4128044150075488557724150703,
  4.4248526077786331931218126607,
  4.4367573696833950978837174226,
  4.4485220755657480390601880108,
  4.4601499825424922251066996387,
  4.4716442354160554434975042364,
  4.4830078717796918071338678728,
  4.4942438268358715824147667492,
  4.5053549379469826935258778603,
  4.5163439489359936825368668713,
  4.5272135141533849868846929582,
  4.5379662023254279976373811303,
  4.5486045001977684231692960239,
  4.5591308159872421073798223397,
  4.5695474826539087740464890064,
  4.5798567610044242379640147796,
  4.5900608426370772991885045755,
  4.6001618527380874001986055856
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_psi_int_impl(const int n, double * result)
{
  if(n <= 0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(n <= PSI_TABLE_NMAX) {
    *result = psi_table[n];
    return GSL_SUCCESS;
  }
  else if(n < 500) {
    double sum = -M_EULER;
    int k;
    for(k=1; k<n; k++) {
      sum += 1.0/k;
    }
    *result = sum;
    return GSL_SUCCESS;
  }
  else {
    /* Abramowitz+Stegun 6.3.18 */
    const double c2 = -1.0/12.0;
    const double c3 =  1.0/120.0;
    const double c4 = -1.0/252.0;
    const double ni2 = (1.0/n)*(1.0/n);
    const double ser = ni2 * (c2 + ni2 * (c3 + c4 * ni2));
    *result   = log(n) - 0.5/n + ser;
    return GSL_SUCCESS;
  }
}


int gsl_sf_psi_impl(const double x, double * result)
{
  double y = fabs(x);
  double xbig  = 1.0  / GSL_SQRT_MACH_EPS;    /* XBIG  = 1.0/SQRT(R1MACH(3)) */
  double dxrel = 10.0 * GSL_SQRT_MACH_EPS;    /* DXREL = SQRT (R1MACH(4))    */
  
  if(y >= 2.0) {
    double aux = 0.;
    if(y < xbig) aux = gsl_sf_cheb_eval(&apsi_cs, 8.0/(y*y)-1.0);
    if(x < 0.0) {
      /* *result = log(y) - 0.5/x + aux - M_PI * cot(M_PI*x); */
      *result = log(y) - 0.5/x + aux - M_PI * cos(M_PI*x)/sin(M_PI*x);
    }
    else {
      *result = log(y) - 0.5/x + aux;
    }
    return GSL_SUCCESS;
  }
  else { /* y < 2.0 */
    if(x == 0.0) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EDOM;
    }
    else {
      double ans;
      int n = x;
      if(x < 0.0) --n;
      y = x - n;
      --n;
      ans = gsl_sf_cheb_eval(&psi_cs, 2.0*y-1.0);
      if(n == 0) {
	*result = ans;
	return GSL_SUCCESS;
      }

      n = -n;

      if(x < 0.0 && x+n-2 == 0.) {
      	/* x is a negative integer */
	*result = 0.; /* FIXME: should be Inf */
	return GSL_EDOM;
      }
      else {
	int i;
	for(i=0; i<n; i++) {
          ans -= 1.0/(x + i);
      	}
	*result = ans;
	if(x < -0.5 && fabs((x-(int)(x-0.5))/x) < dxrel) {
      	  /* loss of precision: x near a negative integer */
	  return GSL_ELOSS;
      	}
	else {
	  return GSL_SUCCESS;
	}
      }
    }
  }
}


int
gsl_sf_psi_1piy_impl(const double y, double * result)
{
  double ay = fabs(y);

  if(ay > 1000.0) {
    /* [Abramowitz+Stegun, 6.3.19] */
    double yi2 = 1.0/(ay*ay);
    double ln  = log(y);
    double sum = yi2 * (1.0/12.0 + 1.0/120.0 * yi2 + 1.0/252.0 * yi2*yi2);
    *result = ln + sum;
    return GSL_SUCCESS;
  }
  else if(ay > 10.0) {
    /* [Abramowitz+Stegun, 6.3.19] */
    double yi2 = 1.0/(ay*ay);
    double ln  = log(y);
    double sum = yi2 * (1.0/12.0 +
                   yi2 * (1.0/120.0 +
		     yi2 * (1.0/252.0 +
                       yi2 * (1.0/240.0 +
		         yi2 * (1.0/132.0 + 691.0/32760.0 * yi2)))));
    *result = ln + sum;
    return GSL_SUCCESS;
  }
  else if(ay > 1.0){
    double y2 = ay*ay;
    double x = (2.0*ay - 11.0)/9.0;
    double c = gsl_sf_cheb_eval(&r1py_cs, x);
    *result = c - M_EULER + y2*(1.0/(1.0+y2) + 0.5/(4.0+y2));
    return GSL_SUCCESS;
  }
  else {
    /* [Abramowitz+Stegun, 6.3.17]
     *
     * -M_EULER + y^2 Sum[1/n 1/(n^2 + y^2), {n,1,M}]
     *   +     Sum[1/n^3, {n,M+1,Infinity}]
     *   - y^2 Sum[1/n^5, {n,M+1,Infinity}]
     *   + y^4 Sum[1/n^7, {n,M+1,Infinity}]
     *   - y^6 Sum[1/n^9, {n,M+1,Infinity}]
     *   + O(y^8)
     *
     * We take M=50 for at least 15 digit precision.
     */
    const int M = 50;
    const double y2 = y*y;
    const double c0 = 0.00019603999466879846570;
    const double c2 = 3.8426659205114376860e-08;
    const double c4 = 1.0041592839497643554e-11;
    const double c6 = 2.9516743763500191289e-15;
    const double p  = c0 + y2 *(-c2 + y2*(c4 - y2*c6));
    double sum = 0.0;
    
    int n;
    for(n=1; n<=M; n++) {
      sum += 1.0/(n * (n*n + y*y));
    }
    
    *result = -M_EULER + y2 * (sum + p);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_psi_int_e(const int n, double * result)
{
  int status = gsl_sf_psi_int_impl(n, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_psi_int_e", status);
  }
  return status;
}

int gsl_sf_psi_e(const double x, double * result)
{
  int status = gsl_sf_psi_impl(x, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_psi_e", status);
  }
  return status;
}

int gsl_sf_psi_1piy_e(const double x, double * result)
{
  int status = gsl_sf_psi_1piy_impl(x, result);
  if(status != GSL_SUCCESS){
    GSL_ERROR("gsl_sf_psi_1piy_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_psi_int(const int n)
{
  double y;
  int status = gsl_sf_psi_int_impl(n, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_psi_int", status);
  }
  return y;
}

double gsl_sf_psi(const double x)
{
  double y;
  int status = gsl_sf_psi_impl(x, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_psi", status);
  }
  return y;
}

double gsl_sf_psi_1piy(const double x)
{
  double y;
  int status = gsl_sf_psi_1piy_impl(x, &y);
  if(status != GSL_SUCCESS){
    GSL_WARNING("gsl_sf_psi_1piy", status);
  }
  return y;
}
