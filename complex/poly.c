/* poly.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Jorma Olavi Tähtinen
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Muller-Newton root finder */
int
gsl_complex_poly_solve(gsl_complex *coeff, int n,
                       gsl_complex *Roots, double *maxerr);

/***********************************************************************
 * Muller-Newton polynom root solver
 ***********************************************************************/

// Muller's root Estimation algorithm (defs)
#define GSL_MULLER_ITERMAX     150
#define GSL_MULLER_CONVERGENCE 100 
#define GSL_MULLER_MAXDIST     1e3 
#define GSL_MULLER_FACTOR      1e5 
#define GSL_MULLER_KITERMAX    1e3 
#define GSL_MULLER_FVALUE      1e36 
#define GSL_MULLER_BOUND1      1.01
#define GSL_MULLER_BOUND2      0.99 
#define GSL_MULLER_BOUND3      0.01  
#define GSL_MULLER_BOUND4      sqrt(DBL_MAX)/1e4 
#define GSL_MULLER_BOUND6      log10(GSL_MULLER_BOUND4)-4 
#define GSL_MULLER_BOUND7      1e-5 
#define GSL_MULLER_NOISESTART  DBL_EPSILON*1e2 
#define GSL_MULLER_NOISEMAX    5     

// Newton method Calculates root of polynomial with complex coefficients
// Intial estimate is the root determined by Muller's method
#define GSL_NEWTON_ITERMAX  20 
#define GSL_NEWTON_FACTOR   5    
#define GSL_NEWTON_FVALUE   1E36 
#define GSL_NEWTON_BOUND    sqrt(DBL_EPSILON)
#define GSL_NEWTON_NOISEMAX 5 

struct gsl_complex_poly_workspace {
  // Workspace for Mullar Algorithm
  gsl_complex x0,x1,x2; // Common Points [x0,f(x0)=P(x0)], ... [x2,f(x2)]
  gsl_complex f0,f1,f2; // Of Parabola And Polynomial
  gsl_complex h1,h2;    // Distance Between X2 And x1
  gsl_complex q2;       // Smaller Root Of Parabola
  int iter;             // Iteration Counter
};

/*----------------------------------------------------------------------
  make monic (1*x^n+coeff(n-1)*x^(n-1)+coeff(n-2)*x^(n-2)+...)
*/
void gsl_complex_poly_makemonic(gsl_complex *coeff, int N) {
  register int n;
  while(--N>=0 && gsl_complex_abs(&coeff[N])==0.0);
  if (N>=0) {
    double factor=1./gsl_complex_abs(&coeff[N]);
    if (factor!=1.0) {
      for(n=0;n<=N;n++) gsl_complex_mul_d(&coeff[n],&coeff[n],factor);
    }
  }
}

/*----------------------------------------------------------------------
  Calc the roots(R[i]) of a quadratic equation: 
  coeff[2]x^2+coeff[1]x+coeff[0]==0
*/
void gsl_complex_poly_quadratic(gsl_complex *coeff,gsl_complex *R) {
  gsl_complex D,N;
  gsl_complex_mul(&D,&coeff[1],&coeff[1]);
  gsl_complex_mul(&N,&coeff[2],&coeff[0]);
  gsl_complex_mul_d(&N,&N,4.0);
  gsl_complex_sub(&D,&D,&N);
  gsl_complex_sqrt(&D,&D);      /* D=Sqrt(coeff[1]^2-4.0*coeff[2]*coeff[0]); */
  gsl_complex_mul_d(&N,&coeff[2],2.0); /* N=2.0*coeff[2] */
  gsl_complex_mul_d(&R[0],&coeff[1],-1.0);
  gsl_complex_set(&R[1],&R[0]);
  gsl_complex_add(&R[0],&R[0],&D); 
  gsl_complex_sub(&R[1],&R[1],&D);
  gsl_complex_div(&R[0],&R[0],&N); /* R[0]=(-coeff[1]+D)/N; */
  gsl_complex_div(&R[1],&R[1],&N); /* R[1]=(-coeff[1]-D)/N; */
}

/*----------------------------------------------------------------------
  Calc The Roots Of Linear/Quadratic Equation (Lin: a[1]x+a[0]==0)
*/
int gsl_complex_poly_linquad(gsl_complex *coeff,int n,gsl_complex *R) {
  if (n==2) { /* Linear Equation (really hard to solve) */
    gsl_complex_mul_d(&R[0],&coeff[0],-1.0);
    gsl_complex_div(&R[0],&R[0],&coeff[1]); /* R[0]=-coeff[0]/coeff[1]; */
    return TRUE;
  } else if (n==3) { // Quadratic Use Predefined formulae
    gsl_complex_poly_quadratic(coeff,R);
    return TRUE;
  }  
  return FALSE; /* not Linear or Quadratic */
}

/*----------------------------------------------------------------------
  Horner Method To Deflate One Root
*/
void gsl_complex_poly_hornc(gsl_complex *coeff,int n,gsl_complex *x0,int flag) {
  register int k;
  if (flag&0x1) { /* complex coefficients */
    gsl_complex tmp;
    for(k=n-2;k>=0;--k) {
      gsl_complex_mul(&tmp,&coeff[k+1],x0);
      gsl_complex_add(&coeff[k],&coeff[k],&tmp);
    }
  } else { /* real coefficients (minor speedup) */    
    for(k=n-2;k>=0;--k) {
      GSL_COMPLEX_P_REAL(&coeff[k])+=GSL_COMPLEX_P_REAL(x0)*GSL_COMPLEX_P_REAL(&coeff[k+1]);
    }
  }
}

/*----------------------------------------------------------------------
  Horner Method To Deflate Two Roots
*/
void gsl_complex_poly_Horncd(gsl_complex *coeff,int n,double a,double b) {
  register int k;
  GSL_REAL(coeff[n-2])+=GSL_REAL(coeff[n-1])*a;
  for(k=n-3;k>=0;--k) {
    GSL_REAL(coeff[k])+=(a*GSL_REAL(coeff[k+1])+b*GSL_REAL(coeff[k+2]));
  }
}

/*----------------------------------------------------------------------
  Main Routine To Deflate Polynomial
*/
int gsl_complex_poly_deflate(gsl_complex *coeff,int n,gsl_complex *R,int flag) {
  double a,b;     /* coefficients of the quadratic polynomial x^2+ax+b */
  gsl_complex x0; /* root to be deflated */
  gsl_complex_set(&x0,&R[n-2]);

  if (GSL_IMAG(x0)!=0.0) flag|=2; /* x0 is complex */
  if (flag==2) { /* real coefficients and complex root => deflate x0 and Conjg(x0) */
    a=2*GSL_REAL(x0);
    b=-GSL_COMPLEX_ABS2(x0);
    gsl_complex_conjugate(&R[n-3],&x0);
    gsl_complex_poly_Horncd(coeff,n,a,b);
    return 2; /* 2 roots deflated */               
  } else { /* deflate only one root */
    gsl_complex_poly_hornc(coeff,n,&x0,flag);
    return 1; /* 1 root deflated */
  }
}

/*----------------------------------------------------------------------
  Muller's root Estimation algorithm 
  L  : List of predeflated polynomial coefficients
  xb : xb Best optained x-value approximation
  e  : bound for |q2| (epsilon)
*/
void gsl_complex_poly_muller_init(struct gsl_complex_poly_workspace *w,
                                  gsl_complex *xb,double *epsilon) 
{
  GSL_SET_COMPLEX(&w->x0,0.,0.);                   /* orig: x0 = 0 + j*1     */
  GSL_SET_COMPLEX(&w->x1,-1./sqrt(2),-1./sqrt(2)); /* orig: x1 = 0 - j*1     */
  GSL_SET_COMPLEX(&w->x2,1./sqrt(2),1./sqrt(2));   /* x2 = (1 + j*1)/sqrt(2) */
  
  gsl_complex_sub(&w->h1,&w->x1,&w->x0); /* h1 = x1 - x0 */
  gsl_complex_sub(&w->h2,&w->x2,&w->x1); /* h2 = x2 - x1 */
  gsl_complex_div(&w->q2,&w->h2,&w->h1); /* q2 = h2/h1   */

  gsl_complex_set(xb,&w->x2);         /* best initial x-value = zero */
  *epsilon=GSL_MULLER_FACTOR*DBL_EPSILON; /* accuracy for determined root*/
  w->iter=0;                          /* reset iteration counter     */
}

/* Polynomial Value f=P(x0) Horner's metods (df!=NULL -> df=P'(x0)) */
void gsl_complex_poly_fdvalue(gsl_complex *coeff,int n,
                              gsl_complex *f,gsl_complex *df,gsl_complex *x0) 
{
  int i;
  gsl_complex tmp;
  if (n<=0) return;
  
  gsl_complex_set(f,&coeff[n-1]); /* f=coeff[n-1]; */
  if (df) {
    GSL_SET_COMPLEX(df,0.0,0.0);
    for(i=n-2;i>=0;--i) {
      gsl_complex_mul(&tmp,df,x0);       /* tmp=df*x   */
      gsl_complex_add(df,&tmp,f);        /* df=tmp+f   */
      gsl_complex_mul(&tmp,f,x0);        /* tmp=f*x    */
      gsl_complex_add(f,&tmp,&coeff[i]); /* f=tmp+a[i] */
    }
  } else {
    for(i=n-2;i>=0;--i) {
      gsl_complex_mul(&tmp,f,x0);        /* tmp=f*x;    */
      gsl_complex_add(f,&tmp,&coeff[i]); /* f=tmp+a[i]; */
    }
  }
}

// Calculate smaller root of Muller's parabola
void gsl_complex_poly_muller_rootofparabola(struct gsl_complex_poly_workspace *w) {
  gsl_complex A2,B2,C2; /* variables to get q2 */
  gsl_complex N1,N2;    /* denominators of q2  */
  gsl_complex D;        /* discriminante       */
  
  /* A2 = q2(f2 - (1+q2)f1 + f0q2) */
  /* B2 = q2[q2(f0-f1) + 2(f2-f1)] + (f2-f1) */
  /* C2 = (1+q2)f[2] */
  gsl_complex_add_d(&B2,&w->q2,1.0);  /* B2=1.0+q2 (B2 as a temp var) */
  gsl_complex_mul(&B2,&B2,&w->f1);    /* B2=f1*(1.0+q2) */
  gsl_complex_mul(&A2,&w->q2,&w->f0); /* A2=q2*f0 */
  gsl_complex_add(&A2,&A2,&w->f2);
  gsl_complex_sub(&A2,&A2,&B2);       /* A2=f2+q2*f0-f1*(1.0+q2)      */
  gsl_complex_mul(&A2,&A2,&w->q2);    /* A2=q2*(f2+q2*f0-f1*(1.0+q2)) */
  gsl_complex_sub(&C2,&w->f2,&w->f1); /* C2(temp var)=f2-f1 */
  gsl_complex_mul_d(&B2,&C2,2.0);     /* B2=2.0*(f2-f1) */
  gsl_complex_sub(&N1,&w->f0,&w->f1); /* N1(temp var)=f0-f1 */
  gsl_complex_mul(&N1,&N1,&w->q2);    /* N1=q2(f0-f1) */
  gsl_complex_add(&B2,&B2,&N1);       /* B2=q2(f0-f1)+2.0*(f2-f1) */
  gsl_complex_mul(&B2,&B2,&w->q2);    /* B2=q2*(q2(f0-f1)+2.0*(f2-f1)) */
  gsl_complex_add(&B2,&B2,&C2);       /* B2=q2*(q2(f0-f1)+2.0*(f2-f1))+(f2-f1) */
  gsl_complex_add_d(&C2,&w->q2,1.0);  /* C2=1.0+q2 */
  gsl_complex_mul(&C2,&C2,&w->f2);    /* C2=f2(1+q2) */

  /* D = B2^2 - 4*A2*C2 */
  gsl_complex_mul(&D,&B2,&B2);    /* D=B2*B2 */
  gsl_complex_mul(&N1,&A2,&C2);   /* N1(temp var)=A2*C2 */
  gsl_complex_mul_d(&N1,&N1,4.0); 
  gsl_complex_sub(&D,&D,&N1);     /* D=B2*B2-4.0*A2*C2;      */
  gsl_complex_sqrt(&D,&D);        /* D=sqrt(B2*B2-4.0*A2*C2; */

  /* Denominators of q2 */
  gsl_complex_sub(&N1,&B2,&D); /* N1=B2-Sqrt(D); */
  gsl_complex_add(&N2,&B2,&D); /* N2=B2+Sqrt(D); */

  /* Choose Denominater With Largest Modulus */
  if (gsl_complex_abs2(&N1)>gsl_complex_abs2(&N2) && 
      gsl_complex_abs(&N1)>DBL_EPSILON)
    {
      gsl_complex_mul_d(&w->q2,&C2,-2.0);
      gsl_complex_div(&w->q2,&w->q2,&N1); /* q2=(-2.0*C2)/N1; */
    } else if (gsl_complex_abs(&N2)>DBL_EPSILON) {
      gsl_complex_mul_d(&w->q2,&C2,-2.0);
      gsl_complex_div(&w->q2,&w->q2,&N2); /* q2=(-2.0*C2)/N2; */
    } else {
      GSL_SET_COMPLEX(&w->q2,cos((double)w->iter),sin((double)w->iter));
    }
}

/* Main Iteration Equation: x2 = h2*q2 + x2  */
/* h2abs: Absolute value of the old distance */
void gsl_complex_poly_muller_iterequ(struct gsl_complex_poly_workspace *w,
                                     double *h2abs) 
{
  double h2absnew; /* Absolute value of the new h2 */
  double help;     /* help variable */
  
  gsl_complex_mul(&w->h2,&w->h2,&w->q2); /* h2*=q2; */
  h2absnew=gsl_complex_abs(&w->h2); /* distance between old and new x2 */
  
  if (h2absnew>(*h2abs*GSL_MULLER_MAXDIST)) { /* maximum relative change */
    help=GSL_MULLER_MAXDIST/h2absnew;
    gsl_complex_mul_d(&w->h2,&w->h2,help);
    gsl_complex_mul_d(&w->q2,&w->q2,help);
  } 
  *h2abs=h2absnew; /* actualize old distance for next iteration */
  gsl_complex_add(&w->x2,&w->x2,&w->h2); /* x2+=h2; */  
}

// Suppress Overflow
// nred the highest exponent of the deflated polynomial
void gsl_complex_poly_muller_suppressoverflow(struct gsl_complex_poly_workspace *w,int nred) {
  int f=FALSE; /* Loop Ready */
  int kiter;   /* Internal Iteration Counter */
  double help; /* Help Variable */
  
  kiter=0; // Reset Iteration Counter
  do { 
    f=FALSE;                      /* Initial Estimation: No Overflow */
    help=gsl_complex_abs(&w->x2); /* help = |x2| */
    if (help>1.0 && fabs(nred*log10(help))>GSL_MULLER_BOUND6) {
      kiter++; /* if |x2|>1 and |x2|^nred>10^BOUND6 */
      if (kiter<GSL_MULLER_KITERMAX) { /* Then Halve The Distance Between */
        gsl_complex_mul_d(&w->h2,&w->h2,0.5);  /* h2*=0.5; new and old x2 */
        gsl_complex_mul_d(&w->q2,&w->q2,0.5);  /* q2*=0.5; */
        gsl_complex_sub(&w->x2,&w->x2,&w->h2); /* x2-=h2; */
        f=TRUE;
      } else 
        kiter=0;  
    }
  } while(f);
}

void gsl_complex_poly_muller_computefunc(struct gsl_complex_poly_workspace *w,
                                         gsl_complex *coeff,int n,
                                         double f1absq,double *f2absq,double epsilon)
{
  double tmp;
  int overflow; /* overflow flag */
  do {
    overflow=FALSE; /* initial estimation: no overflow */
    gsl_complex_poly_muller_suppressoverflow(w,n); /* suppress overflow */
    
    /* Calculate New Value => Result In f2 */
    gsl_complex_poly_fdvalue(coeff,n,&w->f2,NULL,&w->x2);
    
    /* Check Of Too Big Function Values f2absq=|f2|^2 */
    tmp=fabs(GSL_COMPLEX_P_REAL(&w->f2))+fabs(GSL_COMPLEX_P_IMAG(&w->f2));
    if (tmp>GSL_MULLER_BOUND4) { /* limit |f2|^2, when |f2.r|+|f2.i|>BOUND4 */
      *f2absq=tmp;
    } else {
      *f2absq=gsl_complex_abs2(&w->f2); /* |f2|^2 = f2.r^2+f2.i^2 */
    }
    
    /* Increase Iterationcounter */
    w->iter++;
    
    /* Muller's modification to improve convergence */
    if ((*f2absq>(GSL_MULLER_CONVERGENCE*f1absq)) && 
        (gsl_complex_abs(&w->q2)>epsilon) && 
        (w->iter<GSL_MULLER_ITERMAX)) 
      {
        gsl_complex_mul_d(&w->q2,&w->q2,0.5);  /* q2*=0.5; in case of overflow: */
        gsl_complex_mul_d(&w->h2,&w->h2,0.3);  /* h2*=0.3; halve q2 and h2; compute new x2 */
        gsl_complex_sub(&w->x2,&w->x2,&w->h2); /* x2-=h2; */
        overflow=TRUE;
      }
  } while(overflow);
}

/*----------------------------------------------------------------------
  calc muller approximation for newton iteration 
*/
void gsl_complex_poly_muller_calc(struct gsl_complex_poly_workspace *w,
                                  gsl_complex *coeff,int n,gsl_complex *xb) 
{  
  gsl_complex df,tmp;
  double f1absq=GSL_MULLER_FVALUE;  /* f1absq=|f1|^2 */
  double f2absq=GSL_MULLER_FVALUE;  /* f2absq=|f2|^2 */
  double f2absqb=GSL_MULLER_FVALUE; /* f2absqb=|P(xb)|^2 */
  double h2abs;                     /* h2abs=|h2| */
  double epsilon;                   /* bound for |q2| */
  int seconditer=0;                 /* second iteration, when root is too bad */
  int noise=0;                      /* noise counter */
  int rootd=FALSE;                  /* Root determined */
  
  /* Initializing Routine */
  gsl_complex_poly_muller_init(w,xb,&epsilon);

  gsl_complex_poly_fdvalue(coeff,n,&w->f0,NULL,&w->x0);
  gsl_complex_poly_fdvalue(coeff,n,&w->f1,NULL,&w->x1);
  gsl_complex_poly_fdvalue(coeff,n,&w->f2,NULL,&w->x2);
  
  do { /* loop for possible second iteration */
    do { /* main iteration loop calculate the roots of the parabola */
      gsl_complex_poly_muller_rootofparabola(w);
      
      /* store values for the next iteration */
      gsl_complex_set(&w->x0,&w->x1); /* x0=x1; */
      gsl_complex_set(&w->x1,&w->x2); /* x1=x2; */
      
      h2abs=gsl_complex_abs(&w->h2); /* Distance Between x2 and x1 */
      gsl_complex_poly_muller_iterequ(w,&h2abs); /* main iteration-equation */
      
      /* store values for the next iteration */
      gsl_complex_set(&w->f0,&w->f1); /* f0=f1; */
      gsl_complex_set(&w->f1,&w->f2); /* f1=f2; */
      f1absq=f2absq;
      
      /* compute P(x2) and make some checks */
      gsl_complex_poly_muller_computefunc(w,coeff,n,f1absq,&f2absq,epsilon);
      
      /* is the new x-value (x2) the best approximation? */
      if ((f2absq<=(GSL_MULLER_BOUND1*f1absq)) && 
          (f2absq>=(GSL_MULLER_BOUND2*f1absq))) 
        {
          /* function-value changes slowly */
          if (gsl_complex_abs(&w->h2)<GSL_MULLER_BOUND3) { 
            /* if |h[2]| is small enough => double q2 and h[2] */
            gsl_complex_mul_d(&w->q2,&w->q2,2.0); /* q2*=2.0; */
            gsl_complex_mul_d(&w->h2,&w->h2,2.0); /* h2*=2.0; */
          } else {  
            /* otherwise: |q2| = 1 and h[2]=h[2]*q2 */
            GSL_SET_COMPLEX(&w->q2,cos((double)w->iter),sin((double)w->iter));
            gsl_complex_mul(&w->h2,&w->h2,&w->q2); /* h2*=q2; */
          }
        } else if (f2absq<f2absqb) {
          f2absqb=f2absq;             /* the new function value is the */
          gsl_complex_set(xb,&w->x2); /* xb=x2 (best approximation) */
          noise=0;                    /* reset noise counter */
          
          gsl_complex_sub(&tmp,&w->x2,&w->x1);
          gsl_complex_div(&tmp,&tmp,&w->x2); /* tmp=(x2-x1)/x2 */
          if ((sqrt(f2absq)<epsilon) && (gsl_complex_abs(&tmp)<epsilon)) 
            rootd=TRUE;  /* root determined */
        }
      
      /* increase noise counter */
      if (fabs((gsl_complex_abs(xb)-
                gsl_complex_abs(&w->x2))/gsl_complex_abs(xb))<GSL_MULLER_NOISESTART) 
        noise++;
    } while ((w->iter<=GSL_MULLER_ITERMAX) && 
             (!rootd) && 
             (noise<=GSL_MULLER_NOISEMAX));
    
    seconditer++; /* increase seconditer */
     
    /* check, if determined root is good enough */
    if ((seconditer==1) && (f2absqb>0)) { 
      gsl_complex_poly_fdvalue(coeff,n,&w->f2,&df,xb); /* f2=P(x0), df=P'(x0) */      
      if (gsl_complex_abs(&w->f2)/(gsl_complex_abs(&df)*gsl_complex_abs(xb))>GSL_MULLER_BOUND7) {
        /* start second iteration with new initial estimations */
        GSL_SET_COMPLEX(&w->x0,1.0,0.0);
        GSL_SET_COMPLEX(&w->x1,-1.0,0.0);
        GSL_SET_COMPLEX(&w->x2,0.0,0.0);

        gsl_complex_poly_fdvalue(coeff,n,&w->f0,NULL,&w->x0);
        gsl_complex_poly_fdvalue(coeff,n,&w->f1,NULL,&w->x1);
        gsl_complex_poly_fdvalue(coeff,n,&w->f2,NULL,&w->x2);
        
        w->iter=0;    /* reset iteration counter */
        seconditer++; /* increase seconditer */
        rootd=FALSE;  /* no root determined */
        noise=0;      /* reset noise counter */
      }
    }
  } while (seconditer==2);  
}
 
/*----------------------------------------------------------------------
  Newton method Calculates root of polynomial with complex coefficients
  Intial estimate is the root determined by Muller's method
  xmin : best x determined in newton()
*/
void gsl_complex_poly_newton(gsl_complex *xmin,gsl_complex *L,int n,
                             gsl_complex *ns,double *dxabs,int flag) 
{
  int iter=0;  /* counter */
  int noise=0; /* noisecounter */
  double fabsmin=GSL_NEWTON_FVALUE; /* fabsmin = |P(xmin)| */
  double eps=DBL_EPSILON;
  gsl_complex x0;   /* iteration variable for x-value     */
  gsl_complex f;    /* f       = P(x0)                    */
  gsl_complex df;   /* df      = P'(x0)                   */
  gsl_complex dx;   /* dx      = P(x0)/P'(x0)             */
  gsl_complex dxh;  /* help variable dxh = P(x0)/P'(x0)   */

  gsl_complex_set(&x0,ns);      /* Initial Estimation = root determined with Muller method */
  gsl_complex_set(xmin,&x0);    /* initial estimation for the best x-value */
  GSL_SET_COMPLEX(&dx,1.0,0.0); /* initial value: P(x0)/P'(x0)=1+j*0 */  
  *dxabs=gsl_complex_abs(&dx);  /* initial value: |P(x0)/P'(x0)|=1 */
  
  for(iter=0;iter<GSL_NEWTON_ITERMAX;iter++) {
    gsl_complex_poly_fdvalue(L,n,&f,&df,&x0); /* f=P(x0), df=P'(x0) */
    if (gsl_complex_abs(&f)<fabsmin) {  /* the new x0 is a better          */
      gsl_complex_set(xmin,&x0);        /* approximation than the old xmin */
      fabsmin=gsl_complex_abs(&f);      /* store new xmin and fabsmin      */
      noise=0;                          /* reset noise counter             */
    }
    
    if (gsl_complex_abs(&df)!=0.0) { /* calculate new dx */
      gsl_complex_div(&dxh,&f,&df); /* dxh=f/df; */
      if (gsl_complex_abs(&dxh)<*dxabs*GSL_NEWTON_FACTOR) { /* new dx small enough? */
        gsl_complex_set(&dx,&dxh); /* dx=dxh, store new dx for next iteration */
        *dxabs=gsl_complex_abs(&dx);                       
      }
    }
    
    if (gsl_complex_abs(xmin)!=0.0) {
      if (*dxabs/gsl_complex_abs(xmin)<eps || noise ==GSL_NEWTON_NOISEMAX) { // routine ends 
        if (fabs(GSL_COMPLEX_P_IMAG(xmin))<GSL_NEWTON_BOUND && (!flag)) {
          /* define determined root as real if imag. part<BOUND */
          GSL_SET_IMAG(xmin,0.0); /* xmin.I=0 */
        }
        *dxabs=*dxabs/gsl_complex_abs(xmin); /* return relative error */
        return; /* return best approximation (xmin) */
      }
    } 
    gsl_complex_sub(&x0,&x0,&dx); /* main iteration: x0 = x0 - P(x0)/P'(x0) */
    noise++; /* Increase Noise Counter */
  }
    
  if (fabs(GSL_COMPLEX_P_IMAG(xmin))<GSL_NEWTON_BOUND && (!flag)) { 
    // define determined root as real, if imag. part<BOUND
    GSL_SET_IMAG(xmin,0.0); /* xmin.I=0 */
  }

  if (gsl_complex_abs(xmin)!=0.0) 
    *dxabs=*dxabs/gsl_complex_abs(xmin); /* return relative error */
  
  /*
    maximum number of iterations exceeded: xmin has best xmin until now,
    !!! maybe should warn here !!!
  */
}

int
gsl_complex_poly_solve(gsl_complex *coeff, int n,
                            gsl_complex *Roots, double *maxerr,
                            gsl_complex *coeff_deflated,
                            struct gsl_complex_poly_workspace *w,
                            int compcoef)
{
  int k,l=-1,ndef,rdef,off;
  double newerr;
  gsl_complex ns;

  if (n<=1) {
    GSL_ERROR("polynom order (n) must be >1",GSL_EINVAL);
  }
  if (gsl_complex_abs(&coeff[n-1])==0.0) {
    GSL_ERROR ("leading term of polynomial must be non-zero", GSL_EINVAL);
  }

  *maxerr=0.0;

  for(k=0;k<n;++k) gsl_complex_set(&coeff_deflated[k],&coeff[k]);
  for(k=0;k<n&&(gsl_complex_abs(&coeff[k])==0.0);++k);
  for(l=0;l<k;++l) GSL_SET_COMPLEX(&Roots[n-l-2],0.0,0.0); /* x(... k ...(x(x^def+..))) */
  ndef=n-k; off=k;

  /* polynomial is linear or quadratic */
  if (gsl_complex_poly_linquad(&coeff[off],ndef,Roots)) {  
    *maxerr=DBL_EPSILON; 
    return TRUE;
  }
  
  /* Convert To Monic */
  gsl_complex_poly_makemonic(&coeff[off],ndef);

  do {      
    /* First Muller estimate */
    gsl_complex_poly_muller_calc(w,&coeff_deflated[off],ndef,&ns);

    /* Newton method */
    gsl_complex_poly_newton(&Roots[ndef-2],coeff,n,&ns,&newerr,compcoef);

    if (newerr>*maxerr) *maxerr=newerr;
    
    /* deflate polynomial */
    rdef=gsl_complex_poly_deflate(&coeff_deflated[off],ndef,Roots,compcoef?0x1:0x0);

    off+=rdef;
    ndef-=rdef; /* reduce degree of polynomial */

  } while(ndef>2); // Is Quadratic or Linear

  gsl_complex_poly_linquad(&coeff_deflated[off],ndef,Roots);
  if (ndef==2) {
    gsl_complex_poly_newton(&Roots[1],coeff,n,&Roots[1],&newerr,compcoef); 
    if (newerr>*maxerr) *maxerr=newerr;
  }    
  gsl_complex_poly_newton(&Roots[0],coeff,n,&Roots[0],&newerr,compcoef); 
  if (newerr>*maxerr) *maxerr=newerr;

  return GSL_SUCCESS;
}

int
gsl_complex_poly_solve(gsl_complex *coeff, int n,
                       gsl_complex *Roots, double *maxerr)
{
  double dummy;
  int res,k,compcoef=FALSE;  
  gsl_complex *defspace;
  struct gsl_complex_poly_workspace w; /* Work space for muller */

  for(k=0;k<n;++k) if (GSL_IMAG(coeff[k])!=0.0) { compcoef=TRUE; break; }
  
  defspace=(gsl_complex *)malloc(n*sizeof(gsl_complex));  
  
  if (!maxerr) maxerr=&dummy;
  res=gsl_complex_poly_solve(coeff,n,Roots,maxerr,defspace,&w,compcoef);
                             
  free(defspace);
  return res;
}
