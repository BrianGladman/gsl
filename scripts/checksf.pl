#!/usr/bin/perl

while(<DATA>) { $prog .= $_; }
print $prog;

sub bye { die; };

$SIG{'INT'} = \&bye;

# provide input file which lists functions, 1 per line in the form
#   gsl_sf_zeta_e(x, &result);

while (<>) {
    next if /^#/;
    chomp;
    ($f,$a) = /(\w+)\((.*)\)/;
    &hunt($f,$a);
};
exit;

sub hunt {
    my ($fn, $args) = @_;

    my $file="/tmp/bughunt.$$.c";
    my $exe="/tmp/bughunt.$$";
    open(FILE, ">$file");
    $p = $prog;
    $p =~ s/FN/$fn/g;
    $p =~ s/ARGS/$args/g;
    print FILE $p;
    close(FILE);
    print "compiling for $fn($args)....\n";
    sys("gcc -g -O2 -o $exe $file -lgsl -lgslcblas -lm");
    sys("GSL_IEEE_MODE=double-precision,mask-all $exe > $fn.dat;");
    unlink ("$fn.dat") if ! -s "$fn.dat";
}

exit;

sub sys {
    my ($s) = @_;
    print $s, "\n";
    system($s);
}

__END__
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_ieee_utils.h>

inline static double
extrapolate (double x0, double y0, double x1, double y1, double x2, double y2, double x)
{
    double u = x1 - x0;
    double v = x2 - x0;
    double w = x - x0;

    double du = y1 - y0;
    double dv = y2 - y0;

    double dw = w*(du*v*(w - v) + dv*u*(u - w)) / (u*v*(u - v));

    return y0 + dw;
}

double
oscillating (double x0, double x1)
{
    int i, s =0, s1=0,s2=0, c=0;
    int n = 100;
    for (i = -3*n; i <= 4*n; i++)
    {
	gsl_sf_result result;
	double x = x0 + (x1-x0) * ((double)i)/n;
	int status = FN(ARGS);	
	s2 = s1;
	s1 = s;
	s = result.val > 0 ? 1 : -1;
	if (s != s1 || s1 != s2) c++;
    }
    return (double)c/(double)n;
}

int
display (double x0, double x1)
{
    int i;
    for (i = -1000; i <= 1000; i++)
    {
	gsl_sf_result result;
	double x = x0 + (x1-x0) * ((double)i)/1000;
	int status = FN(ARGS);	
	if (status == 0)
	    printf("%.18e %.18e %g\n", x, result.val, result.err);
    }
}

int
main ()
{
  double x_min = 1e-100, x_max = 1e100;
  double x_0, x_1, x_2, x_3;
  double x, dx = 0.001;

  size_t count = 0;
  double y_1, y_2, y_3, y;

  gsl_ieee_env_setup();
  gsl_set_error_handler_off();

  fprintf(stderr, "scanning at resolution %g\n", dx);

  x = -x_max;

  while (x < x_max)
  {
      gsl_sf_result result;

      if (x < -x_min)
	  x /= 1.0+dx;
      else if (x > 0)
	  x *= 1.0 + dx;
      else if (x < 0) {
	  x = 0.0;
	  count = 0;
      } else
	  x = x_min;
      
      result.val = count % 2;
      result.err = count % 2;

      {
	  int status = FN(ARGS);
	  if (status)
	      continue;
      }

      if (count % 100000 == 0)
	  fprintf(stderr, "x = %g delta = %g\n", x, x*dx);

      y_3 = y_2; x_3 = x_2;
      y_2 = y_1; x_2 = x_1;
      y_1 = y; x_1 = x_0;
      y = result.val; x_0 = x;
      
      {
	  double y_ext = extrapolate (x_3, y_3, x_2, y_2, x_1, y_1, x);
	  double delta = fabs(y_1-y_2) + fabs(y_2-y_3) + GSL_DBL_EPSILON * (fabs(y_1)+fabs(y_2)+fabs(y_3));

	  if (count > 2 && fabs(y_ext - y) > 100*delta)
      {
	  if (oscillating (x_1, x) > 0.1)
	      continue;
          fprintf(stderr,"FAIL -1: %g %.18e\n", x_2, y_2);
          fprintf(stderr,"FAIL  0: %g %.18e\n", x_1, y_1);
          fprintf(stderr, "FAIL +1: %g %.18e\n", x_0, y);
	  display(x_1,x_0);
	  exit(1);
      }
      }
      count++;
  }
}


