$fn="gsl_sf_hyperg_1F1_e";
$args = "-2.05, 1.0, x";

$file="/tmp/bughunt.$$.c";
$exe="/tmp/bughunt.$$";
open(FILE, ">$file");
while (<DATA>)
{
    s/FN/$fn/;
    s/ARGS/$args/;
    print FILE $_;
}
close(FILE);
system("cd /tmp; gcc -v -g -O2 $file -o $exe");
system("$exe");

exit;

__END__
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_sf.h>

int
main ()
{
  double x_min = 1e-100, x_max = 1e100;
  
  double x, dx = 0.001;

  size_t count = 0;
  double y_1, y_2, y;
  
  printf("scanning at resolution %g\n", dx);
  
  for (x = x_min; x < x_max; x *= 1.0+dx)
  {
      gsl_sf_result result;
      int status = FN(ARGS,&result);
      y_2 = y_1;
      y_1 = y;
      y = result.val;
      
      if (count > 2 && fabs(y - y_1) > 100*fabs(y_1 - y_2))
      {
          printf("%g %g %g %g\n", x, y_2, y_1, y);
      }
      
      count++;
  }
}
