static unsigned long int lcg_seed (unsigned long int * seed) ;

static unsigned long int lcg_seed (unsigned long int * seed)
{
  long int x = * seed ;
  
  const long int m = 2147483647, a = 16807, q = 127773, r = 2836 ;
  const long int h  = x / q;    

  long int t = a * (x - h * q) - h * r;

  if (t < 0) 
    {
      t += m;
    }
  
  * seed = t ;

  return t ;
}
