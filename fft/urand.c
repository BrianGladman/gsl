double urand (void);

double urand (void) { 
  static unsigned long int x = 1;
  x = (1103515245 * x + 12345) & 0x7fffffffUL ;
  return x / 2147483648.0  ;
}
