static int little_endian_p (void) ;

static void sprint_nybble(int i, char *s) ;
static void sprint_byte(int i, char *s) ;

static void make_float_bigendian (float * x) ;
static void make_double_bigendian (double * x) ;

static int determine_ieee_type (int non_zero, int exponent, int max_exponent);

