static int little_endian_p (void) ;
static int endianness (void) ;
static int double_endianness (void) ;
static void setup_dynamic_endianness(int *b0,int *b1,int *b2,int *b3,
			      int *b4,int *b5,int *b6,int *b7);

static void sprint_nybble(int i, char *s) ;
static void sprint_byte(int i, char *s) ;

static void make_float_bigendian (float * x) ;
static void make_double_bigendian (double * x) ;

static int determine_ieee_type (int non_zero, int exponent, int max_exponent);

