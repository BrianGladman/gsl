int little_endian_p (void) ;
int endianness (void) ;
int double_endianness (void) ;
void setup_dynamic_endianness(int *b0,int *b1,int *b2,int *b3,
			      int *b4,int *b5,int *b6,int *b7);

void sprint_nybble(int i, char *s) ;
void sprint_byte(int i, char *s) ;

void make_float_bigendian (float * x) ;
void make_double_bigendian (double * x) ;

int determine_ieee_type (int non_zero, int exponent, int max_exponent);

