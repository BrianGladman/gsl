int
FUNCTION(compare_complex,results) (const char *name_a, const BASE a[],
				   const char *name_b, const BASE b[],
				   size_t stride, size_t n, 
				   const double allowed_ticks);

int
FUNCTION(compare_real,results) (const char *name_a, const BASE a[],
				   const char *name_b, const BASE b[],
				   size_t stride, size_t n, 
				   const double allowed_ticks);
