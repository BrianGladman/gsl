
double gsl_stats_ivariance (int *array, int size);
double gsl_stats_dvariance (double *array, int size);

double gsl_stats_isd (int *array, int size);
double gsl_stats_dsd (double *array, int size);

double gsl_stats_iest_variance (int *array, int size);
double gsl_stats_dest_variance (double *array, int size);

double gsl_stats_iest_sd (int *array, int size);
double gsl_stats_dest_sd (double *array, int size);

double gsl_stats_iipvariance(int *array1, int *array2, int size1, int size2);
double gsl_stats_ddpvariance(double *array1, double *array2, int size1, int size2);

int gsl_stats_imax (int *array, int size);
double gsl_stats_dmax (double *array, int size);

int gsl_stats_imin (int *array, int size);
double gsl_stats_dmin (double *array, int size);
