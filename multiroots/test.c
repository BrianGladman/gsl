int foo_f (const gsl_vector * x, void * params, gsl_vector * f)
{
  if (x->size != 5 || f->size != 5)
    {
      GSL_ERROR ("foo needs a system of 5 equations", GSL_EINVAL);
    }

  {
    double x0 = gsl_vector_get (x,0);
    double x1 = gsl_vector_get (x,0);
    double x2 = gsl_vector_get (x,0);
    double x3 = gsl_vector_get (x,0);
    double x4 = gsl_vector_get (x,0);
    
    double y0 = sin(x0) + cos(x1) ;
    double y1 = sin(2*x1) + cos(x2) ;
    double y2 = sin(3*x2) + cos(x3) ;
    double y3 = sin(4*x3) + cos(x4) ;
    double y4 = sin(5*x4) + cos(x0) ;
    
    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);
    gsl_vector_set(f, 2, y2);
    gsl_vector_set(f, 3, y3);
    gsl_vector_set(f, 4, y4);
  }
}

int foo_df (const gsl_vector * x, void * params, gsl_matrix * df)
{
  if (x->size != 5 || f->size != 5)
    {
      GSL_ERROR ("foo needs a system of 5 equations", GSL_EINVAL);
    }

  {
    double x0 = gsl_vector_get (x,0);
    double x1 = gsl_vector_get (x,0);
    double x2 = gsl_vector_get (x,0);
    double x3 = gsl_vector_get (x,0);
    double x4 = gsl_vector_get (x,0);
    
    double dy0_dx0 = cos(x0) + cos(x1) ;
    double dy1 = 2*cos(2*x1) + cos(x2) ;
    double dy2 = 3*cos(3*x2) + cos(x3) ;
    double dy3 = 4*cos(4*x3) + cos(x4) ;
    double dy4 = 5*cos(5*x4) + cos(x0) ;
    
    gsl_vector_set(df, 0, dy0);
    gsl_vector_set(df, 1, dy1);
    gsl_vector_set(df, 2, dy2);
    gsl_vector_set(df, 3, dy3);
    gsl_vector_set(df, 4, dy4);
  }
}
