static void
sort_results (gsl_integration_workspace * workspace);

static void
sort_results (gsl_integration_workspace * workspace)
{
  size_t i;
  
  double * elist = workspace->elist ;
  size_t * order = workspace->order ;

  size_t nint = workspace->size;

  for (i = 0; i < nint; i++)
    {
      size_t i1 = order[i];
      double e1 = elist[i1];
      size_t i_max = i1;
      size_t j;

      for (j = i + 1; j < nint; j++)
	{
	  size_t i2 = order[j];
	  double e2 = elist[i2];

	  if (e2 >= e1)
	    {
	      i_max = i2;
	      e1 = e2;
	    }
	}

      if (i_max != i1)
	{
	  order[i] = order[i_max];
	  order[i_max] = i1;
	}
    }

  workspace->i = order[0] ;
}


