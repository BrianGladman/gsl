/* The smallest interval has the largest error.  Before bisecting
   decrease the sum of the errors over the larger intervals
   (error_over_large_intervals) and perform extrapolation. */

static int
increase_nrmax (gsl_integration_workspace * workspace);

static int
increase_nrmax (gsl_integration_workspace * workspace)
{
  int k;
  int id = workspace->nrmax;
  int jupbnd;

  const size_t * level = workspace->level;
  const size_t * order = workspace->order;

  size_t limit = workspace->limit ;
  size_t last = workspace->size - 1 ;

  if (last > (1 + limit / 2))
    {
      jupbnd = limit + 1 - last;
    }
  else
    {
      jupbnd = last;
    }
  
  for (k = id; k <= jupbnd; k++)
    {
      size_t i_max = order[workspace->nrmax];
      
      workspace->i = i_max ;

      if (level[i_max] < workspace->maximum_level)
	{
	  return 1;
	}

      workspace->nrmax++;

      /* FIXME, shouldn't we adjust workspace->i here to match nrmax? */
    }
  return 0;
}

static int
large_interval (gsl_integration_workspace * workspace)
{
  size_t i = workspace->i ;
  const size_t * level = workspace->level;
  
  if (level[i] < workspace->maximum_level)
    {
      return 1 ;
    }
  else
    {
      return 0 ;
    }
}

