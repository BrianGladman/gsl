static inline void 
qpsrt (gsl_integration_workspace * workspace);

static inline
void qpsrt (gsl_integration_workspace * workspace)
{
  const size_t last = workspace->size - 1;
  const size_t limit = workspace->limit;

  double * elist = workspace->elist;
  size_t * order = workspace->order;

  double errmax ;
  double errmin ;
  int i, k, top;

  size_t i_nrmax = workspace->nrmax;
  size_t i_maxerr = order[i_nrmax] ;
  
  /* Check whether the list contains more than two error estimates */

  if (last < 2) 
    {
      order[0] = 0 ;
      order[1] = 1 ;
      workspace->i = i_maxerr ;
      return ;
    }

  errmax = elist[i_maxerr] ;

  /* This part of the routine is only executed if, due to a difficult
     integrand, subdivision increased the error estimate. In the normal
     case the insert procedure should start after the nrmax-th largest
     error estimate. */

  while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]]) 
    {
      order[i_nrmax] = order[i_nrmax - 1] ;
      /* printf("NR copied order[%d] to order[%d]\n", i_nrmax-1, i_nrmax); */
      i_nrmax-- ;
    } 

  /* Compute the number of elements in the list to be maintained in
     descending order. This number depends on the number of
     subdivisions still allowed. */
  
  if(last < (limit/2 + 2)) 
    {
      top = last ;
    }
  else
    {
      top = limit - last + 1;
    }

  /* printf("top = %d\n", top) ; */

  
  /* Insert errmax by traversing the list top-down, starting
     comparison from the element elist(order(i_nrmax+1)). */
  
  i = i_nrmax + 1 ;
  
  while (errmax < elist[order[i]] && i < top)
    {
      order[i-1] = order[i] ;
      /* printf("DOWN copied order[%d] to order[%d]\n", i, i-1); */
      i++ ;
    }
  
  order[i-1] = i_maxerr ;
  /* printf("I-1 put %d into order[%d]\n", i_maxerr, i-1) ; */
  
  /* Insert errmin by traversing the list bottom-up */
  
  errmin = elist[last] ;
  
  k = top - 1 ;
  
  while (errmin >= elist[order[k]] && k > i - 2)
    {
      order[k+1] = order[k] ;
      /* printf("UP copied order[%d] to order[%d]\n", k, k+1); */
      k-- ;
    }
  
  order[k+1] = last ;
  /* printf("K put %d into order[%d]\n", last, k+1) ; */

  /* Set i_max and e_max */

  i_maxerr = order[i_nrmax] ;
  
  workspace->i = i_maxerr ;
  workspace->nrmax = i_nrmax ;
}


