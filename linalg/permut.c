/* In-place Permutations 

   permute:    OUT[i]       = IN[perm[i]]     i = 0 .. N-1
   invpermute: OUT[perm[i]] = IN[i]           i = 0 .. N-1

   PERM is an index map, i.e. a vector which contains a permutation of
   the integers 0 .. N-1.

   From Knuth "Sorting and Searching", Volume 3 (3rd ed), Section 5.2
   Exercise 10 (answers), p 617

   FIXME: these have not been fully tested.

*/

static int permute (const gsl_vector_int * perm, gsl_vector * v);
static int invpermute (const gsl_vector_int * perm, gsl_vector * v);

static int
permute (const gsl_vector_int * perm, gsl_vector * v)
{
  int i, k, pk;
  double r1, t;

  size_t N = v->size;

  for (i = 0; i < N; i++)
    {
      k = gsl_vector_int_get(perm,i);

      while (k > i) 
        k = gsl_vector_int_get(perm,k);

      if (k < i)
        continue ;
      
      /* Now have k == i, i.e the least in its cycle */

      pk = gsl_vector_int_get(perm,k);

      if (pk == i)
        continue ;

      /* shuffle the elements of the cycle */

      t = gsl_vector_get(v,i);

      while (pk != i)
        {
          r1 = gsl_vector_get(v,pk);
          gsl_vector_set(v,k,r1);
          k = pk;
          pk = gsl_vector_int_get(perm,k);
        };

      gsl_vector_set(v,k,t);
    }

  return GSL_SUCCESS;
}

static int
invpermute (const gsl_vector_int * perm, gsl_vector * v)
{
  int i, k, pk;
  double r1, t;

  size_t N = v->size;

  for (i = 0; i < N; i++)
    {
      k = gsl_vector_int_get(perm,i);

      while (k > i) 
        k = gsl_vector_int_get(perm,k);

      if (k < i)
        continue ;
      
      /* Now have k == i, i.e the least in its cycle */

      pk = gsl_vector_int_get(perm,k);

      if (pk == i)
        continue ;

      /* shuffle the elements of the cycle in the inverse direction */

      t = gsl_vector_get(v,k);

      while (pk != i)
        {
          r1 = gsl_vector_get(v,pk);
          gsl_vector_set(v,pk,t);
          k = pk;
          pk = gsl_vector_int_get(perm,k);
          t = r1;
        };

      gsl_vector_set(v,pk,t);
    }

  return GSL_SUCCESS;
}
