#include <stdio.h>
#include <gsl/gsl_multiset.h>

void
print_all_multisets ( size_t n, size_t k );

int
main (void)
{
  gsl_multiset * c;
  size_t i;

  printf("all multisets of {0,1,2,3} by size (lex. order)\n") ;
  for(i = 0; i <= 4; i++)
    {
      c = gsl_multiset_calloc (4, i);
      do
        {
          printf("{");
          gsl_multiset_fprintf (stdout, c, " %u");
          printf(" }\n");
        }
      while (gsl_multiset_next(c) == GSL_SUCCESS);
      gsl_multiset_free(c);
    }
  printf("all multisets of {1,2,3,4} by size (reverse lex. order)\n") ;
  for(i = 0; i <= 4; i++)
    {
      c = gsl_multiset_alloc (4, i) ;
      gsl_multiset_init_last(c) ;
      do
        {
          printf("{");
          gsl_multiset_fprintf (stdout, c, " %u");
          printf(" }\n");
        }
      while (gsl_multiset_prev(c) == GSL_SUCCESS);
      gsl_multiset_free(c);
    }
  printf("\n");

  print_all_multisets(5, 3);
  print_all_multisets(5, 0);
  print_all_multisets(5, 5);
  print_all_multisets(1, 1);
  print_all_multisets(3, 1);

  return 0;
}

void
print_all_multisets (size_t n, size_t k)
{
  gsl_multiset * c = gsl_multiset_calloc (n, k);

  printf("multisets %u choose %u (with replacement)\n", n, k);
  do
    {
      gsl_multiset_fprintf (stdout, c, " %u");
      printf("\n");
    }
  while (gsl_multiset_next(c) == GSL_SUCCESS);
  while (gsl_multiset_next(c) == GSL_SUCCESS);
  do
    {
      gsl_multiset_fprintf (stdout, c, " %u");
      printf("\n");
    }
  while (gsl_multiset_prev(c) == GSL_SUCCESS);
  printf("\n");
  gsl_multiset_free(c);
}
