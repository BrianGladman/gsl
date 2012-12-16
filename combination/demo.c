#include <stdio.h>
#include <gsl/gsl_combination.h>


void
print_all_combinations ( size_t n, size_t k );

int 
main (void) 
{
  gsl_combination * c;
  size_t i;

  printf("all subsets of {0,1,2,3} by size (lex. order)\n") ;
  for(i = 0; i <= 4; i++)
    {
      c = gsl_combination_calloc (4, i);
      do
        {
          printf("{");
          gsl_combination_fprintf (stdout, c, " %u");
          printf(" }\n");
        }
      while (gsl_combination_next(c) == GSL_SUCCESS);
      gsl_combination_free(c);
    }
  printf("all subsets of {1,2,3,4} by size (reverse lex. order)\n") ;
  for(i = 0; i <= 4; i++)
    {
      c = gsl_combination_alloc (4, i) ;
      gsl_combination_init_last(c) ;
      do
        {
          printf("{");
          gsl_combination_fprintf (stdout, c, " %u");
          printf(" }\n");
        }
      while (gsl_combination_prev(c) == GSL_SUCCESS);
      gsl_combination_free(c);
    }
  printf("\n");

  print_all_combinations(5, 3);
  print_all_combinations(5, 0);
  print_all_combinations(5, 5);
  print_all_combinations(1, 1);
  print_all_combinations(3, 1);

  return 0;
}

void
print_all_combinations (size_t n, size_t k)
{
  gsl_combination * c = gsl_combination_calloc (n, k);

  printf("combinations %u choose %u\n", n, k);
  do
    {
      gsl_combination_fprintf (stdout, c, " %u");
      printf("\n");
    }
  while (gsl_combination_next(c) == GSL_SUCCESS);
  do
    {
      gsl_combination_fprintf (stdout, c, " %u");
      printf("\n");
    }
  while (gsl_combination_prev(c) == GSL_SUCCESS);
  printf("\n");
  gsl_combination_free(c);
}
