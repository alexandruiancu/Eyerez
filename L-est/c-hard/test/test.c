#include <stdio.h>
#include <check.h>
#include "../problem.h"
#include <gsl/gsl_matrix.h>

START_TEST (test_computeTraces)
{
  int i;
  int n = 10;
  for (i = 1; i <= n; i++) {
    gsl_matrix *m = gsl_matrix_alloc(i, i);
    gsl_matrix_set_identity(m);
    fail_unless(trace(m, i) == i);
    gsl_matrix_free(m);
  }
}
END_TEST

START_TEST (test_noErrorPriors)
{
  ellLogPrior(30, 30.7, 2);

  gsl_matrix *A = gsl_matrix_alloc(3, 3);
  gsl_matrix_set_identity(A);
  aLogPrior(A);
  gsl_matrix_free(A);

  gsl_vector *shift = gsl_vector_alloc(3);
  gsl_vector_set_all(shift, 1/sqrt(3));
  shiftLogPrior(shift, 3);
  gsl_vector_free(shift);
}
END_TEST

Suite *
prior_suite(void) {
  Suite *s = suite_create("priors");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_computeTraces);
  tcase_add_test(tc_core, test_noErrorPriors);
  
  suite_add_tcase (s, tc_core);

  return s;
}

int 
main (void) {
  int number_failed;
  Suite *s = prior_suite ();
  SRunner *sr = srunner_create (s);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
