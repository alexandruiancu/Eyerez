#include <stdio.h>
#include <check.h>
#include "../problem.h"
#include "../data.h"
#include <gsl/gsl_matrix.h>

START_TEST (test_computeTraces)
{
  int i;
  int n = 10;
  for (i = 1; i <= n; i++) {
    gsl_matrix *m = gsl_matrix_alloc(i, i);
    gsl_matrix_set_identity(m);
    int tr = trace(m, i);
    fail_unless(tr == i, 
                "Trace of %d-dimensional matrix should be %d (was %d)", 
                i, i, tr);
    gsl_matrix_free(m);
  }
}
END_TEST

START_TEST (test_noErrorPriors)
{
  fail_unless(ellLogPrior(30, 30.7, 2),
              "ellLogPrior failed to return result (or result == 0)");

  gsl_matrix *A = gsl_matrix_alloc(3, 3);
  gsl_matrix_set_identity(A);
  fail_unless(aLogPrior(A),
              "aLogPrior failed to return result (or result == 0)");
  gsl_matrix_free(A);

  gsl_vector *shift = gsl_vector_alloc(3);
  gsl_vector_set_all(shift, 1/sqrt(3));
  fail_unless(shiftLogPrior(shift, 3),
              "shiftLogPrior failed to return result (or result == 0)");
  gsl_vector_free(shift);
}
END_TEST

START_TEST (test_countlines)
{
  FILE *f = tmpfile();
  int i = 0;
  int N = 100;
  for (i = 0; i < N; i++) {
    fprintf (f, "%f,%f,%f\n", 1.0, 2.0, 3.0);
  }
  rewind(f);
  int count = countlines(f);
  fail_unless(N == count,
              "countlines counted %d lines in a %d-line file", count, N);

  fclose(f);
}
END_TEST

START_TEST (test_readObservations)
{
  obsmat os = allocObsCSV(fopen("e2.csv", "r"));
  fail_unless(os->size2 == 1000,
              "readObservationsCSV did not count 1000 lines in e2.csv (%d instead)",
              os->size2);

  gsl_matrix_free(os);
}
END_TEST

START_TEST (test_toFromRealSpace)
{
  obsmat os = allocObsCSV(fopen("e2.csv", "r"));
  size_t n = os->size1;
  size_t m = os->size2;
  obsmat temp = gsl_matrix_alloc(n, m);
  gsl_matrix_memcpy(temp, os);

  toRealSpace(temp, 30.4, -1);
  toLaserSpace(temp, 30.4, -1);
  
  size_t checks = 0;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      double a = gsl_matrix_get(temp, i, j);
      double b = gsl_matrix_get(os, i, j);
      checks += abs((long)(a-b)) < 0.001 ? 1 : 0;
    }
  }
  fail_unless(checks == n*m,
              "Not all elements were equal after passing through real space (%d/%d correct)",
              checks, n*m);

  gsl_matrix_free(temp);
  gsl_matrix_free(os);
}
END_TEST

START_TEST (test_logLik)
{
  obsmat os = allocObsCSV(fopen("e2.csv", "r"));
  gsl_matrix *A = gsl_matrix_alloc(3,3);
  gsl_vector *c = gsl_vector_alloc(3);
  double L, sdError;

  gsl_matrix_set_identity(A);
  gsl_vector_set_basis(c, 0);
  sdError = 2.0;
  L = 30.4;

  double out = logLik(os, L, c, A, sdError);
  fail_unless(out, "logLik failed to return value or value was 0 (value = %f)", out);

  gsl_matrix_free(A);
  gsl_matrix_free(os);
  gsl_vector_free(c);
  
}
END_TEST

Suite *
prior_suite(void) {
  Suite *s = suite_create("priors");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_computeTraces);
  tcase_add_test(tc_core, test_noErrorPriors);
  tcase_add_test(tc_core, test_countlines);
  tcase_add_test(tc_core, test_readObservations);
  tcase_add_test(tc_core, test_toFromRealSpace);
  tcase_add_test(tc_core, test_logLik);
  suite_add_tcase (s, tc_core);

  return s;
}

Suite *
data_suite(void) {
  Suite *s = suite_create("data methods");

  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_countlines);
  tcase_add_test(tc_core, test_readObservations);
  tcase_add_test(tc_core, test_toFromRealSpace);
  tcase_add_test(tc_core, test_logLik);
  suite_add_tcase (s, tc_core);

  return s;
}

int 
main (void) {
  int number_failed;
  SRunner *sr = srunner_create( prior_suite() );
  srunner_add_suite( sr, data_suite() );

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
