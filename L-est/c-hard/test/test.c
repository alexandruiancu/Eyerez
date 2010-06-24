#include "../problem.h"
#include "../data.h"
#include "../mcmc.h"

#include <stdio.h>
#include <check.h>
#include <gsl/gsl_matrix.h>

START_TEST (test_computeTraces)
{
  int i;
  int n = 10;
  for (i = 1; i <= n; i++) {
    gsl_matrix *m = gsl_matrix_alloc(i, i);
    gsl_matrix_set_identity(m);
    int tr = trace(m);
    fail_unless(tr == i, 
                "Trace of %d-dimensional matrix should be %d (was %d)", 
                i, i, tr);
    gsl_matrix_free(m);
  }
}
END_TEST

START_TEST (test_noErrorPriors)
{
  double out;
  out = ellLogPrior(30, 2);
  fail_if(isnan(out),
          "ellLogPrior failed to return reasonable result (result = %f).", out);

  gsl_matrix *A = gsl_matrix_alloc(3, 3);
  gsl_matrix_set_identity(A);
  out = aLogPrior(A);
  fail_if(isnan(out),
              "aLogPrior failed to return reasonable result (result = %f).", out);
  gsl_matrix_free(A);

  gsl_vector *shift = gsl_vector_alloc(3);
  gsl_vector_set_all(shift, 1/sqrt(3));
  out = shiftLogPrior(shift, 3);
  fail_if(isnan(out),
              "shiftLogPrior failed to return reasonable result (result = %f).", out);
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
  fail_if(isnan(out), "logLik failed to return a sensible value (value = %f)", out);

  gsl_matrix_free(A);
  gsl_matrix_free(os);
  gsl_vector_free(c);
  
}
END_TEST

START_TEST (test_assembleA)
{
  gsl_vector *v = gsl_vector_alloc(6);
  gsl_vector_set(v, 0, 1);
  gsl_vector_set(v, 1, 1);
  gsl_vector_set(v, 2, 1);
  gsl_vector_set(v, 3, 2);
  gsl_vector_set(v, 4, 3);
  gsl_vector_set(v, 5, 4);

  gsl_matrix *A = gsl_matrix_alloc(3, 3);
  gsl_matrix *M = gsl_matrix_alloc(3, 3);
  gsl_matrix_memcpy(M, A);

  gsl_matrix_set(M, 0, 0, 1);
  gsl_matrix_set(M, 1, 1, 1);
  gsl_matrix_set(M, 2, 2, 1);
  gsl_matrix_set(M, 0, 1, 2);
  gsl_matrix_set(M, 1, 0, 2);
  gsl_matrix_set(M, 0, 2, 3);
  gsl_matrix_set(M, 2, 0, 3);
  gsl_matrix_set(M, 1, 2, 4);
  gsl_matrix_set(M, 2, 1, 4);

  assembleA(v, A);

  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      fail_unless(gsl_matrix_get(A, i, j) == gsl_matrix_get(M, i, j),
                  "Failed to correctly assemble the matrix A");
    }
  }

  gsl_vector_free(v);
  gsl_matrix_free(A);
}
END_TEST

START_TEST (test_findMode)
{
  obsmat os = allocObsCSV(fopen("e1.csv", "r"));
  size_t N = 13;
  gsl_vector *inits, *mode;
  inits = gsl_vector_alloc(N);
  mode = gsl_vector_alloc(N);

  gsl_vector_set(inits, 0, 29);
  gsl_vector_set(inits, 1, 1.0);
  gsl_vector_set(inits, 2, 0.0);
  gsl_vector_set(inits, 3, 0.0);
  gsl_vector_set(inits, 4, 0.0);
  gsl_vector_set(inits, 5, 1.0);
  gsl_vector_set(inits, 6, 3.1);
  gsl_vector_set(inits, 7, 3.1);
  gsl_vector_set(inits, 8, 3.1);
  gsl_vector_set(inits, 9, 0.001);
  gsl_vector_set(inits, 10, 0.001);
  gsl_vector_set(inits, 11, 0.001);
  gsl_vector_set(inits, 12, 0.01);

  printf("invLogPost --> %f\n\n", invLogPost(inits, (void *)os));

  int i;
  gsl_vector *grad = gsl_vector_alloc(N);
  dfInvLogPost(inits, (void *)os, grad);
  for (i = 0; i < 13; i++) {
    fprintf(stdout, "grad[%d]: %f \n", i, gsl_vector_get(grad, i));
  }  
  gsl_vector_free(grad);

  for (i = 0; i < 13; i++) {
    fprintf(stdout, "inits[%d]: %f \n", i, gsl_vector_get(inits, i));
  }

  int status = findMode(os, inits, mode);
  fail_if(status, "Failed to converge finding mode.");

  gsl_vector_free(inits);
  gsl_vector_free(mode);
  gsl_matrix_free(os);
}
END_TEST

START_TEST (test_curvature)
{
  obsmat os = allocObsCSV(fopen("e1.csv", "r"));

  size_t N = 13;
  gsl_vector *mode = gsl_vector_alloc(N);
  gsl_matrix *sigma = gsl_matrix_alloc(N, N);

  gsl_vector_set(mode, 0, 30.4);
  gsl_vector_set(mode, 1, 1.0);
  gsl_vector_set(mode, 2, 0.0);
  gsl_vector_set(mode, 3, 0.0);
  gsl_vector_set(mode, 4, 0.0);
  gsl_vector_set(mode, 5, 0.1);
  gsl_vector_set(mode, 6, 3.0);
  gsl_vector_set(mode, 7, 3.0);
  gsl_vector_set(mode, 8, 3.0);
  gsl_vector_set(mode, 9, 0.0);
  gsl_vector_set(mode, 10, 0.0);
  gsl_vector_set(mode, 11, 0.0);
  gsl_vector_set(mode, 12, 0.2);

  printf("logPost --> %f\n\n", logPostVectorized(mode, (void *)os));

  varianceEst(mode, (void *)os, sigma);

  size_t i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      printf("% 9.3f ", gsl_matrix_get(sigma, i, j));
    }
    printf("\n");
  }

  gsl_vector_free(mode);
  gsl_matrix_free(sigma);
  gsl_matrix_free(os);
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

Suite *
mcmc_suite(void) {
  Suite *s = suite_create("MCMC");

  TCase *tc_core = tcase_create("Core");
  tcase_set_timeout(tc_core, 30);
  tcase_add_test(tc_core, test_assembleA);
  tcase_add_test(tc_core, test_findMode);
  tcase_add_test(tc_core, test_curvature);

  suite_add_tcase (s, tc_core);

  return s;
}


int 
main (void) {
  int number_failed;
  SRunner *sr = srunner_create( prior_suite() );
  srunner_add_suite( sr, data_suite() );
  srunner_add_suite( sr, mcmc_suite() );

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
