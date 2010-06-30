#include "../problem.h"
#include "../data.h"

#include <gtest/gtest.h>

#include <stdio.h>
#include <gsl/gsl_matrix.h>

TEST (mathy, compute_traces) {
  int i;
  int n = 10;
  for (i = 1; i <= n; i++) {
    gsl_matrix *m = gsl_matrix_alloc(i, i);
    gsl_matrix_set_identity(m);
    int tr = trace(m);
    ASSERT_EQ(i, tr);
    gsl_matrix_free(m);
  }
}


TEST (priors, all_priors_run)
{
  double out;
  out = ellLogPrior(30, 2);
  if (isnan(out)) FAIL() << \
                    "ellLogPrior failed to return reasonable result";

  gsl_matrix *A = gsl_matrix_alloc(3, 3);
  gsl_matrix_set_identity(A);
  out = aLogPrior(A);
  if(isnan(out)) FAIL() << \
                   "aLogPrior failed to return reasonable result";
  gsl_matrix_free(A);

  gsl_vector *shift = gsl_vector_alloc(3);
  gsl_vector_set_all(shift, 1/sqrt(3));
  out = shiftLogPrior(shift, 3);
  if(isnan(out)) FAIL() << \
                   "shiftLogPrior failed to return reasonable result";
  gsl_vector_free(shift);
}


TEST (data, countlines_manual)
{
  FILE *f = tmpfile();
  int i = 0;
  int N = 100;
  for (i = 0; i < N; i++) {
    fprintf (f, "%f,%f,%f\n", 1.0, 2.0, 3.0);
  }
  rewind(f);
  int count = countlines(f);
  ASSERT_EQ(N, count);

  fclose(f);
}

TEST (data, observations_read_csv)
{
  gsl_matrix *os = allocObsCSV(fopen("e2.csv", "r"));
  ASSERT_EQ(1000, os->size2);

  gsl_matrix_free(os);
}

TEST (data, to_and_from_real_space)
{
  gsl_matrix *os = allocObsCSV(fopen("e2.csv", "r"));
  size_t n = os->size1;
  size_t m = os->size2;
  gsl_matrix *temp = gsl_matrix_alloc(n, m);
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
  ASSERT_EQ(n*m, checks);

  gsl_matrix_free(temp);
  gsl_matrix_free(os);
}

TEST (mathy, test_logLik)
{
  gsl_matrix *os = allocObsCSV(fopen("e2.csv", "r"));
  gsl_matrix *A = gsl_matrix_alloc(3,3);
  gsl_vector *c = gsl_vector_alloc(3);
  double L, sdError;

  gsl_matrix_set_identity(A);
  gsl_vector_set_basis(c, 0);
  sdError = 2.0;
  L = 30.4;

  double out = logLik(os, L, c, A, sdError);
  if(isnan(out)) FAIL() << \
                   "logLik failed to return a sensible value";

  gsl_matrix_free(A);
  gsl_matrix_free(os);
  gsl_vector_free(c);
  
}

TEST (mathy, assemble_A)
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
      ASSERT_EQ(gsl_matrix_get(M, i, j), gsl_matrix_get(A, i, j));
    }
  }

  gsl_vector_free(v);
  gsl_matrix_free(A);
}

// double curvatureCheckFn(const gsl_vector *x, void *params) {
//   return gsl_ran_bivariate_gaussian_pdf(gsl_vector_get(x, 0),
//                                         gsl_vector_get(x, 1),
//                                         1.0, 1.0, 0.0);
// }

// TEST (mathy, curvature_simple)
// {
//   size_t N = 2;
//   gsl_vector *mode = gsl_vector_alloc(N);
//   gsl_matrix *sigma = gsl_matrix_alloc(N, N);

//   gsl_vector_set(mode, 0, 0.0);
//   gsl_vector_set(mode, 1, 0.0);

//   printf("curvatureSimple --> %f\n\n", curvatureCheckFn(mode, NULL));
//   varianceEst(curvatureCheckFn, mode, NULL, sigma);

//   size_t i, j;
//   for (i = 0; i < N; i++) {
//     for (j = 0; j < N; j++) {
//       printf("% 7.3f ", gsl_matrix_get(sigma, i, j));
//     }
//     printf("\n");
//   }
// }

// TEST (mathy, curvature)
// {
//   gsl_matrix *os = allocObsCSV(fopen("e1.csv", "r"));

//   size_t N = 13;
//   gsl_vector *mode = gsl_vector_alloc(N);
//   gsl_matrix *sigma = gsl_matrix_alloc(N, N);

//   gsl_vector_set(mode, 0, 30.4);
//   gsl_vector_set(mode, 1, 0.1);
//   gsl_vector_set(mode, 2, 0.0);
//   gsl_vector_set(mode, 3, 0.0);
//   gsl_vector_set(mode, 4, 0.0);
//   gsl_vector_set(mode, 5, 0.1);
//   gsl_vector_set(mode, 6, 1/9.0);
//   gsl_vector_set(mode, 7, 1/9.0);
//   gsl_vector_set(mode, 8, 1/9.0);
//   gsl_vector_set(mode, 9, 0.0);
//   gsl_vector_set(mode, 10, 0.0);
//   gsl_vector_set(mode, 11, 0.0);
//   gsl_vector_set(mode, 12, 0.2);

//   printf("logPost --> %f\n\n", logPostVectorized(mode, (void *)os));

//   varianceEst(logPostVectorized, mode, (void *)os, sigma);

//   size_t i, j;
//   for (i = 0; i < N; i++) {
//     for (j = 0; j < N; j++) {
//       printf("% 13.3f ", gsl_matrix_get(sigma, i, j));
//     }
//     printf("\n");
//   }

//   gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
//   int err = gsl_linalg_cholesky_decomp(sigma);
//   EXPECT_NE(GSL_EDOM, err) << "Variance is not positive-definite.";
//   gsl_set_error_handler(old_handler);

//   gsl_vector_free(mode);
//   gsl_matrix_free(sigma);
//   gsl_matrix_free(os);
// }
