#include "problem.h"
#include "data.h"
#include "mcmc.h"

#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>

// class MetroJump {
// public: 
//   gsl_vector *state;
//   gsl_metrix *covariance;
// };


/*
  The state of the jump is a gsl_vector[14] st.
  s[0] = ell;
  s[1] = sdEll;
  s[2-4] = shift;
  s[5] = sdShift;
  s[6-8] = { A(11), A(22), A(33) };
  s[9-11] = { A(12), A(13), A(23) };
  s[12] = sdError;
 */

double 
logPostVectorized(const gsl_vector *x, void *param_obs) {
  
  // Convert from the void pointer to get observations.
  const obsmat obs = (const obsmat)param_obs;
  
  gsl_matrix *A = gsl_matrix_alloc(3, 3);
  assembleA(&gsl_vector_const_subvector(x, 6, 6).vector, A);

  double logP = logPost(obs,
                        gsl_vector_get(x, 0),
                        gsl_vector_get(x, 1),
                        &gsl_vector_const_subvector(x, 2, 3).vector,
                        gsl_vector_get(x, 5),
                        A,
                        gsl_vector_get(x, 12));

  gsl_matrix_free(A);

  if (isnan(logP)) {
    return GSL_NAN;
  } else {
    return logP;
  }
}

double 
invLogPost(const gsl_vector *x, void *param_obs) {
  return -logPostVectorized(x, param_obs);
}

void
assembleA(const gsl_vector *Abits, gsl_matrix *A) {
  // Main diagonal
  gsl_matrix_set(A, 0, 0, gsl_vector_get(Abits, 0));
  gsl_matrix_set(A, 1, 1, gsl_vector_get(Abits, 1));
  gsl_matrix_set(A, 2, 2, gsl_vector_get(Abits, 2));

  // Covariances
  double D = gsl_vector_get(Abits, 3);
  double E = gsl_vector_get(Abits, 4);
  double F = gsl_vector_get(Abits, 5);
  gsl_matrix_set(A, 0, 1, D);
  gsl_matrix_set(A, 1, 0, D);
  gsl_matrix_set(A, 0, 2, E);
  gsl_matrix_set(A, 2, 0, E);
  gsl_matrix_set(A, 1, 2, F);
  gsl_matrix_set(A, 2, 1, F);
}

void
dfInvLogPost(const gsl_vector * x, void * params, gsl_vector * grad) {
  size_t N = x->size;
  size_t i;
  double h = 0.01;
  gsl_vector *xcpy = gsl_vector_alloc(N);
  gsl_vector_memcpy(xcpy, x);
  
  double a, b, c, d;
  double xi0, dfdx;
  for (i = 0; i < N; i++) {
    xi0 = gsl_vector_get(x, i);
    
    // A
    gsl_vector_set(xcpy, i, xi0 + 2*h);
    a = invLogPost(xcpy, params);

    // B
    gsl_vector_set(xcpy, i, xi0 + h);
    b = invLogPost(xcpy, params);

    // C
    gsl_vector_set(xcpy, i, xi0 - h);
    c = invLogPost(xcpy, params);

    // D
    gsl_vector_set(xcpy, i, xi0 - 2*h);
    d = invLogPost(xcpy, params);

    dfdx = (-a + 8*b - 8*c + d)/(12*h);
    gsl_vector_set(grad, i, dfdx);

    gsl_vector_set(xcpy, i, xi0);
  }

  gsl_vector_free(xcpy);
}

void
fdfInvLogPost(const gsl_vector * x, void * params, double * f, gsl_vector * g) {
  *f = invLogPost(x, params);
  dfInvLogPost(x, params, g);
}

void
varianceEst(//const double (*f)(const gsl_vector * x, void * params), 
            const gsl_vector * x, 
            void * params,
            gsl_matrix *sigma) {
  /* Compute an estimate of the variance matrix as the inverse of the
     negative curvature at the point, which should be a mode. This
     assumes a quadratic log likelihood approximation. The curvature
     is computed numerically using a 2D 5-point stencil.
  */

  double toli = 0.0001, tolj = 0.0002;
  size_t i, j, N = x->size;
  gsl_vector *xcpy = gsl_vector_alloc(N);
  gsl_vector_memcpy(xcpy, x);
  
  // Compute the curvature
  double xi0, xj0;
  double a, b, c, d, curv_ij;
  for (i = 0; i < N; i++) {
    for (j = i; j < N; j++) {
      xi0 = gsl_vector_get(x, i);
      xj0 = gsl_vector_get(x, j);

      // Perturb the mode
      gsl_vector_set(xcpy, i, xi0 + toli);
      gsl_vector_set(xcpy, j, xj0 + tolj);
      a = logPostVectorized(xcpy, params);

      gsl_vector_set(xcpy, i, xi0 - toli);
      gsl_vector_set(xcpy, j, xj0 + tolj);
      b = logPostVectorized(xcpy, params);

      gsl_vector_set(xcpy, i, xi0 + toli);
      gsl_vector_set(xcpy, j, xj0 - tolj);
      c = logPostVectorized(xcpy, params);

      gsl_vector_set(xcpy, i, xi0 - toli);
      gsl_vector_set(xcpy, j, xj0 - tolj);
      d = logPostVectorized(xcpy, params);

      curv_ij = (a - b - c + d)/(4*toli*tolj);

      // Restore the mode
      gsl_vector_set(xcpy, i, xi0);
      gsl_vector_set(xcpy, j, xj0);
      
      // Store in curvature
      gsl_matrix_set(sigma, i, j, curv_ij);
      gsl_matrix_set(sigma, j, i, curv_ij);
    }
  }
  gsl_vector_free(xcpy);

  // Compute variance as V = [-L''(x_mode)]^-1
  gsl_matrix *lu = gsl_matrix_alloc(N, N);
  gsl_permutation *perm = gsl_permutation_alloc(N);
  int signum = 1;
  gsl_matrix_memcpy(lu, sigma);
  gsl_permutation_init(perm);

  gsl_matrix_scale(lu, -1.0);
  gsl_linalg_LU_decomp(lu, perm, &signum); // lu is now the LU-decomposition
  
  // gsl_linalg_LU_invert(lu, perm, sigma); // sigma is now the variance matrix
  // and we are done.

  gsl_matrix_free(lu);
  gsl_permutation_free(perm);
}

int
findMode(const obsmat obs, const gsl_vector *inits, gsl_vector *mode) {
  size_t N = 13;
  gsl_vector_memcpy(mode, inits);

  const gsl_multimin_fdfminimizer_type *T =
    gsl_multimin_fdfminimizer_vector_bfgs2;

  gsl_multimin_fdfminimizer *s = 
    gsl_multimin_fdfminimizer_alloc(T, N);

  gsl_multimin_function_fdf logPostF;
  logPostF.n      = N;
  logPostF.f      = &invLogPost;
  logPostF.df     = &dfInvLogPost;
  logPostF.fdf    = &fdfInvLogPost;
  logPostF.params = (void *)obs;

  gsl_multimin_fdfminimizer_set(s, &logPostF, inits, 0.01, 0.001);

  size_t iter = 0;
  int status;
  do {
    printf("+");
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);
    if (status == GSL_ENOPROG) { 
      printf("No progress...\n");
      int i;
      for (i = 0; i < 13; i++) {
        fprintf(stdout, "grad[%d]: %f \n", i, gsl_vector_get(s->gradient, i));
      }  
      mode = gsl_multimin_fdfminimizer_x(s);
      break;
    }

    printf("-");

    mode = gsl_multimin_fdfminimizer_x(s);

    status = gsl_multimin_test_gradient(s->gradient, 1e-3);

    if (status == GSL_SUCCESS) { break; }

  } 
  while (status == GSL_CONTINUE &&  iter < 1000);

  int i;
  for (i = 0; i < 13; i++) {
    fprintf(stdout, "mode[%d]: %f \n", i, gsl_vector_get(mode, i));
  }

  gsl_multimin_fdfminimizer_free(s);

  return status;
}
