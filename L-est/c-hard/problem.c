#include "problem.h"
#include "data.h"

#include <math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

const double SCALE = 1;

double 
logPost(const obsmat os, 
        const double ell,
        const double sdEll,
        const gsl_vector *c,
        const double sdShift,
        const gsl_matrix *A, 
        const double sdError) 
{
  double logP = 0;
  logP += logLik(os, ell, c, A, sdError);
  logP += aLogPrior(A);
  logP += ellLogPrior(ell, sdEll);
  logP += shiftLogPrior(c, sdShift);
  return logP;
}

double 
logLik(const obsmat os, 
       const double ell, 
       const gsl_vector *c, 
       const gsl_matrix *A, 
       const double sdError) 
{
  obsmat oscpy = gsl_matrix_alloc(os->size1, os->size2);
  gsl_matrix_memcpy(oscpy, os);
  toRealSpace(oscpy, ell, 1);

  size_t N = os->size2;
  size_t i;
  double logL = 0;

  for (i = 0; i < N; i++) {
    gsl_vector *o = gsl_vector_alloc(3);
    gsl_vector_const_view view = gsl_matrix_const_column(oscpy, i);
    gsl_vector_memcpy(o, &(&view)->vector);
    
    double d = computeDistance(o, c, A);
    d = log(gsl_ran_gaussian_pdf(d, 10));
    
    logL += SCALE*d;
  }

  gsl_matrix_free(oscpy);
  return logL;
}

double
computeDistance(const gsl_vector *ob, const gsl_vector *c, const gsl_matrix *A) {
  gsl_vector *obcpy = gsl_vector_alloc(3);
  gsl_vector *half = gsl_vector_alloc(3);

  gsl_vector_memcpy(obcpy, ob);
  gsl_vector_add(obcpy, c);

  // put Ax in half
  gsl_blas_dsymv(CblasUpper, 1.0, A, obcpy, 0.0, half);
  
  double dist;
  gsl_blas_ddot(obcpy, half, &dist);

  gsl_vector_free(obcpy);
  gsl_vector_free(half);
  
  return (dist - 1.0);
}

double
ellLogPrior(const double ell, const double sdEll) {
  double ALIGNMU = 30.4055916;
  double ALIGNSD = 3; 
  
  double logL = 0;

  // log P(L | mu, sd) ...
  logL += log(gsl_ran_gaussian_pdf(ALIGNMU - ell, sdEll));

  // ... + log P(sd)
  logL -= log(sdEll);

  return SCALE*logL;
}

double
trace(const gsl_matrix * X) {
  size_t n = X->size1;
  double tr = 0;
  size_t i;
  for (i = 0; i < n; i++) {
    tr += gsl_matrix_get(X, i, i);
  }
  return tr;
}

double 
aLogPrior(const gsl_matrix * A) {
  double EXP_R = 3;
  double DOF = 30;
  double logL = 0;
  
  gsl_permutation * perm = gsl_permutation_calloc(3);
  int s = 0;
  
  gsl_matrix *Acpy = gsl_matrix_alloc(3, 3);
  gsl_matrix *lu = gsl_matrix_alloc(3, 3);
  gsl_matrix_memcpy(Acpy, A);
  gsl_matrix_memcpy(lu, A);

  gsl_linalg_LU_decomp(lu, perm, &s);
  gsl_matrix_scale(Acpy, 1/EXP_R);

  double tr = trace(Acpy);
  double lndet = gsl_linalg_LU_lndet(lu);

  logL += (DOF - 4)/2 * lndet;
  logL += (-1/2) * tr;

  gsl_matrix_free(Acpy);
  gsl_matrix_free(lu);

  return SCALE*logL;
}

double
shiftLogPrior(const gsl_vector * shift, const double sdShift) {
  double logL = 0;
  size_t i;
  for (i = 0; i < 3; i++) {
    logL += log(gsl_ran_gaussian_pdf(gsl_vector_get(shift, i), sdShift));
  }
  
  logL -= log(sdShift);
  return SCALE*logL;
}
