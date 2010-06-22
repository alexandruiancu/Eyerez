#include "problem.h"

#include <math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

double 
logLik(gsl_matrix * data, double L, gsl_vector c, gsl_matrix * A) {
  // TO IMPLEMENT
}

double
ellLogPrior(double ell, double muEll, double sdEll) {
  double ALIGNMU = 30.4055916;
  double ALIGNSD = 3; 
  
  double logL = 0;

  // log P(L | mu, sd) ...
  logL += log(gsl_ran_gaussian_pdf(muEll - ell, sdEll));

  // ... + log P(mu)
  logL += log(gsl_ran_gaussian_pdf(ALIGNMU - muEll, ALIGNSD));

  // ... + log P(sd)
  logL -= log(sdEll);

  return logL;
}


double
trace(gsl_matrix * X, int n) {
  double tr = 0;
  size_t i;
  for (i = 0; i < n; i++) {
    tr += gsl_matrix_get(X, i, i);
  }
  return tr;
}

double 
aLogPrior(gsl_matrix * A) {
  double EXP_R = 3;
  double DOF = 30;
  double logL = 0;
  

  gsl_permutation * perm = gsl_permutation_calloc(3);
  int s = 0;
  
  gsl_matrix *lu = gsl_matrix_alloc(3, 3);
  gsl_matrix *Ai = gsl_matrix_alloc(3, 3);
  gsl_matrix_memcpy(lu, A);
  gsl_linalg_LU_decomp(lu, perm, &s);
  gsl_matrix_memcpy(Ai, A);
  gsl_linalg_LU_invert(lu, perm, Ai);
  double tr = trace(Ai, 3);
  double lndet = gsl_linalg_LU_lndet(lu);

  logL += (3 - DOF/2) * log(EXP_R/DOF);
  logL += (DOF - 4)/2 * lndet;
  logL += (-1/2) * pow(EXP_R/DOF, 3) * tr;

  gsl_matrix_free(Ai);
  gsl_matrix_free(lu);

  return logL;
}

double
shiftLogPrior(gsl_vector * shift, double sdShift) {
  double logL = 0;
  size_t i;
  for (i = 0; i < 3; i++) {
    logL += log(gsl_ran_gaussian_pdf(gsl_vector_get(shift, i), sdShift));
  }
  
  logL -= log(sdShift);
}
