#include "data.h"

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

double 
logLik(const obsmat os, 
       const double L, 
       const gsl_vector *c, 
       const gsl_matrix *A, 
       const double sdError);

double
computeDistance(const gsl_vector *ob, const gsl_vector *c, const gsl_matrix *A);

double 
ellLogPrior(double ell, double muEll, double sdEll);

double 
trace(gsl_matrix * X, size_t n);

double 
aLogPrior(gsl_matrix * A);

double
shiftLogPrior(gsl_vector * shift, double sdShift);
