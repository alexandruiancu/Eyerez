#include "data.h"

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

double 
logPost(const obsmat os, 
        const double ell,
        const double sdEll,
        const gsl_vector *c,
        const double sdShift,
        const gsl_matrix *A, 
        const double sdError);

double 
logLik(const obsmat os, 
       const double ell, 
       const gsl_vector *c, 
       const gsl_matrix *A, 
       const double sdError);

double
computeDistance(const gsl_vector *ob, const gsl_vector *c, const gsl_matrix *A);

double 
ellLogPrior(const double ell, const double sdEll);

double 
trace(const gsl_matrix * X);

double 
aLogPrior(const gsl_matrix * A);

double
shiftLogPrior(const gsl_vector * shift, const double sdShift);
