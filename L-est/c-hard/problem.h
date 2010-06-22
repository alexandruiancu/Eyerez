#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

double 
distanceLogLik(gsl_matrix * data, double L, gsl_vector c, gsl_matrix * A);

double 
ellLogPrior(double ell, double muEll, double sdEll);

double 
trace(gsl_matrix * X, int n);

double 
aLogPrior(gsl_matrix * A);

double
shiftLogPrior(gsl_vector * shift, double sdShift);
