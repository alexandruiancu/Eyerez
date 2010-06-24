#include "problem.h"
#include "data.h"

#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>

double 
logPostVectorized(const gsl_vector *x, void *param_obs);

double 
invLogPost(const gsl_vector *x, void *param_obs);

void
assembleA(const gsl_vector *Abits, gsl_matrix *A);

void
dfInvLogPost(const gsl_vector * x, void * params, gsl_vector * grad);

void
fdfInvLogPost(const gsl_vector * x, void * params, double * f, gsl_vector * g);

void
varianceEst(//double (*f)(const gsl_vector * x, void * params), 
            const gsl_vector * x, 
            void * params,
            gsl_matrix *sigma);

int
findMode(const obsmat obs, const gsl_vector *inits, gsl_vector *mode);
