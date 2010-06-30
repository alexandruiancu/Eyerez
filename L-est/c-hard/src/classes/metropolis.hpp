#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <classes/observations.hpp>
#include <classes/variance.hpp>
#include <classes/problem.hpp>
#include <classes/random.hpp>

#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

class Metropolis {
  // Resampling period
  const static size_t  countMax;

  // the log likelihood function
  double (*logL)(const gsl_vector *x, Observations *obs);

  Random       *R;              // A random number instance
  size_t        N;              // size of state vector
  size_t        count;          // counter till next covariance estimation
  Observations *obs;            // matrix of observations; params for logL
  int           accepted;       // Was the last sample accepted or rejected?

  int           burning;        // boolean: is in the burnin mode?
  gsl_matrix   *history;        // last count vectors (count x N, tall style)
  gsl_vector   *acceptReject;   // History of acceptances/rejections

public:
  double        rescale;        // Value of scale expansion from sample variance
  VarianceM    *var;            // Current variance estimate
  gsl_vector   *state;          // The current position of the sampler

  Metropolis(double (*logL)(const gsl_vector *x, Observations *obs), 
             const gsl_vector *state_init,
             Observations *obs_in,
             const double sigma = 0.1);

  ~Metropolis();

  ///  Performs a metropolis jump updating the parameter vector x
  ///  according to the metropolis-hastings algorithm.
  void jump(double temp = 1.0);

  /// Ends the burnin stage and proceeds to full sampling.
  void freeze();

private:
  /// Sample a new proposal point from mvNorm(state, variance).
  gsl_vector *allocProposal();

  void storeSample();

  /// Check sampler state and rescale if necessary.
  void checkForRescale();
};


/** Compute an estimate of the variance matrix as the inverse of the
    negative curvature at the point, which should be a mode. This
    assumes a quadratic log likelihood approximation. The curvature
    is computed numerically using a 2D 5-point stencil.
*/  
int varianceEst(double (*f)(const gsl_vector * x, const gsl_matrix *obs),
                const gsl_vector * x, 
                Observations *obs,
                gsl_matrix *sigma);

#endif
