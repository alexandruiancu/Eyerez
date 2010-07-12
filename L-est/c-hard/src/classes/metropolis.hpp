#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <classes/state.hpp>
#include <classes/sampler.hpp>

class Observations;
class VarianceM;
class Random;

class Metropolis : public Sampler {
  // Resampling period
  const static size_t  countMax;

  // the log likelihood function
  double (*logL)(const gsl_vector *x, Observations *obs);

  Random       *R;              // A random number instance
  size_t        D;              // size of state vector
  size_t        count;          // counter till next covariance estimation
  Observations *obs;            // matrix of observations; params for logL
  bool          accepted;       // Was the last sample accepted or rejected?

  gsl_matrix   *history;        // last count vectors (count x N, tall style)
  gsl_vector   *acceptReject;   // History of acceptances/rejections

  void init (size_t D,
             view_f view, update_f update, init_f guess0,
             const double sigma = 1);

public:
  double        rescale;        // Value of scale expansion from sample variance
  VarianceM    *var;            // Current variance estimate
  gsl_vector   *state;          // The current position of the sampler

  /// Initialize the Metropolis sampler with the appropriate State
  /// manipulation functions.
  Metropolis (const State *st,
              view_f view, update_f update, init_f guess0,
              const double sigma = 1);

  /// Initialize the Metropolis sampler without having an example state.
  Metropolis (size_t D,
              view_f view, update_f update, init_f guess0,
              const double sigma = 1);

  ~Metropolis ();

  /// Clone the metropolis state
  Metropolis *clone ();

  /// Perform a Metropolis jump from the state x0.
  void jump (State *x0);

private:
  /// Update the passed State object as a new proposed state.
  void updateWithProposal(State *x);

  /// Store the viewed part of the new state and perform adaptive
  /// update; disabled when frozen.
  void storeSample(State *x, bool accepted);

  /// Check sampler state and rescale if necessary.
  void checkForRescale();

};

// /** Compute an estimate of the variance matrix as the inverse of the
//     negative curvature at the point, which should be a mode. This
//     assumes a quadratic log likelihood approximation. The curvature
//     is computed numerically using a 2D 5-point stencil.
// */  
// int varianceEst(double (*f)(const gsl_vector * x, const gsl_matrix *obs),
//                 const gsl_vector * x, 
//                 Observations *obs,
//                 gsl_matrix *sigma);

#endif
