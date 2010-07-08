#ifndef MCMC_H
#define MCMC_H

#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

class State;
class Sampler;
class Chain;

/** The `MCMC` class is a master class for MCMC sampling. Given a
    state and a sampler it initializes and manages a number of chains
    which each updates on their state independently. Then it merges
    all of the sampled points to a final output matrix dependent on
    R_hat convergence of each chain. */

class MCMC {
public:

  size_t D, nsamples, nchains, nchainsamples;

  /// Initialize the MCMC with its sampler and state pair then run it
  /// to completion.
  MCMC (Sampler *samp, State *st, 
        const size_t nsamples, 
        const size_t nchains, 
        const size_t nadapt = 2000);

  /// Free the MCMC and each of its chains.
  ~MCMC ();

  /// Check convergence output
  gsl_vector *convergenceRhats ();

  /// Get the output matrix for chain N (ceil(nsamples/nchains) by D).
  gsl_matrix *samplesFromChain (size_t n);

  /// Get the merged output matrices (approx. nsamples by D).
  gsl_matrix *copyAllSamples ();

private:

  std::vector<Chain *> chainpool;
};

#endif
