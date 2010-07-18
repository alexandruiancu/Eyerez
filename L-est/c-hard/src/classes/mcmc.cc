#include "mcmc.hpp"

#include <classes/sampler.hpp>
#include <classes/state.hpp>
#include <classes/chain.hpp>

#include <cmath>

MCMC::MCMC (Sampler *samp, State *st,
            const size_t nsamples_in,
            const size_t nchains_in,
            const size_t nadapt) {
  D             = st->D;
  nsamples      = nsamples_in;
  nchains       = nchains_in;
  nchainsamples = (size_t)ceil(nsamples/nchains);

  // Initialize some things
  chainpool.reserve(5);
  for (size_t i = 0; i < nchains; i++) {
    samp->initState(st);
    Sampler *samp_c = samp->clone();
    State *st_c     = st->clone();
    Chain *c = new Chain(samp_c, st_c, nadapt, nchainsamples);
    chainpool.push_back(c);
  }

  // Cycle through each chain performing updates.
  const size_t BLOCKSIZE = 100;
  bool allDone;
  do {
    allDone = true; // Our initial hope...
    for (size_t i = 0; i < nchains; i++) {
      chainpool[i]->run(BLOCKSIZE);
      allDone &= chainpool[i]->isDone(); // is slowly dashed...
    }
  } while (!allDone); // until it's saved.

}

// TODO: Bus error!
gsl_matrix *
MCMC::copyAllSamples () {
  gsl_matrix *out = gsl_matrix_alloc(nsamples, D);

  size_t chaini, sampi;
  for (size_t i = 0; i < nsamples; i++) {
    chaini = i/nchainsamples; // What is the chainpool index?
    for (size_t j = 0; j < D; j++) {
      sampi = i%nchainsamples;
      gsl_matrix_set(out, i, j,
                     gsl_matrix_get(chainpool[chaini]->samples, 
                                    sampi, j));
    }
  }

  return out;
}
