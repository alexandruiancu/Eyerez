#include <classes/chain.hpp>

#include <stdio.h>

#include <classes/state.hpp>
#include <classes/sampler.hpp>

Chain::Chain (State *st_in, const size_t nposterior) {
  st = st_in;
  initOutput(nposterior);
}

Chain::Chain (Sampler *samp_in, State *st_in, 
              const size_t nadaptive, 
              const size_t nposterior) {
  to_adapt = nadaptive;
  st = st_in;
  samp = samp_in;
  totalJumps = 0;
  samp->initState(st);
  initOutput(nposterior);
}

Chain::~Chain () {
  gsl_matrix_free(samples);
}

void
Chain::newOutput (const size_t nposterior) {
  gsl_matrix_free(samples);
  initOutput(nposterior);
}

void
Chain::initOutput (const size_t nposterior) {
  size_t D = st->D;
  n_curr = 0;
  samples = gsl_matrix_alloc(nposterior, D);
}

void
Chain::run (const size_t niter) {
  if (!st->isFullyInitialized()) {
    samp->initState(st);
  }
  for (size_t i = 0; i < niter; i++) {
    step();
    if (isDone()) { break; }
  }
}

void
Chain::step () {
  if (!st->isFullyInitialized()) {
    fprintf(stderr, "State is not fully initialized: failing...\n");
    throw;
  }

  if (!samp->isFrozen() && to_adapt > 0) {
    // We're in adaptive mode.
    fprintf(stderr, "a");
    samp->jump(st);
    to_adapt--;

  } else if (!samp->isFrozen()) {
    // We should start attempting to freeze it.
    fprintf(stderr, "f");
    samp->freeze();
    samp->jump(st);
    
  } else {
    // It's totally frozen. Now samples should be stored.
    fprintf(stderr, "+");
    samp->jump(st);
    storeState(st);
  }
}

bool
Chain::isDone () {
  return (n_curr >= samples->size1);
}

void 
Chain::storeState (State *x) {
  gsl_vector_view v = gsl_matrix_row(samples, n_curr);
  st->pickle(&v.vector);
  n_curr++;
}
