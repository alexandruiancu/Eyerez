#ifndef CHAIN_H
#define CHAIN_H

#include <gsl/gsl_matrix.h>

// Forward declarations.
class Sampler;
class State;

class Chain {
public:

  /// The chain's ouput matrix of posterior samples.
  gsl_matrix *samples;

  /// Initialize a Chain and clear enough memory for the output
  /// matrix.
  Chain (State *st, const size_t nposterior);

  /// Create a new chain object which encapsulates a sampler and a
  /// state. During initialization the state is initialized as well by
  /// using the Sampler. Also, the history matrix should be
  /// initialized to prove that the memory exists. No jumps are
  /// performed yet.
  Chain (Sampler *samp, State *st, 
         const size_t nadaptive, 
         const size_t nposterior);

  /// Clean up the chain. It is not responsible for the sampler or the
  /// state it was initialized with. This will be primarily used for
  /// cleaning the history matrix.
  ~Chain ();

  /// Run the sampler for niter steps. These steps may be adaptive or
  /// posterior sampling steps. If the chain completes before jumping
  /// the full niter steps it will do so silently and should be
  /// checked at a later time using `isDone`.
  void run (const size_t niter);

  /// Returns true if the sampler still needs to be `run` more before
  /// it completely fills its history matrix.
  bool isDone ();

  /// Link in a new output matrix. The sampler will forget and free
  /// its previous `nposterior`-sized matrix, create a new matrix, and
  /// seek to fill it instead.
  void newOutput (const size_t nposterior);

private:

  State *st;
  Sampler *samp;
  size_t to_adapt, N_total, n_curr;

  /// Initialize the output matrix and the counters to maintain
  /// position within it.
  void initOutput (const size_t nposterior);

  /// Stores a state in the current position of the output samples
  /// matrix and updates its index.
  void storeState (State *x);

  /// The total number of jumps performed so far.
  size_t totalJumps;

  /// Take a single step, adaptive or otherwise as appropriate.
  void step ();

};

#endif
