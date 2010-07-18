#ifndef SEQUENTIALLY_H
#define SEQUENTIALLY_H

#include <vector>

#include <classes/sampler.hpp>

/** Updates the state by running through a collection of samplers
    sequentially. */

typedef std::vector <Sampler *> samp_vec;

class Sequentially : public Sampler {

public:
  /// The vector of constituent samplers.
  samp_vec samplers;
  
  /// The number of constituent samplers.
  samp_vec::size_type S;

  /// Initialize the Sequentially sampler with a std::vector of
  /// Samplers. The Sequentially sampler takes ownership of each of
  /// its constituents and will free them upon deletion.
  Sequentially (samp_vec samplers_in);
  
  /// Frees each constituent sampler.
  ~Sequentially ();
  
  /// Clone each constituent sampler to form a new Sequentially
  /// sampler.
  Sequentially *clone ();

  /// Running through the samplers vector, jump the state in each
  /// sampler then return.
  void jump (State *);

  /// Call freeze on every constituent sampler.
  void freeze ();

  /// The union of the `isFrozen` states of each constituent sampler.
  bool isFrozen ();

  /// Initialize the state by using the each constituent sampler
  /// sequentially. Overridden since Sequentially does not have its
  /// own `view`, `update`, or `guess0` functions.
  void initState (State *st);
};

#endif
