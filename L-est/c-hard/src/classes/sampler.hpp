#ifndef SAMPLER_INTERFACE_H
#define SAMPLER_INTERFACE_H

#include <cstdio>
#include <classes/state.hpp>

class Sampler {
protected:
  view_f view;
  update_f update;
  init_f guess0;
  bool frozen;

public:

  /// Initialize the sampler with `view` and `update` functions to
  /// abstract handling the `State` object.
  Sampler (view_f view_in, update_f update_in, init_f guess0_in) {
    view = view_in;
    update = update_in;
    guess0 = guess0_in;
    frozen = false;
  }

  /// Default constructor for abstract and meta-Samplers.
  Sampler () {
    frozen = false;
  }

  /// Clone the sampler state.
  virtual Sampler *clone () = 0;

  /// Update a state in place by performing the appropriate jump.
  virtual void jump (State *) = 0;

  /// Copy the passed state and then update it.
  State *jumpClone (const State *s1) {
    State *s2 = s1->clone();

    jump(s2);
    return s2;
  }

  /// Request that the current sampler end its adaptive mode; note
  /// that this is only a request and until `isFrozen` confirms the
  /// freeze the sampler must still be considered adaptive. Steps to
  /// go into frozen mode include preventing further sampling method
  /// changes and jettisoning any memory needed for adaptive jumps.
  virtual void freeze () {
    forceFreeze();
  }

  /// Force the sampler to freeze wherever it is. This should be
  /// called only if `freeze` fails to work after some time.
  virtual void forceFreeze () {
    frozen = true;
  }

  /// Check to see if the current sampler is frozen.
  virtual bool isFrozen () {
    return frozen;
  }
  
  /// Initialize the state with whatever paramters are sampled over by
  /// this sampler. This is never responsible for initializing any
  /// other part of the sampler.
  virtual void initState (State *st) {
    update(guess0(), st);
  };

};

#endif
