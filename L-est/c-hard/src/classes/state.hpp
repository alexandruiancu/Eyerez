#ifndef STATE_INTERFACE_H
#define STATE_INTERFACE_H

#include <cstdio>
#include <gsl/gsl_vector.h>

/** Represents an overall state in a Markov chain. These objects are
    passed into samplers along with externally implemented functions
    to draw vectors representing subspaces of the motile parameters of
    the overall state and the complementary function which will take a
    similar vector and update the current state with the new
    parameters.

    If ST is a subclass of State, the externally implemented functions
    should have prototypes similar to:

    gsl_vector *view_i(const State *);
    int update_i(const gsl_vector *, State *);

    Individual states are responsible for keeping up with changes that
    occur via those external functions. For instance, the update_i
    functions should let the state know what components have changed
    so that later when the `repair` method is called the State object
    can update any internal states with respect to these new changes.

    This functionality allows a state to, for instance, recompute a
    particularly time intensive part of the likelihood function only
    when its dependent parameters are updated.

 */

class State {
public:

  /// Number of publically facing, traced parameters.
  size_t D;

  /// Clones the state without relying on built-in C++ copy
  /// machinery. This is for typechecking.
  virtual State *clone () const = 0;

  /// Test to see if the state has been fully initialized. Should be
  /// called before views to prevent unintialized values from slipping
  /// into the samplers.
  virtual bool isFullyInitialized () = 0;

  /// isValid implements any range truncation necessary to trim the
  /// proposal distribution to better match the expected posterior
  /// distribution. States considered invalid will be rejected from
  /// samplers without ticking the jump counter and therefore a false
  /// isValid is not equivalent to a -INFINITY likelihood.
  virtual bool isValid () {
    return true;
  }

  /// Computes the log likelihood of the current state. Necessary for
  /// Metropolis jumps.
  virtual double logLik () = 0;

  /// `repair` gives a state a chance to recalculate any internal
  /// variables after some outward facing ones have been updated during
  /// a jump. This MUST be called after any atomic update operation to
  /// ensure that the internal state is consistent.
  virtual void repair () = 0;

  /// Serialize the state into a gsl_vector (length D) so that it can
  /// be stored in a matrix tracking the state's progression.
  virtual void pickle (gsl_vector *out) = 0;

  /// Deserialize a pickled gsl_vector (length D) by fully
  /// initializing the current state with its values.
  virtual void unpickle (const gsl_vector *in) = 0;

};

// The two interface function types.
typedef gsl_vector *(*view_f)(const State *);
typedef int         (*update_f)(const gsl_vector *, State *);

/// A function that creates a random initial guess vector in order to
/// have the sampler initialize some part of the state. Output vector
/// must correspond to the input vector for the `update` function.
typedef gsl_vector *(*init_f)();


#endif
