#ifndef PROBLEM_H
#define PROBLEM_H

#include <cmath>
#include <iostream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <classes/state.hpp>

// Forward declarations
class Observations;
class Random;

namespace MyUtilities {
  class parameter {

    double self;
    bool initialized;
    bool edited;
    
  public:

    parameter ();

    /// Getter and setter.
    double operator() () const;
    void operator() (double new_value);

    /// Determine if the slot has been initialized.
    bool isInitialized ();

    /// Determine if the slot has been changed since the last repair.
    bool isEdited ();

    /// Clean the value after a repair.
    void clean ();

  };
}

namespace LEst {

  class Estimate : public State {

    const static size_t NSAMPLES = 100;
    
    Observations *observations;

    // Priors
    double logLik, priorQ, priorLambda, priorErr;
    double priorC, priorCErr, priorL;

    double computeLogLik();
    double computePriorQ();
    double computePriorLambda();
    double computePriorErr();
    double computePriorC();
    double computePriorCErr();
    double computePriorL();

    void repair ();

    

    /// Compute the value t, used in finding the shortest line between the
    /// ellipsoid surface and an arbitrary point. The ellipsoid surface is
    /// defined to be rectangular and zero-centered with major axis
    /// lengths of (a, b, c).
    double 
    findT(const double a, const double b, const double c,
          const double u, const double v, const double w);
    
    /// Compute the difference vector (y-x) for any point y
    /// where x is the nearest point on the rectangular, zero-centered
    /// ellipse with major axis lengths (a, b, c). This vector is stored
    /// in the passed output parameter out.
    void
    computeDelta(const double a, const double b, const double c,
                 const double u, const double v, const double w,
                 gsl_vector *out);
    
    
    /// Calculate log error likelihood for a given measurement against 
    /// a fixed ellipsoid.
    double
    ellipsoidLogPdf(const double a, const double b, const double c,
                    const double u, const double v, const double w,
                    const double sdR, 
                    const double sdZ);
    
    /// observationLogLik without directly allocating inner work
    /// vectors.
    double
    _observationLogLik(const gsl_vector *ob,
                       const gsl_vector *lambda,
                       const gsl_matrix *Q,
                       const gsl_vector *c,
                       const double sdR,
                       const double sdZ,
                       gsl_vector *uvw,
                       gsl_vector *ytmp);
    
    /// Compute the log likelihood of the independent observation matrix
    /// given that they are measurements of a given ellipsoid with normal,
    /// independent errors. The ellipsoid is defined by the decomposition
    /// of the matrix A where ellipsoid points are x st. t(x) A x - 1 == 0
    /// and the vector c which pinpoints its center.
    /// The positive-definite matrix A is eigendecomposed to A = Q L Qt
    /// where Q are the eigenvectors corresponding to the eigenvalues
    /// diag(L). In this way, Qt rotates vectors into the eigenspace of
    /// the ellipsoid (ie. the space where the major axes of the ellipsoid
    /// are aligned with the basis).
    double
    observationLogLik(const gsl_vector *ob,
                      const gsl_vector *lambda,
                      const gsl_matrix *Q,
                      const gsl_vector *c,
                      const double sdR,
                      const double sdZ);
    
    /// Computes the log likelihood of the complete dataset of
    /// observations according to the model that they are
    /// rectangular-space points which have been transformed in
    /// polar-space by [r <- L-r] for a known L; these (x, y, z) points
    /// are assumed to fit a known ellipsoid via observationLogLik.
    double
    datasetLogLik(Observations *obs,
                  const double L,
                  const gsl_vector *lambda,
                  const gsl_matrix *Q,
                  const gsl_vector *c,
                  const double sdR,
                  const double sdZ);

  public:

    /// Public facing parameter interface
    MyUtilities::parameter L;
    MyUtilities::parameter A, B, C;
    MyUtilities::parameter rot1, rot2, rot3;
    MyUtilities::parameter cx, cy, cz;
    MyUtilities::parameter sdCx, sdCy, sdCz;
    MyUtilities::parameter sdR, sdZ;    

    /// Initialize a state with its observations. The state will
    /// maintain the object and delete it when necessary.
    Estimate (Observations *obs);

    /// Cleans up the `Observations` object.
    ~Estimate ();

    /// Clone this estimate
    Estimate *clone () const;

    /// Determines if every slot has been initialized by samplers.
    bool isFullyInitialized ();

    /// Returns the result of boundchecking for sampler truncations.
    bool isValid ();

    double logPost ();

    void pickle (gsl_vector *out) const;

    void unpickle (const gsl_vector *in);
  };

  gsl_vector *view_ell(const State *);
  gsl_vector *view_lambda(const State *);
  gsl_vector *view_rot(const State *);
  gsl_vector *view_shift(const State *);
  gsl_vector *view_cerr(const State *);
  gsl_vector *view_err(const State *);

  void update_ell(const gsl_vector *, State *);
  void update_lambda(const gsl_vector *, State *);
  void update_rot(const gsl_vector *, State *);
  void update_shift(const gsl_vector *, State *);
  void update_cerr(const gsl_vector *, State *);
  void update_err(const gsl_vector *, State *);

  gsl_vector *guess0_ell();
  gsl_vector *guess0_lambda();
  gsl_vector *guess0_rot();
  gsl_vector *guess0_shift();
  gsl_vector *guess0_cerr();
  gsl_vector *guess0_err();

  // view_f   view_ell, view_lambda, view_rot, view_shift, view_cerr, view_err;
  // update_f update_ell, update_lambda, update_rot, update_shift, update_cerr, update_err;
  // init_f   guess0_ell, guess0_lambda, guess0_rot, guess0_shift, guess0_cerr, guess0_err;
}

#endif
