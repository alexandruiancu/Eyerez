#ifndef PROBLEM_H
#define PROBLEM_H

#include <cmath>
#include <iostream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

#include <classes/observations.hpp>
#include <classes/random.hpp>

namespace Problem {

  /** Defines a problem, a function of the type 
      double (*)(const gsl_vector *x, Observations *obs_in)
      which is the logLikelihood of the problem space.

      For this problem the units are in um.
  */

  double
  logPost(const gsl_vector *x, Observations *obs_in);

  class State {
  public:

    Observations *obs;

    double muL;
    double sdL;
    double L;
    double A;
    double B;
    double C;
    double rot1;
    double rot2;
    double rot3;
    double cx;
    double cy;
    double cz;
    double sdCx;
    double sdCy;
    double sdCz;
    double sdR;
    double sdZ;

    /// Create a random state without initializing on Observations.
    State();

    /// Create a state representation from the vector x.
    State(const gsl_vector *x, Observations *obs_in);

    /// Create a vector representation of the state; useful for
    /// manually generating state vectors.
    gsl_vector *allocVectorized();

    /// Check truncations and state whether this state has a positive
    /// likelihood of existence.
    int isValid();

    /// Calculate the log likelihood of this state.
    double logLik();

    double priors();
  };
  
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
}

#endif
