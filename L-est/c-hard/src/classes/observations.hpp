#ifndef OBSERVATIONS_H
#define OBSERVATIONS_H

#include <cmath>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <csv_parser/csv_parser.hpp>

enum CoordSystem { 
  PolarCoord,
  RectCoord
};

/// Wrapper class over gsl_matrix to perform some operations related
/// to dataset manipulations.
class Observations {
  double L; /// The internal parameter representing the distance
            /// between the laser geometry and the rotational axis of
            /// the eye.

  void
  init(const gsl_matrix *obs, const CoordSystem system_in);

public:
  gsl_matrix *data;
  size_t N, D;
  CoordSystem system;
  
  /// Create a new Observations class from an input gsl_matrix of
  /// observations in tall form (Nx3).
  Observations(const gsl_matrix *obs, const CoordSystem system_in = PolarCoord);

  /// Create a new Observations class from observations read from a
  /// csv file in tall form (Nx3).
  Observations(FILE *f, const CoordSystem system_in = PolarCoord);

  ~Observations();

  /// Set the L parameter for rsXXX readouts.
  void
  setL(const double L = 30.4);

  /// Writes the value of observation N to the passed vector.
  void
  writeObservation(const size_t n, gsl_vector *o);

  /// Allocates a new vector copied from a given row of observations
  /// already converted to real space.
  gsl_vector *
  allocAsObservation(const size_t n);
};

namespace Util {
  size_t countlines(FILE *f);
}

#endif
