#ifndef VARIANCE_H
#define VARIANCE_H

#include <cmath>
#include <iostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

using namespace std;


class VarianceM {
public:
  size_t      N;
  gsl_matrix *self;             // The NxN variance matrix itself
  gsl_matrix *cholesky;         // The lower diagonal square root of the variance

  /// Initialize the matrix and then perform this->is(sigma).
  VarianceM(const size_t N_in, const double sigma);

  /// Initialize the matrix using the array of doubles diag as the
  /// diagonal entries.
  VarianceM(const size_t N_in, const double *diag);

  /// Initialize the matrix and then perform this->is(diag).
  VarianceM(const gsl_vector *diag);

  /// Initialize the matrix and then perform this->is(diag).
  VarianceM(const gsl_matrix *var);

  /// Copy the variance matrix var. Cholesky is not refactored and is
  /// therefore referenced exactly as in var.
  VarianceM(const VarianceM &var);

  /// Assign this variance to be equal to variance var. Cholesky is
  /// assumed consistent and thus not refactored.
  VarianceM& operator=(const VarianceM &var);

  ~VarianceM();

  friend ostream &
  operator<< (ostream &os, const VarianceM &var);

  /// Set the variance matrix as the parameter sigma; variance is
  /// assumed to be spherical with scale sigma.
  void is(const double sigma);

  /// Set the variance matrix so that the array of doubles diag is the
  /// main diagonal; assumes complete independence of variances;
  void is(const double *diag);

  /// Set the variance matrix as the parameter diag; assumes complete
  /// independence.
  void is(const gsl_vector *diag);

  /// Set the variance matrix as the parameter var. Var must be
  /// symmetric and positive semi-definite and thus
  /// cholesky-factorable.
  void is(const gsl_matrix *var);

  /// Compute the sample variance of a PxN matrix of P N-dimensional
  /// observations. If scale is specified the matrix is rescaled by
  /// that factor.
  void isSampleVariance(const gsl_matrix *obs, const double scale = 1.0);
  
  /// Scale the vector x as though it were a N(0,1) draw such that it
  /// follows mvN(0, self).
  void scaleNormal(gsl_vector *x);

protected:

  /// Factor self into the lower-diagonal cholesky matrix X = A At.
  void factor();
};

ostream &operator<< (ostream &os, const VarianceM &var);

#endif
