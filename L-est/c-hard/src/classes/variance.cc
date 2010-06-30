#include "variance.hpp"

#include <cmath>
#include <iostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

VarianceM::VarianceM(const size_t N_in, const double sigma) {
  N = N_in;
  self = gsl_matrix_alloc(N, N);
  cholesky = gsl_matrix_alloc(N, N);
  is(sigma);
}

VarianceM::VarianceM(const size_t N_in, const double *diag) {
  N = N_in;
  self = gsl_matrix_alloc(N, N);
  cholesky = gsl_matrix_alloc(N, N);
  is(diag);
}

VarianceM::VarianceM(const gsl_vector *diag) {
  N = diag->size;
  self = gsl_matrix_alloc(N, N);
  cholesky = gsl_matrix_alloc(N, N);
  is(diag);
}

VarianceM::VarianceM(const gsl_matrix *var) {
  N = var->size1;
  self = gsl_matrix_alloc(N, N);
  cholesky = gsl_matrix_alloc(N, N);
  is(var);
}

VarianceM::VarianceM(const VarianceM &var) {
  N = var.N;
  gsl_matrix_memcpy(self, var.self);
  gsl_matrix_memcpy(cholesky, var.cholesky);
}

VarianceM& 
VarianceM::operator=(const VarianceM &var) {
  N = var.N;
  gsl_matrix_memcpy(self, var.self);
  gsl_matrix_memcpy(cholesky, var.cholesky);
  return *this;
}

VarianceM::~VarianceM() {
  gsl_matrix_free(self);
  gsl_matrix_free(cholesky);
}

void 
VarianceM::is(const double sigma) {
  gsl_matrix_set_zero(self);
  for (size_t i = 0; i < N; i++) {
    gsl_matrix_set(self, i, i, sigma);
  }
  factor();
}

void 
VarianceM::is(const double *diag) {
  gsl_matrix_set_zero(self);
  for (size_t i = 0; i < N; i++) {
    gsl_matrix_set(self, i, i, diag[i]);
  }
  factor();
}

void 
VarianceM::is(const gsl_vector *diag) {
  gsl_matrix_set_zero(self);
  for (size_t i = 0; i < N; i++) {
    gsl_matrix_set(self, i, i, gsl_vector_get(diag, i));
  }
  factor();
}

void 
VarianceM::is(const gsl_matrix *var) {
  gsl_matrix_memcpy(self, var);
  factor();
}

void 
VarianceM::isSampleVariance(const gsl_matrix *obs, const double scale) {
  size_t P = obs->size1;
  size_t Nobs = obs->size2;
  gsl_matrix *obscpy = gsl_matrix_alloc(P, N);
  if (Nobs != N) throw;
  
  gsl_matrix_memcpy(obscpy, obs);
  
  // Compute the center
  double acc;
  gsl_vector *center = gsl_vector_alloc(N);
  for (size_t j = 0; j < N; j++) {
    acc = 0;
    for (size_t i = 0; i < P; i++) {
      acc += gsl_matrix_get(obs, i, j);
    }
    gsl_vector_set(center, j, acc/P);
  }

  // Add up the variances
  gsl_matrix *sampleVar = gsl_matrix_alloc(N, N);
  gsl_matrix_set_zero(sampleVar);

  gsl_vector_view view;
  for (size_t q = 0; q < P; q++) {
    view = gsl_matrix_row(obscpy, q);
    
    // Center the obs matrix
    gsl_vector_sub(&view.vector, center);
    
    gsl_blas_dger(1.0/(P - 1), &view.vector, &view.vector, sampleVar);
  }

  gsl_matrix_scale(sampleVar, scale);
  is(sampleVar);
  
  gsl_matrix_free(sampleVar);
  gsl_vector_free(center);
  gsl_matrix_free(obscpy);
}

void 
VarianceM::scaleNormal(gsl_vector *x) {
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, cholesky, x);
}

void 
VarianceM::factor() {
  gsl_matrix *temp = gsl_matrix_alloc(N, N);
  gsl_matrix_memcpy(temp, self);
  
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  int error = gsl_linalg_cholesky_decomp(temp);
  gsl_set_error_handler(old_handler);
  
  if (error != GSL_EDOM) {
    // Remove the upper triangular part
    for (size_t i = 0; i < N; i++) {
      for (size_t j = i; j < N; j++) {
        if (i != j) { gsl_matrix_set(temp, i, j, 0.0); }
      }
    }      
    gsl_matrix_memcpy(cholesky, temp);
  } else { 
    cerr << "Failed to factor sample covariance!\n";
    throw error;
  }    
  
  gsl_matrix_free(temp);
}


ostream &
operator<< (ostream &os, const VarianceM &var) {
  size_t N = var.N;
  
  os << "Variance is: \n";
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      os.precision(4);
      os.width(12);
      os << gsl_matrix_get(var.self, i, j) << " ";
    }
    os << "\n";
  }
  os << "\n";
  
  os << "Cholesky is: \n";
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      os.precision(4);
      os.width(12);
      os << gsl_matrix_get(var.cholesky, i, j) << " ";
    }
    os << "\n";
  }  
  
  return os;
}
