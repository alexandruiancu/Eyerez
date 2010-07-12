#include "metropolis.hpp"

#include <cmath>
#include <iostream>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>

#include <classes/observations.hpp>
#include <classes/variance.hpp>
#include <classes/random.hpp>
#include <classes/problem.hpp>

const size_t Metropolis::countMax = 50;



Metropolis::Metropolis (const State *st,
                        view_f view, update_f update, init_f guess0,
                        const double sigma) 
  : Sampler (view, update, guess0) {

  // Find dimensionality from the view_f
  gsl_vector *vw = view(st);
  D = vw->size;
  gsl_vector_free(vw);

  init(D, view, update, guess0, sigma);
}

Metropolis::Metropolis (size_t D,
                        view_f view, update_f update, init_f guess0,
                        const double sigma)
  : Sampler (view, update, guess0) {

  init(D, view, update, guess0, sigma);

}

void
Metropolis::init (size_t D,
                  view_f view, update_f update, init_f guess0,
                  const double sigma) {

  // Initialize some support objects
  R = new Random ();
  var = new VarianceM (D, sigma);
  
  // Allocate memory
  history      = gsl_matrix_alloc(countMax, D);
  acceptReject = gsl_vector_alloc(countMax);

  // Initialize some variables
  count = 0;
  rescale = 1;

  // TODO: Estimate the variance from the curvature at initial state.
  //       Should replace the blanket, uncorrelated state from sigma.
}

Metropolis::~Metropolis() {
  delete R;
  delete var;
  gsl_matrix_free(history);
  gsl_vector_free(acceptReject);
}

Metropolis *
Metropolis::clone() {
  Metropolis *m = new Metropolis(D, view, update, guess0);
  m->var = var;
  return m;
}

void
Metropolis::jump(State *x) {
  State *xp = x->clone();
  updateWithProposal(xp);

  double alpha = xp->logLik() - x->logLik();
  double u = R->drawUniform();

  if (alpha >= 0 || u < exp(alpha)) {
    update(view(xp), x);
    accepted = true;
  } else {
    accepted = false;
  }
  delete xp;

  storeSample(x, accepted);
}

void
Metropolis::updateWithProposal(State *x) {
  gsl_vector *vw = view(x);
  gsl_vector *xp = gsl_vector_alloc(D);

  do {
    // Complete a draw
    for (size_t i = 0; i < D; i++) { 
      gsl_vector_set(xp, i, R->drawGaussian());
    }
    
    var->scaleNormal(xp);
    gsl_vector_add(xp, vw);

    // But reject until the state is valid
  } while (!x->isValid());

  gsl_vector_free(xp);
}

void 
Metropolis::storeSample(State *x, bool accepted) {
  if (!frozen) {
    // Copy state into history[count, :]
    gsl_vector *vw = view(x);
    for (size_t i = 0; i < D; i++) {
      // Add just a bit of noise to make covariance easier to
      // calculate.
      double perturb = R->drawGaussian(0.0, 0.001);
      gsl_matrix_set(history, count, i, perturb+gsl_vector_get(vw, i));
    }
    
    // Update the acceptance history
    gsl_vector_set(acceptReject, count, accepted);

    // Shift to the next row
    count++;
    
    checkForRescale();
  }
}

void
Metropolis::checkForRescale() {
  if (count >= countMax) {
    // Initialize a resample
    count = 0;
    
    var->isSampleVariance(history, rescale*2.4*sqrt(D));
    
    // Update the rescale factor
    double ar = 0;
    for (size_t i = 0; i < countMax; i++) {
      ar += gsl_vector_get(acceptReject, i);
    }
    ar = ar/countMax;

    // In 1-D, optimal acceptance ratio is 0.44.
    // In 5-D or larger it drops to 0.23.
    rescale = rescale*(1 + 4.2*(ar - 0.23));

    // Reset the history
    gsl_matrix_set_zero(history);

  }
}


// int
// varianceEst(double (*f)(const gsl_vector * x, Observations *obs), 
//             const gsl_vector * x, 
//             Observations *obs,
//             gsl_matrix *sigma) {

//   double toli = 0.001, tolj = 0.001;
//   size_t i, j, D = x->size;
//   gsl_vector *xcpy = gsl_vector_alloc(D);
//   gsl_vector_memcpy(xcpy, x);
  
//   // Compute the curvature
//   gsl_vector *edi, *edj;
//   edi = gsl_vector_alloc(D);
//   edj = gsl_vector_alloc(D);
//   double xi0, xj0;
//   double a, b, c, d, curv_ij;
//   for (i = 0; i < D; i++) {
//     for (j = i; j < D; j++) {

//       // Perturb the mode
//       if (i == j) {
//         // Simplified method
//         b = (*f)(x, obs);
//         c = b;

//         gsl_vector_set(xcpy, i, gsl_vector_get(x, i) + 2*toli);
//         a = (*f)(xcpy, obs);

//         gsl_vector_set(xcpy, i, gsl_vector_get(x, i) - 2*toli);
//         d = (*f)(xcpy, obs);

//         gsl_vector_memcpy(xcpy, x);

//       } else {
//         xi0 = gsl_vector_get(x, i);
//         xj0 = gsl_vector_get(x, j);

//         gsl_vector_set(xcpy, i, xi0 + toli);
//         gsl_vector_set(xcpy, j, xj0 + tolj);
//         a = (*f)(xcpy, obs);
        
//         gsl_vector_set(xcpy, i, xi0 - toli);
//         gsl_vector_set(xcpy, j, xj0 + tolj);
//         b = (*f)(xcpy, obs);
        
//         gsl_vector_set(xcpy, i, xi0 + toli);
//         gsl_vector_set(xcpy, j, xj0 - tolj);
//         c = (*f)(xcpy, obs);
        
//         gsl_vector_set(xcpy, i, xi0 - toli);
//         gsl_vector_set(xcpy, j, xj0 - tolj);
//         d = (*f)(xcpy, obs);

//         // Restore mode
//         gsl_vector_memcpy(xcpy, x);
//       }
      
//       curv_ij = (a - b - c + d)/(4*toli*tolj);

//       // Store in curvature
//       gsl_matrix_set(sigma, i, j, curv_ij);
//       gsl_matrix_set(sigma, j, i, curv_ij);
//     }
//   }
//   gsl_vector_free(xcpy);

//   // Compute variance as V = [-L''(x_mode)]^-1
//   gsl_matrix *lu = gsl_matrix_alloc(D, D);
//   gsl_permutation *perm = gsl_permutation_alloc(D);
//   int signum = 1;
//   gsl_matrix_memcpy(lu, sigma);
//   gsl_permutation_init(perm);

//   gsl_matrix_scale(lu, -1.0);
//   gsl_linalg_LU_decomp(lu, perm, &signum); // lu is now the LU-decomposition
  
//   int error = gsl_linalg_LU_invert(lu, perm, sigma); // sigma is now the variance matrix

//   return error;
//   // and we are done.

//   gsl_matrix_free(lu);
//   gsl_permutation_free(perm);
// }
