#include "metropolis.hpp"

const size_t Metropolis::countMax = 50;

Metropolis::Metropolis(double (*logL_in)(const gsl_vector *x, Observations *obs_in), 
                       const gsl_vector *state_init,
                       Observations *obs_in,
                       const double sigma) {

  // Find dimensionality from initial state
  N = state_init->size;

  R = new Random();
  var = new VarianceM(N, sigma);
  
  // Allocate memory
  obs          = obs_in;
  state        = gsl_vector_alloc(N);
  history      = gsl_matrix_alloc(countMax, N);
  acceptReject = gsl_vector_alloc(countMax);

  // Initialize some variables
  count = 0;
  burning = 1;
  rescale = 1;

  logL = logL_in;

  gsl_vector_memcpy(state, state_init);

  // Estimate the variance from the curvature at initial state.
  //
  //
  // TODO.
  //
  //
  //
}

Metropolis::~Metropolis() {
  delete R;
  delete var;
  gsl_vector_free(state);
  gsl_matrix_free(history);
  gsl_vector_free(acceptReject);
}

void 
Metropolis::jump(double temp) {
  gsl_vector *xp = allocProposal();
  accepted = 0;

  double alpha = (*logL)(xp, obs) - (*logL)(state, obs);
  double u = R->drawUniform();
  
  if (alpha >= 0 || u < exp(alpha)) {
    // accept proposal
    gsl_vector_memcpy(state, xp);
    accepted = 1;
  }

  storeSample();
  gsl_vector_free(xp);
}

void
Metropolis::freeze() {
  burning = 0;
}

gsl_vector *
Metropolis::allocProposal() {
  gsl_vector *xp = gsl_vector_alloc(N);

  // Draw N uncorrelated N(0, 1)
  for (size_t i = 0; i < N; i++) { 
    gsl_vector_set(xp, i, R->drawGaussian());
  }

  var->scaleNormal(xp);
  gsl_vector_add(xp, state);

  return xp;
}

void 
Metropolis::storeSample() {
  if (burning) {
    // Copy state into history[count, :]
    for (size_t i = 0; i < N; i++) {
      // Add just a bit of noise to make covariance easier to
      // calculate.
      double perturb = R->drawGaussian(0.0, 0.001);
      gsl_matrix_set(history, count, i, perturb+gsl_vector_get(state, i));
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
    
    var->isSampleVariance(history, rescale*2.4*sqrt(N));
    
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


int
varianceEst(double (*f)(const gsl_vector * x, Observations *obs), 
            const gsl_vector * x, 
            Observations *obs,
            gsl_matrix *sigma) {

  double toli = 0.001, tolj = 0.001;
  size_t i, j, N = x->size;
  gsl_vector *xcpy = gsl_vector_alloc(N);
  gsl_vector_memcpy(xcpy, x);
  
  // Compute the curvature
  gsl_vector *edi, *edj;
  edi = gsl_vector_alloc(N);
  edj = gsl_vector_alloc(N);
  double xi0, xj0;
  double a, b, c, d, curv_ij;
  for (i = 0; i < N; i++) {
    for (j = i; j < N; j++) {

      // Perturb the mode
      if (i == j) {
        // Simplified method
        b = (*f)(x, obs);
        c = b;

        gsl_vector_set(xcpy, i, gsl_vector_get(x, i) + 2*toli);
        a = (*f)(xcpy, obs);

        gsl_vector_set(xcpy, i, gsl_vector_get(x, i) - 2*toli);
        d = (*f)(xcpy, obs);

        gsl_vector_memcpy(xcpy, x);

      } else {
        xi0 = gsl_vector_get(x, i);
        xj0 = gsl_vector_get(x, j);

        gsl_vector_set(xcpy, i, xi0 + toli);
        gsl_vector_set(xcpy, j, xj0 + tolj);
        a = (*f)(xcpy, obs);
        
        gsl_vector_set(xcpy, i, xi0 - toli);
        gsl_vector_set(xcpy, j, xj0 + tolj);
        b = (*f)(xcpy, obs);
        
        gsl_vector_set(xcpy, i, xi0 + toli);
        gsl_vector_set(xcpy, j, xj0 - tolj);
        c = (*f)(xcpy, obs);
        
        gsl_vector_set(xcpy, i, xi0 - toli);
        gsl_vector_set(xcpy, j, xj0 - tolj);
        d = (*f)(xcpy, obs);

        // Restore mode
        gsl_vector_memcpy(xcpy, x);
      }
      
      curv_ij = (a - b - c + d)/(4*toli*tolj);

      // Store in curvature
      gsl_matrix_set(sigma, i, j, curv_ij);
      gsl_matrix_set(sigma, j, i, curv_ij);
    }
  }
  gsl_vector_free(xcpy);

  // Compute variance as V = [-L''(x_mode)]^-1
  gsl_matrix *lu = gsl_matrix_alloc(N, N);
  gsl_permutation *perm = gsl_permutation_alloc(N);
  int signum = 1;
  gsl_matrix_memcpy(lu, sigma);
  gsl_permutation_init(perm);

  gsl_matrix_scale(lu, -1.0);
  gsl_linalg_LU_decomp(lu, perm, &signum); // lu is now the LU-decomposition
  
  int error = gsl_linalg_LU_invert(lu, perm, sigma); // sigma is now the variance matrix

  return error;
  // and we are done.

  gsl_matrix_free(lu);
  gsl_permutation_free(perm);
}
