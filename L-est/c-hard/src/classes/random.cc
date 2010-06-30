#include "random.hpp"

#include <iostream>
#include <gsl/gsl_randist.h>


Random::Random() {
  gsl_rng_env_setup();
  rngT = gsl_rng_default;
  rng = gsl_rng_alloc(rngT);
}

Random::Random(const Random &r) {
  cerr << "Random number generator copying not implemented yet.";
  throw;
}

Random &
Random::operator= (const Random &r) {
  cerr << "Random number generator assignment not implemented yet.";
  throw;
}

Random::~Random() {
  gsl_rng_free(rng);
}


double 
Random::drawGaussian(double mu, double sigma) {
  double x0 = gsl_ran_gaussian(rng, sigma);
  return (mu + x0);
}

double 
Random::drawUniform(double min, double max) {
  return gsl_ran_flat(rng, min, max);
}

double
Random::drawGamma(double a, double b) {
  return gsl_ran_gamma(rng, a, b);
}
