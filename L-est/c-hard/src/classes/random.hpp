#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include <gsl/gsl_randist.h>

using namespace std;


class Random {
  const gsl_rng_type *rngT;
  gsl_rng *rng;

public:

  // Create a new random generator.
  Random();

  // These two aren't implemented yet.
  Random(const Random &r);
  Random &operator= (const Random &r);

  // Free it.
  ~Random();

  // Sample a random gaussian variate N(mu, sigma).
  double drawGaussian(double mu = 0.0, double sigma = 1.0);

  // Sample a random uniform variate U(min = 0, max = 1). Should be
  // constrained st. u ~ U, u in [min, max).
  double drawUniform(double min = 0.0, double max = 1.0);

  // Sample a random gamma variance.
  double drawGamma(double a = 1.0, double b = 1.0);
  
  // Sample a random integer between [0, maxIndex-1] including the
  // entire size of an unsigned long integer as its range.
  unsigned long int
  drawUniformIndex(unsigned long int maxIndex);
};

#endif
