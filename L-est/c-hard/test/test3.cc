#include <classes/problem.hpp>
#include <iostream>

int main () {
  Observations *os = new Observations(fopen("e1.csv", "r"), RectCoord);

  gsl_vector *lambda = gsl_vector_alloc(3);
  gsl_vector_set_all(lambda, 1.0/9.0);

  gsl_matrix *Q = gsl_matrix_alloc(3, 3);
  gsl_matrix_set_identity(Q);

  gsl_vector *c = gsl_vector_alloc(3);
  gsl_vector_set_zero(c);

  double ll = Problem::datasetLogLik(os, 30.4, lambda, Q, c, 0.1, 0.1);
  cerr << ll;
}
