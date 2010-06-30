#include <classes/metropolis.hpp>
#include <classes/problem.hpp>

#include <iostream>

int main () {

  Observations *obs = new Observations(fopen("e-rect.csv", "r"), RectCoord);

  Problem::State s;
  s.obs = obs;
  gsl_vector *init = s.allocVectorized();

  Metropolis *m1 = new Metropolis (Problem::logPost, init, obs, 0.01);

  // Burn-in stage.
  double temp;
  for (size_t i = 0; i < 2000; i++) {
    m1->jump();
    if (rand()%100 < 10) {
      cerr << "Rescale is: " << m1->rescale << "\n";
    }
  }

  m1->freeze();
  for (size_t i = 0; i < 4000; i++) {
    m1->jump();
    for (size_t i = 0; i < 17; i++) {
      if (i == 16) {
        fprintf(stdout, "%f  ", gsl_vector_get(m1->state, i));
      } else {
        fprintf(stdout, "%f, ", gsl_vector_get(m1->state, i));
      }
    }
    fprintf(stdout, "\n");
  }

  delete m1;
  gsl_vector_free(init);
}

