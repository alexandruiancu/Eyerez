#include <classes/metropolis.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <iostream>

using namespace std;

gsl_vector *
allocRandomInit() {
  Random *R = new Random();

  gsl_vector *x0 = gsl_vector_alloc(2);

  gsl_vector_set(x0, 0, R->drawUniform(0, 1));
  gsl_vector_set(x0, 1, R->drawUniform(0, 1));

  delete R;
  return x0;
}

namespace simpleEx {

  struct StateOf {
    double p1;
    double p2;
    
    StateOf(const gsl_vector *x) {
      p1 = gsl_vector_get(x, 0);
      p2 = gsl_vector_get(x, 1);
    }

    gsl_vector *
    allocVectorized() {
      gsl_vector *out = gsl_vector_alloc(2);
      gsl_vector_set(out, 0, p1);
      gsl_vector_set(out, 1, p2);

      return out;
    }

    void
    writeIntoVec(gsl_vector *out) {
      gsl_vector_set(out, 0, p1);
      gsl_vector_set(out, 1, p2);
    }
  };
}


double
logPost(const gsl_vector *x0, const gsl_matrix *obs) {
  simpleEx::StateOf st = simpleEx::StateOf(x0);
  
  double M1, N1, M2, N2, logP;
  M1 = 100;
  N1 = 30;
  M2 = 230;
  N2 = 45;

  if (st.p1 < 0 || st.p2 < 0 || st.p1 > 1 || st.p2 > 1) {
    logP = -INFINITY;
  } else {
    logP = 
      N1 * log(st.p1) +
      (M1-N1) * log(1 - st.p1) +
      N2 * log(st.p2) +
      (M2-N2) * log(1 - st.p2);
  }

  return logP;
}

int main () {

  gsl_vector *x0 = allocRandomInit();
  gsl_matrix *obs = gsl_matrix_alloc(1, 1);
  Metropolis *m1 = new Metropolis(logPost, x0, obs, 1);

  for (size_t i = 0; i < 2000; i++) { m1->jump(); }
  m1->freeze();

  for (size_t i = 0; i < 4000; i++) {
    m1->jump();
    for (size_t i = 0; i < 2; i++) {
      if (i == 12) {
        fprintf(stdout, "%f  ", gsl_vector_get(m1->state, i));
      } else {
        fprintf(stdout, "%f, ", gsl_vector_get(m1->state, i));
      }
    }
    fprintf(stdout, "\n");
  }

  delete m1;
  gsl_matrix_free(obs);
  gsl_vector_free(x0);
}
