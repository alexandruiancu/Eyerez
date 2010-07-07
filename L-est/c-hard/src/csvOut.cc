#include "csvOut.hpp"

void outputRow(gsl_vector *m, size_t N) {
    for (size_t i = 0; i < N; i++) {      
      fprintf(stdout, "%f", gsl_vector_get(m, i));
      if (i == N-1) {
        fprintf(stdout, " ");
      } else {
        fprintf(stdout, ", ");
      }
    }
    fprintf(stdout, "\n");
}

int main (int argc, char *argv[]) {

  FILE * fh = fopen(argv[1], "r");

  if (!fh) {
    fprintf(stderr, "Could not find CSV file!\n");
    throw;
  }

  fprintf(stderr, "Grabbed file, starting sims... \n");

  Observations *obs = new Observations(fh, PolarCoord);
  
  Problem::State s;
  s.obs = obs;
  gsl_vector *init = s.allocVectorized();
  
  Metropolis *m1 = new Metropolis (Problem::logPost, init, obs, 0.01);
  
  
  for (size_t i = 0; i < 2000; i++) { 
    m1->jump(); 
    fprintf(stderr, "+");
  }

  fprintf(stderr, "Burnt-in, starting output...\n");
  
  m1->freeze();
  for (size_t i = 0; i < 4000; i++) { 
    m1->jump();   
    outputRow(m1->state, 17);
    fprintf(stderr, ".");
  }
  
  delete m1;
  gsl_vector_free(init);    
  delete obs;
}
