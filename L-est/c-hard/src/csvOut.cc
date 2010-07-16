#include "csvOut.hpp"

#include <classes/observations.hpp>
#include <classes/problem.hpp>
#include <classes/metropolis.hpp>
#include <classes/mcmc.hpp>

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
  fclose(fh);

  LEst::Estimate *st = new LEst::Estimate (obs);

  Metropolis *samp = new Metropolis (st,
                                     LEst::view_ell,
                                     LEst::update_ell,
                                     LEst::guess0_ell);

  MCMC *mcmc = new MCMC (samp, st, 4000, 3);

  delete obs;
}
