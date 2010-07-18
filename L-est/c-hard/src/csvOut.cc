#include "csvOut.hpp"

#include <classes/observations.hpp>
#include <classes/problem.hpp>
#include <classes/metropolis.hpp>
#include <classes/sequentially.hpp>
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

  samp_vec samplers;
  samplers.reserve(6);
  samplers.push_back(new Metropolis (1,
                                     LEst::view_ell,
                                     LEst::update_ell,
                                     LEst::guess0_ell));

  samplers.push_back(new Metropolis (3,
                                     LEst::view_lambda,
                                     LEst::update_lambda,
                                     LEst::guess0_lambda));

  samplers.push_back(new Metropolis (3,
                                     LEst::view_rot,
                                     LEst::update_rot,
                                     LEst::guess0_rot));

  samplers.push_back(new Metropolis (3,
                                     LEst::view_shift,
                                     LEst::update_shift,
                                     LEst::guess0_shift));

  samplers.push_back(new Metropolis (3,
                                     LEst::view_cerr,
                                     LEst::update_cerr,
                                     LEst::guess0_cerr));

  samplers.push_back(new Metropolis (2,
                                     LEst::view_err,
                                     LEst::update_err,
                                     LEst::guess0_err));
  
  Sequentially *samp = new Sequentially(samplers);

  MCMC *mcmc = new MCMC (samp, st, 100, 3, 1);
  gsl_matrix *m = mcmc->copyAllSamples();

  // for (size_t i = 0; i < 100; i++) {
  //   gsl_vector_view v = gsl_matrix_row(m, i);
  //   outputRow(&v.vector, 15);
  // }

  delete obs;
}
