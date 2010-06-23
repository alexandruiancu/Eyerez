#include "data.h"

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <csv_parser/csv_parser.hpp>

/*
  A data matrix of observations, os, is a 3xN matrix with each column
  vector representing a single observation.
 */

obsmat 
allocObsCSV(FILE *f) {
  int n = countlines(f);
  obsmat os = gsl_matrix_alloc(3, n);
  
  csv_parser file_parser;

  file_parser.init(f);
  file_parser.set_field_term_char(',');
  file_parser.set_line_term_char('\n');

  int i = 0;

  while (file_parser.has_more_rows()) {
    csv_row row = file_parser.get_row();
    
    for (unsigned int j = 0; j < row.size(); j++) {
      gsl_matrix_set(os, i, j, atof(row[j].c_str()));
    }
  }

  return os;
}

int
countlines(FILE *f) {
  /* Count the number of lines in file f and then return the handler
     at its initial position.
  */
  int lines = 0;
  int c;
  
  do {
    c = fgetc(f);
    if (c == '\n') {
      lines++;
    }
  } while (c != EOF);

  rewind(f);
  return lines;
}

void
toRealSpace(obsmat os, double L, int sign) {
  size_t N = os->size2;
  size_t i;
  double hole;

  for (i = 0; i < N; i++) {
    hole = gsl_matrix_get(os, 0, i);
    gsl_matrix_set(os, 0, i, L - hole);
  }
}

void
toLaserSpace(obsmat os, double L, int sign) {
  size_t N = os->size2;
  size_t i;
  double hole;

  for (i = 0; i < N; i++) {
    hole = gsl_matrix_get(os, 0, i);
    gsl_matrix_set(os, 0, i, L - hole);
  }  
}
