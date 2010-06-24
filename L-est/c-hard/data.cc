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
      gsl_matrix_set(os, j, i, atof(row[j].c_str()));
    }

    i++;
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

  double x, y, r, th;
  for (i = 0; i < N; i++) {
    // Read the current observation
    x = gsl_matrix_get(os, 0, i);
    y = gsl_matrix_get(os, 1, i);

    // Convert to polar coordinates and apply the L
    r = L - sqrt(pow(x, 2) + pow(y, 2));
    th = atan2(y, x);

    // Return to rectangular coordinates
    x = r * cos(th);
    y = r * sin(th);

    // Store it in the matrix again
    gsl_matrix_set(os, 0, i, x);
    gsl_matrix_set(os, 1, i, y);
  }
}

void
toLaserSpace(obsmat os, double L, int sign) {
  size_t N = os->size2;
  size_t i;

  double x, y, r, th;
  for (i = 0; i < N; i++) {
    // Read the current observation
    x = gsl_matrix_get(os, 0, i);
    y = gsl_matrix_get(os, 1, i);

    // Convert to polar coordinates and apply the L
    r = L - sqrt(pow(x, 2) + pow(y, 2));
    th = atan2(y, x);

    // Return to rectangular coordinates
    x = r * cos(th);
    y = r * sin(th);

    // Store it in the matrix again
    gsl_matrix_set(os, 0, i, x);
    gsl_matrix_set(os, 1, i, y);
  }
}
