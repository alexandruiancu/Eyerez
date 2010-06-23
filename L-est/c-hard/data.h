#include <stdio.h>
#include <math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <csv_parser/csv_parser.hpp>

typedef gsl_matrix *obsmat;

obsmat 
allocObsCSV(FILE *f);

int 
countlines(FILE *f);

void
toRealSpace(obsmat os, double L, int sign);

void
toLaserSpace(obsmat os, double L, int sign);
