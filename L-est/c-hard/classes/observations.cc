#include "observations.hpp"

void
Observations::init(const gsl_matrix *obs, const CoordSystem system_in) {
  N = obs->size1;
  D = obs->size2;

  data = gsl_matrix_alloc(N, D);
  gsl_matrix_memcpy(data, obs);

  system = system_in;
  setL();
}

Observations::Observations(const gsl_matrix *obs, const CoordSystem system_in) {
  init(obs, system_in);
}

Observations::Observations(FILE *f, const CoordSystem system_in) {
  N = Util::countlines(f);
  D = 3;

  gsl_matrix *obs = gsl_matrix_alloc(N, D);

  csv_parser file_parser;
  file_parser.init(f);
  file_parser.set_field_term_char(',');
  file_parser.set_line_term_char('\n');

  int i = 0;

  while (file_parser.has_more_rows()) {
    csv_row row = file_parser.get_row();
    
    for (unsigned int j = 0; j < row.size(); j++) {
      gsl_matrix_set(obs, i, j, atof(row[j].c_str()));
    }

    i++;
  }

  init(obs, system_in);
}

Observations::~Observations() {
  gsl_matrix_free(data);
}

void
Observations::setL(const double L_in) { L = L_in; }

void
Observations::writeObservation(const size_t n, gsl_vector *o) {
  double x, y, z, r, th;
  if (system == RectCoord) {
    x  = gsl_matrix_get(data, n, 0);
    y  = gsl_matrix_get(data, n, 1);
    z  = gsl_matrix_get(data, n, 2);

    r  = L - sqrt(x*x + y*y);
    th = atan2(y, x);
  } else {
    r  = L - gsl_matrix_get(data, n, 0);
    th = gsl_matrix_get(data, n, 1);
    z  = gsl_matrix_get(data, n, 2);
  }
  gsl_vector_set(o, 0, r * cos(th));
  gsl_vector_set(o, 1, r * sin(th));
  gsl_vector_set(o, 2, z);
}

gsl_vector *
Observations::allocAsObservation(const size_t n) {
  gsl_vector *o = gsl_vector_alloc(D);
  writeObservation(n, o);
  return o;
}

namespace Util {

  size_t countlines(FILE *f) {
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
}
