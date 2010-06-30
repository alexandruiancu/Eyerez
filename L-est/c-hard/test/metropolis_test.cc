#include <metropolis.h>
#include "../problem.h"
#include <gtest/gtest.h>

#include <stdio.h>
#include <gsl/gsl_matrix.h>

class testWithMetropolis : public ::testing::Test {
public:  
  gsl_vector *mode;
  gsl_matrix *os;
  Metropolis m1;

  virtual void SetUp() {
    mode = gsl_vector_alloc(13);
    gsl_vector_set(mode, 0, 30.4);
    gsl_vector_set(mode, 1, 0.1);
    gsl_vector_set(mode, 2, 0.0);
    gsl_vector_set(mode, 3, 0.0);
    gsl_vector_set(mode, 4, 0.0);
    gsl_vector_set(mode, 5, 0.1);
    gsl_vector_set(mode, 6, 1/9.0);
    gsl_vector_set(mode, 7, 1/9.0);
    gsl_vector_set(mode, 8, 1/9.0);
    gsl_vector_set(mode, 9, 0.0);
    gsl_vector_set(mode, 10, 0.0);
    gsl_vector_set(mode, 11, 0.0);
    gsl_vector_set(mode, 12, 0.2);

    os = allocObsCSV(fopen("e2.csv", "r"));
    m1 = Metropolis (logPostVectorized, 13, (void *)os);
  }

  virtual void TearDown() {
    gsl_matrix_free(os);
    gsl_vector_free(mode);
  }
};

TEST_F (testWithMetropolis, core) {
  m1.jump(mode);
}
