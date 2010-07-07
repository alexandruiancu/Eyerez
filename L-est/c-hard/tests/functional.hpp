#include <cxxtest/TestSuite.h>

#include <classes/metropolis.hpp>
#include <classes/problem.hpp>

class MetropolisSetup : public CxxTest::TestSuite {
public:
  
  void testMakeMetropolis() {

    FILE *f = fopen("tests/fixtures/e-rect.csv", "r");
    TS_ASSERT (f);

    Observations *obs = 
      new Observations(f, RectCoord);
    
    Problem::State s;
    s.obs = obs;
    gsl_vector *init = s.allocVectorized();

    Metropolis *m1 = new Metropolis(Problem::logPost, init, obs, 0.01);

    delete m1;
    gsl_vector_free(init);
    delete obs;
  }
};

class MetropolisFunction : public CxxTest::TestSuite {
public:
  Observations *obs;
  Problem::State s;
  gsl_vector *init;
  Metropolis *m1;

  void setUp() { 
    obs = new Observations(fopen("sphereoid-polar-low.csv", "r"), PolarCoord);

    Problem::State s;
    s.obs = obs;
    init = s.allocVectorized();

    m1 = new Metropolis (Problem::logPost, init, obs, 0.01);
  }

  void tearDown() { 
    delete m1;
    gsl_vector_free(init);    
    delete obs;
  }

  void testDidJump() {
    for (size_t i = 0; i < 200; i++) { 
      TS_ASSERT_THROWS_NOTHING(m1->jump()); 
    }
    m1->freeze();
    for (size_t i = 0; i < 500; i++) { 
      TS_ASSERT_THROWS_NOTHING(m1->jump()); 
    }

    double initEle, currentEle;
    for (size_t i = 0; i < 17; i++) {
      initEle = gsl_vector_get(m1->state, i);
      currentEle = gsl_vector_get(init, i);
      TS_ASSERT_DIFFERS(currentEle, initEle);

    }
  }
};
