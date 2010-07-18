#include "problem.hpp"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_blas.h>

#include <classes/observations.hpp>
#include <classes/random.hpp>

namespace MyUtilities {

  parameter::parameter () {
    initialized = false;
    edited = false;
  }
  
  double 
  parameter::operator() () const {
    return self;
  }
 
  void 
  parameter::operator() (const double new_value) {
    self = new_value;
    if (!initialized) { initialized = true; }
    edited = true;
  }

  bool
  parameter::isInitialized () {
    return initialized;
  }

  bool
  parameter::isEdited () {
    return edited;
  }

  void
  parameter::clean () {
    edited = false;
  }

}

namespace LEst {
  Estimate::Estimate (Observations *obs_in) {
    observations = obs_in;
    D = 15;
  }

  Estimate::~Estimate () {
    delete observations;
  }

  Estimate *
  Estimate::clone () const {
    Estimate *out = new Estimate (observations);
    gsl_vector *a_pickle = gsl_vector_alloc(15);
    pickle(a_pickle);
    out->unpickle(a_pickle);
    gsl_vector_free(a_pickle);

    return out;
  }
  
  bool
  Estimate::isFullyInitialized () {
    bool out = true;
    out &= L.isInitialized();
    out &= A.isInitialized();
    out &= B.isInitialized();
    out &= C.isInitialized();
    out &= rot1.isInitialized();
    out &= rot2.isInitialized();
    out &= rot3.isInitialized();
    out &= cx.isInitialized();
    out &= cy.isInitialized();
    out &= cz.isInitialized();
    out &= sdCx.isInitialized();
    out &= sdCy.isInitialized();
    out &= sdCz.isInitialized();
    out &= sdR.isInitialized();
    out &= sdZ.isInitialized();

    return out;
  }

  bool
  Estimate::isValid () {
    bool ok = true;

    // SDs are positive
    ok &= sdCx() > 0;
    ok &= sdCy() > 0;
    ok &= sdCz() > 0;
    ok &= sdR()  > 0;
    ok &= sdZ()  > 0;

    // Rotations are valid
    ok &= (rot1() > -M_PI && rot1() < M_PI);
    ok &= (rot2() > 0     && rot2() < M_PI);
    ok &= (rot3() > -M_PI && rot3() < M_PI);
    
    // L is positive
    ok &= L() > 0;

    // Major axes are positive
    ok &= A() > 0;
    ok &= B() > 0;
    ok &= C() > 0;

    return ok;
  }

  double
  Estimate::logPost () {
    if (isValid()) {
      repair();
      return (logLik + priorQ + priorLambda + priorErr +
              priorC + priorCErr + priorL);
    } else {
      return -INFINITY;
    }
  }

  void
  Estimate::repair () {
    // repair the likelihood
    if (L.isEdited() ||
        A.isEdited() ||
        B.isEdited() ||
        C.isEdited() ||
        rot1.isEdited() ||
        rot2.isEdited() ||
        rot3.isEdited() ||
        cx.isEdited() ||
        cy.isEdited() ||
        cz.isEdited() ||
        sdR.isEdited() ||
        sdZ.isEdited()) {

      logLik = computeLogLik();

    }


    if (rot1.isEdited() ||
        rot2.isEdited() ||
        rot3.isEdited()) {
      
      priorQ = computePriorQ();

    }


    if (A.isEdited() ||
        B.isEdited() ||
        C.isEdited()) {

      priorLambda = computePriorLambda();

    }


    if (sdR.isEdited() || sdZ.isEdited()) {
      priorErr = computePriorErr();
    }

    if (cx.isEdited() || cy.isEdited() || cz.isEdited() ||
        sdCx.isEdited() ||
        sdCy.isEdited() ||
        sdCz.isEdited()) {
      priorC = computePriorC();
    }

    if (sdCx.isEdited() ||
        sdCy.isEdited() ||
        sdCz.isEdited()) {
      priorCErr = computePriorCErr();
    }

    if (L.isEdited()) {
      priorL = computePriorL();
    }

    L.clean(); 
    A.clean();
    B.clean();
    C.clean();
    rot1.clean();
    rot2.clean();
    rot3.clean();
    cx.clean();
    cy.clean();
    cz.clean();
    sdCx.clean();
    sdCy.clean();
    sdCz.clean();
    sdR.clean();
    sdZ.clean();

  }

  void
  Estimate::pickle (gsl_vector *out) const {
    if (out->size != 15) {
      fprintf(stderr, "Could not pickle State: pickling vector is of wrong dimension (!= 15)");
    }

    gsl_vector_set(out, 0, L());
    gsl_vector_set(out, 1, A());
    gsl_vector_set(out, 2, B());
    gsl_vector_set(out, 3, C());
    gsl_vector_set(out, 4, rot1());
    gsl_vector_set(out, 5, rot2());
    gsl_vector_set(out, 6, rot3());
    gsl_vector_set(out, 7, cx());
    gsl_vector_set(out, 8, cy());
    gsl_vector_set(out, 9, cz());
    gsl_vector_set(out, 10, sdCx());
    gsl_vector_set(out, 11, sdCy());
    gsl_vector_set(out, 12, sdCz());
    gsl_vector_set(out, 13, sdR());
    gsl_vector_set(out, 14, sdZ());

  }

  void
  Estimate::unpickle (const gsl_vector *in) {
    if (in->size != 15) {
      fprintf(stderr, "Could not unpickle State: pickled vector is of wrong dimension (!= 15)");
    }

    L(   gsl_vector_get(in, 0));
    A(   gsl_vector_get(in, 1));
    B(   gsl_vector_get(in, 2));
    C(   gsl_vector_get(in, 3));
    rot1(gsl_vector_get(in, 4));
    rot2(gsl_vector_get(in, 5));
    rot3(gsl_vector_get(in, 6));
    cx(  gsl_vector_get(in, 7));
    cy(  gsl_vector_get(in, 8));
    cz(  gsl_vector_get(in, 9));
    sdCx(gsl_vector_get(in, 10));
    sdCy(gsl_vector_get(in, 11));
    sdCz(gsl_vector_get(in, 12));
    sdR( gsl_vector_get(in, 13));
    sdZ( gsl_vector_get(in, 14));

  }



  double
  Estimate::computeLogLik() {
      gsl_vector *lambda = gsl_vector_alloc(3);

      gsl_vector_set(lambda, 0, A());
      gsl_vector_set(lambda, 1, B());
      gsl_vector_set(lambda, 2, C());
      
      gsl_matrix *Q = gsl_matrix_alloc(3, 3);
      double ca = cos(rot1());
      double sa = sin(rot1());
      double cb = cos(rot2());
      double sb = sin(rot2());
      double cc = cos(rot3());
      double sc = sin(rot3());      
      gsl_matrix_set(Q, 0, 0, ca*cc - cb*sa*sc);
      gsl_matrix_set(Q, 0, 1, cc*sa + ca*cb*sc);
      gsl_matrix_set(Q, 0, 2, sb*sc);
      gsl_matrix_set(Q, 1, 0, -cb*cc*sa - ca*sc);
      gsl_matrix_set(Q, 1, 1, ca*cb*cc - sa*sc);
      gsl_matrix_set(Q, 1, 2, cc*sb);
      gsl_matrix_set(Q, 2, 0, sa*sb);
      gsl_matrix_set(Q, 2, 1, -ca*sb);
      gsl_matrix_set(Q, 2, 2, cb);

      gsl_vector *c = gsl_vector_alloc(3);
      gsl_vector_set(c, 0, cx());
      gsl_vector_set(c, 0, cy());
      gsl_vector_set(c, 0, cz());

      double logLik = datasetLogLik(observations, L(), lambda, Q, c, sdR(), sdZ());

      gsl_vector_free(c);
      gsl_matrix_free(Q);
      gsl_vector_free(lambda);

      return logLik;
  }

  double
  Estimate::computePriorQ () {
    double logPrior = 0;

    const double ROT_SD = M_PI/8;
    logPrior += log(gsl_ran_gaussian_pdf(rot1(), ROT_SD));
    logPrior += log(gsl_ran_gaussian_pdf(abs(rot2()), ROT_SD/2));
    logPrior += log(gsl_ran_gaussian_pdf(rot3(), ROT_SD));

    return logPrior;
  }

  double
  Estimate::computePriorLambda () {
    double logPrior = 0;

    const double AXES_MU = 1750;
    const double AXES_SD = 500;
    logPrior += log(gsl_ran_gaussian_pdf(A() - AXES_MU, AXES_SD));
    logPrior += log(gsl_ran_gaussian_pdf(B() - AXES_MU, AXES_SD));
    logPrior += log(gsl_ran_gaussian_pdf(C() - AXES_MU, AXES_SD));
    return logPrior;
  }

  double
  Estimate::computePriorErr () {
    double logPrior = 0;

    logPrior -= log(sdR());
    logPrior -= log(sdZ());
    return logPrior;
  }

  double
  Estimate::computePriorC () {
    double logPrior = 0;

    logPrior += log(gsl_ran_gaussian_pdf(cx(), sqrt(sdCx())));
    logPrior += log(gsl_ran_gaussian_pdf(cy(), sqrt(sdCy())));
    logPrior += log(gsl_ran_gaussian_pdf(cz(), sqrt(sdCz())));

    return logPrior;
  }

  double
  Estimate::computePriorCErr () {
    double logPrior = 0;

    logPrior -= log(sdCx());
    logPrior -= log(sdCy());
    logPrior -= log(sdCz());

    return logPrior;
  }

  double
  Estimate::computePriorL () {
    return log(gsl_ran_gaussian_pdf(L() - 30400, 1500));
  }
  
  double
  Estimate::datasetLogLik(Observations *obs,
                          const double L,
                          const gsl_vector *lambda,
                          const gsl_matrix *Q,
                          const gsl_vector *c,
                          const double sdR,
                          const double sdZ) {
    size_t N = obs->N;
    size_t D = obs->D;

    Random *R = new Random();
    
    double aveLogLik = 0;
    gsl_vector *uvw    = gsl_vector_alloc(D);
    gsl_vector *ytmp   = gsl_vector_alloc(D);
    gsl_vector *ob     = gsl_vector_alloc(D);

    /* Bootstrap with replacement NSAMPLES samples from the
       observation matrix to find an estimate for the average logLik
       then scale it by N/NSAMPLES to get an estimate for the total
       logLik of the dataset.

       This introduces the possibility of instability of the estimate,
       but so long as NSAMPLES is suitable large (TODO: determine how
       large it must be) this should converge well.

       Additionally, if any of the computations fail to converge for
       any reason, do not include them in the sample. I believe this
       will eliminate points that are "near" any of the orthogonal
       planes in the system but I believe "near" must mean
       "numerically zero" and therefore with eyes around 3mm in
       diameter the majority of the points should be included in the
       sampling population.

       TODO: check above assumption.
    */

    size_t nchosen = 0;
    do {
      nchosen++;

      // Get a random observation
      obs->writeObservation(R->drawUniformIndex(N), ob);     

      // Try to find the likelihood, but if it fails forget about it.
      try {
      aveLogLik +=
        _observationLogLik(ob, lambda, Q, c, sdR, sdZ, uvw, ytmp);
      } catch(int) { nchosen--; }
      
    }
    while (nchosen < NSAMPLES);

    gsl_vector_free(ob);
    gsl_vector_free(ytmp);
    gsl_vector_free(uvw);
    delete R;
    
    return exp(log(N) - log(NSAMPLES))*aveLogLik;
  }

  double 
  Estimate::findT(const double a, const double b, const double c,
                  const double u, const double v, const double w) {
    
    /* We're looking for the largest root of F(t),
       
       (t+a^2)^2 (t+b^2)^2 (t+c^2)^2 
       - a^2 u^2 (t+b^2)^2 (t+c^2)^2 
       - b^2 v^2 (t+a^2)^2 (t+c^2)^2 
       - c^2 w^2 (t+a^2)^2 (t+b^2)^2
       
       which is very complex to expand...
       
       a^4 b^4 c^4+2 a^4 b^4 c^2 t+a^4 b^4 t^2+2 a^4 b^2 c^4 t+4 a^4 b^2 c^2 t^2+2 a^4 b^2 t^3+a^4 c^4 t^2+2 a^4 c^2 t^3+a^4 t^4+2 a^2 b^4 c^4 t+4 a^2 b^4 c^2 t^2+2 a^2 b^4 t^3+4 a^2 b^2 c^4 t^2+8 a^2 b^2 c^2 t^3+4 a^2 b^2 t^4+2 a^2 c^4 t^3+4 a^2 c^2 t^4+2 a^2 t^5+b^4 c^4 t^2+2 b^4 c^2 t^3+b^4 t^4+2 b^2 c^4 t^3+4 b^2 c^2 t^4+2 b^2 t^5+c^4 t^4+2 c^2 t^5+t^6
       
       - a^2 b^4 c^4 u^2+2 a^2 b^4 c^2 t u^2+a^2 b^4 t^2 u^2+2 a^2 b^2 c^4 t u^2+4 a^2 b^2 c^2 t^2 u^2+2 a^2 b^2 t^3 u^2+a^2 c^4 t^2 u^2+2 a^2 c^2 t^3 u^2+a^2 t^4 u^2
       
       - a^4 b^2 c^4 v^2+2 a^4 b^2 c^2 t v^2+a^4 b^2 t^2 v^2+2 a^2 b^2 c^4 t v^2+4 a^2 b^2 c^2 t^2 v^2+2 a^2 b^2 t^3 v^2+b^2 c^4 t^2 v^2+2 b^2 c^2 t^3 v^2+b^2 t^4 v^2
       
       - a^4 b^4 c^2 w^2+2 a^4 b^2 c^2 t w^2+a^4 c^2 t^2 w^2+2 a^2 b^4 c^2 t w^2+4 a^2 b^2 c^2 t^2 w^2+2 a^2 c^2 t^3 w^2+b^4 c^2 t^2 w^2+2 b^2 c^2 t^3 w^2+c^2 t^4 w^2
       
       
    */
    
    double a2 = pow(a, 2);
    double a4 = pow(a, 4);
    double b2 = pow(b, 2);
    double b4 = pow(b, 4);
    double c2 = pow(c, 2);
    double c4 = pow(c, 4);
    double u2 = pow(u, 2);
    double v2 = pow(v, 2);
    double w2 = pow(w, 2);
    double mod = a2 + b2 + c2;
    double prod = a2*b2*c2;
    
    double const_term, t1, t2, t3, t4, t5, t6;
    
    const_term  = 1 - u2/a2 - v2/b2 - w2/c2;
    const_term *= pow(prod, 2);
    
    t1  = a2*b2 + a2*c2 + b2*c2;
    t1 -= b2*u2 + c2*u2 + a2*v2 + c2*v2 + a2*w2 + b2*w2;
    t1 *= 2 * prod;
    
    t2  = a4*b4 + a4*c4 + b4*c4;
    t2 += 4 * prod * mod;
    t2 -= 4 * prod * (u2 + v2 + w2);
    t2 -= a2*b4*u2 + a4*b2*v2 + b2*c4*v2 + b4*c2*w2 + a2*c4*u2 + a4*c2*w2;
    
    t3  = 2 * (a4*b2 + a4*c2 + a2*b4 + a2*c4 + b4*c2 + b2*c4);
    t3 += 8 * prod;
    t3 -= 2 * prod * (u2/c2 + u2/b2 + v2/c2 + v2/a2 + w2/b2 + w2/a2);
    
    t4  = a4 + b4 + c4;
    t4 += 4 * (a2*b2 + a2*c2 + b2*c2);
    t4 -= a2*u2 + b2*v2 + c2*w2;
    
    t5  = 2*mod;
    
    t6  = 1;
    
    const size_t N = 7;
    double coefs[N] = { const_term, t1, t2, t3, t4, t5, t6 };
    double z[2*(N-1)] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    gsl_poly_complex_workspace * workspace
      = gsl_poly_complex_workspace_alloc(N);

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    int err = gsl_poly_complex_solve(coefs, N, workspace, z);
    gsl_set_error_handler(old_handler);

    if (err) throw err;
    
    gsl_poly_complex_workspace_free(workspace);
    
    double real, complex, out = 0;
    for (size_t i = 0; i < 6; i++) {
      real = z[2*i];
      complex = z[2*i + 1];
      if (real > out) { out = real; }
      // if (abs(complex) > 0) { throw; } 
    }
    
    return out;
  }
  
  void
  Estimate::computeDelta(const double a, const double b, const double c,
                         const double u, const double v, const double w,
               gsl_vector *out) {
    
    double t = findT(a, b, c, u, v, w);
    
    gsl_vector_set(out, 0, u*t/(t + pow(a, 2)));
    gsl_vector_set(out, 1, v*t/(t + pow(b, 2)));
    gsl_vector_set(out, 2, w*t/(t + pow(c, 2)));
  }
  
  double
  Estimate::ellipsoidLogPdf(const double a, const double b, const double c,
                            const double u, const double v, const double w,
                            const double sdR, 
                            const double sdZ) {
    
    gsl_vector *diff = gsl_vector_alloc(3);
    computeDelta(a, b, c, u, v, w, diff);

    double logp;
    logp  = log(gsl_ran_gaussian_pdf(gsl_vector_get(diff, 0), sdR));
    logp += log(gsl_ran_gaussian_pdf(gsl_vector_get(diff, 1), sdR));
    logp += log(gsl_ran_gaussian_pdf(gsl_vector_get(diff, 2), sdZ));
    
    gsl_vector_free(diff);
    
    return logp;
  }
  
  double
  Estimate::observationLogLik(const gsl_vector *ob,
                              const gsl_vector *lambda,
                              const gsl_matrix *Q,
                              const gsl_vector *c,
                              const double sdR,
                              const double sdZ) {

    size_t D = ob->size;
    
    gsl_vector *uvw    = gsl_vector_alloc(D);
    gsl_vector *ytmp   = gsl_vector_alloc(D);
    
    double logLik =
      _observationLogLik(ob, lambda, Q, c, sdR, sdZ, uvw, ytmp);
    
    gsl_vector_free(ytmp);
    gsl_vector_free(uvw);

    return logLik;
  } 
  
  double
  Estimate::_observationLogLik(const gsl_vector *ob,
                               const gsl_vector *lambda,
                               const gsl_matrix *Q,
                               const gsl_vector *c,
                               const double sdR,
                               const double sdZ,
                               gsl_vector *uvw,
                               gsl_vector *ytmp) {
    
    // !! These are some of the more expensive operations
    // Shift and rotate to the eigenspace
    // Qt (y - c) ---> (u, v, w)
    gsl_vector_memcpy(ytmp, c); // ytmp <- c
    gsl_vector_scale(ytmp, -1); // ytmp <- -o
    gsl_vector_add(ytmp, ob);   // ytmp <- ytmp + (an observation)
                                // [ytmp <- (-c) + (an observation)]
    gsl_vector_set_zero(uvw);
    gsl_blas_dgemv(CblasTrans,
                   1.0, Q, ytmp,
                   0.0, uvw);   // uvw <- Qt * ytmp
                                // [uvs <- Qt ((an obs) - c)]
    
    return
      ellipsoidLogPdf(1.0/sqrt(gsl_vector_get(lambda, 0)), // a
                      1.0/sqrt(gsl_vector_get(lambda, 1)), // b
                      1.0/sqrt(gsl_vector_get(lambda, 2)), // c
                      gsl_vector_get(uvw, 0),              // u
                      gsl_vector_get(uvw, 1),              // v
                      gsl_vector_get(uvw, 2),              // w
                      sdR, sdZ);
  }

  gsl_vector *view_ell(const State *st) {
    Estimate *est = (Estimate *)st;
    gsl_vector *out = gsl_vector_alloc(1);
    gsl_vector_set(out, 0, est->L());
    return out;
  }
  
  gsl_vector *view_lambda(const State *st) {
    Estimate *est = (Estimate *)st;
    gsl_vector *out = gsl_vector_alloc(3);
    gsl_vector_set(out, 0, est->A());
    gsl_vector_set(out, 1, est->B());
    gsl_vector_set(out, 2, est->C());
    return out;
  }
  
  gsl_vector *view_rot(const State *st) {
    Estimate *est = (Estimate *)st;
    gsl_vector *out = gsl_vector_alloc(3);
    gsl_vector_set(out, 0, est->rot1());
    gsl_vector_set(out, 1, est->rot2());
    gsl_vector_set(out, 2, est->rot3());
    return out;
  }
  
  gsl_vector *view_shift(const State *st) {
    Estimate *est = (Estimate *)st;
    gsl_vector *out = gsl_vector_alloc(3);
    gsl_vector_set(out, 0, est->cx());
    gsl_vector_set(out, 1, est->cy());
    gsl_vector_set(out, 2, est->cz());
    return out;
  }
  
  gsl_vector *view_cerr(const State *st) {
    Estimate *est = (Estimate *)st;
    gsl_vector *out = gsl_vector_alloc(3);
    gsl_vector_set(out, 0, est->sdCx());
    gsl_vector_set(out, 1, est->sdCy());
    gsl_vector_set(out, 2, est->sdCz());
    return out;
  }
  
  gsl_vector *view_err(const State *st) {
    Estimate *est = (Estimate *)st;
    gsl_vector *out = gsl_vector_alloc(2);
    gsl_vector_set(out, 0, est->sdR());
    gsl_vector_set(out, 1, est->sdZ());
    return out;
  }
  
  void update_ell(const gsl_vector * vec, State *st) {
    Estimate *est = (Estimate *)st;
    est->L(gsl_vector_get(vec, 0));
  }
  
  void update_lambda(const gsl_vector * vec, State *st) {
    Estimate *est = (Estimate *)st;
    est->A(gsl_vector_get(vec, 0));
    est->B(gsl_vector_get(vec, 1));
    est->C(gsl_vector_get(vec, 2));
  }
  
  void update_rot(const gsl_vector * vec, State *st) {
    Estimate *est = (Estimate *)st;
    est->rot1(gsl_vector_get(vec, 0));
    est->rot2(gsl_vector_get(vec, 1));
    est->rot3(gsl_vector_get(vec, 2));
  }
  
  void update_shift(const gsl_vector * vec, State *st) {
    Estimate *est = (Estimate *)st;
    est->cx(gsl_vector_get(vec, 0));
    est->cy(gsl_vector_get(vec, 1));
    est->cz(gsl_vector_get(vec, 2));
  }
  
  void update_cerr(const gsl_vector * vec, State *st) {
    Estimate *est = (Estimate *)st;
    est->sdCx(gsl_vector_get(vec, 0));
    est->sdCy(gsl_vector_get(vec, 1));
    est->sdCz(gsl_vector_get(vec, 2));
  }
  
  void update_err(const gsl_vector * vec, State *st) {
    Estimate *est = (Estimate *)st;
    est->sdR(gsl_vector_get(vec, 0));
    est->sdZ(gsl_vector_get(vec, 1));
  }
    
  // TODO: Update this to truly random initialization.
  
  gsl_vector *guess0_ell() {
    Random *r = new Random();
    gsl_vector *init = gsl_vector_alloc(1);
    gsl_vector_set(init, 0, r->drawGaussian(30400, 1000));

    delete r;
    return init;
  }
  
  gsl_vector *guess0_lambda() {
    Random *r = new Random();
    gsl_vector *init = gsl_vector_alloc(3);
    gsl_vector_set(init, 0, r->drawGaussian(1750, 250));
    gsl_vector_set(init, 1, r->drawGaussian(1750, 250));
    gsl_vector_set(init, 2, r->drawGaussian(1750, 250));

    delete r;
    return init;
  }
  
  gsl_vector *guess0_rot() {
    Random *r = new Random();
    gsl_vector *init = gsl_vector_alloc(3);
    gsl_vector_set(init, 0, 0.001);
    gsl_vector_set(init, 1, 0.001);
    gsl_vector_set(init, 2, 0.001);

    delete r;
    return init;
  }
  
  gsl_vector *guess0_shift() {
    Random *r = new Random();
    gsl_vector *init = gsl_vector_alloc(3);
    gsl_vector_set(init, 0, r->drawGaussian(0, 1000));
    gsl_vector_set(init, 1, r->drawGaussian(0, 1000));
    gsl_vector_set(init, 2, r->drawGaussian(0, 1000));

    delete r;
    return init;
  }
  
  gsl_vector *guess0_cerr() {
    Random *r = new Random();
    gsl_vector *init = gsl_vector_alloc(3);
    gsl_vector_set(init, 0, r->drawGaussian(1000.0, 200.0));
    gsl_vector_set(init, 1, r->drawGaussian(1000.0, 200.0));
    gsl_vector_set(init, 2, r->drawGaussian(1000.0, 200.0));

    delete r;
    return init;
  }
  
  gsl_vector *guess0_err() {
    Random *r = new Random();
    gsl_vector *init = gsl_vector_alloc(2);
    gsl_vector_set(init, 0, r->drawGaussian(100.0, 20.0));
    gsl_vector_set(init, 1, r->drawGaussian(100.0, 20.0));

    delete r;
    return init;
  }

}
