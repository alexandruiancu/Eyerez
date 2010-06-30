#include "problem.hpp"

using namespace std;
namespace Problem {

  double
  logPost(const gsl_vector *x, Observations *obs_in) {
    State s (x, obs_in);
    return s.logLik();
  }

  State::State() {
    Random *R = new Random();

    // muL  = R->drawGaussian(30.4, 1);
    // sdL  = R->drawGamma(0.01, 1);
    // L    = R->drawGaussian(30.4, 1);
    // A    = R->drawGaussian(1.75, 0.5);
    // B    = R->drawGaussian(1.75, 0.5);
    // C    = R->drawGaussian(1.75, 0.5);
    // rot1 = R->drawGaussian(0, M_PI/8);
    // rot2 = R->drawGaussian(0, M_PI/16);
    // rot3 = R->drawGaussian(0, M_PI/8);
    // sdCx = R->drawGamma(0.01, 1);
    // sdCy = R->drawGamma(0.01, 1);
    // sdCz = R->drawGamma(0.01, 1);
    // cx   = R->drawGaussian(0, 1);
    // cx   = R->drawGaussian(0, 1);
    // cx   = R->drawGaussian(0, 1);
    // sdR  = R->drawGamma(0.01, 1);
    // sdZ  = R->drawGamma(0.01, 1);

    muL  = 30.4;
    sdL  = 0.1;
    L    = 30.4;
    A    = 1.75;
    B    = 1.75;
    C    = 1.75;
    rot1 = 0;
    rot2 = 0;
    rot3 = 0;
    sdCx = 0.1;
    sdCy = 0.1;
    sdCz = 0.1;
    cx   = 0;
    cx   = 0;
    cx   = 0;
    sdR  = 0.1;
    sdZ  = 0.1;

    delete R;
  }

  State::State(const gsl_vector *x, Observations *obs_in) {
    L    = gsl_vector_get(x, 0);  // V1
    muL  = gsl_vector_get(x, 1);  // V2
    sdL  = gsl_vector_get(x, 2);  // V3

    A    = gsl_vector_get(x, 3);  // V4
    B    = gsl_vector_get(x, 4);  // V5
    C    = gsl_vector_get(x, 5);  // V6
    rot1 = gsl_vector_get(x, 6);  // V7
    rot2 = gsl_vector_get(x, 7);  // V8
    rot3 = gsl_vector_get(x, 8);  // V9

    cx   = gsl_vector_get(x, 9);  // V10
    cy   = gsl_vector_get(x, 10); // V11
    cz   = gsl_vector_get(x, 11); // V12
    sdCx = gsl_vector_get(x, 12); // V13
    sdCy = gsl_vector_get(x, 13); // V14
    sdCz = gsl_vector_get(x, 14); // V15
    
    sdR  = gsl_vector_get(x, 15); // V16
    sdZ  = gsl_vector_get(x, 16); // V17

    obs = obs_in;
  }

  gsl_vector *
  State::allocVectorized() {
    gsl_vector *vec = gsl_vector_alloc(17);
    double load[17] = {L, muL, sdL, 
                       A, B, C, 
                       rot1, rot2, rot3, 
                       cx, cy, cz,
                       sdCx, sdCy, sdCz,
                       sdR, sdZ};
    for (size_t i = 0; i < 17; i++) { gsl_vector_set(vec, i, load[i]); }

    return vec;
  }

  int
  State::isValid() {
    int ok = 1;

    // SDs are positive
    ok = ok && sdL  > 0;
    ok = ok && sdCx > 0;
    ok = ok && sdCy > 0;
    ok = ok && sdCz > 0;
    ok = ok && sdR  > 0;
    ok = ok && sdZ  > 0;

    // Rotations are valid
    ok = ok && (rot1 > -M_PI && rot1 < M_PI);
    ok = ok && (rot2 > 0     && rot2 < M_PI);
    ok = ok && (rot3 > -M_PI && rot3 < M_PI);
    
    // L and muL are positive
    ok = ok && L > 0;
    ok = ok && muL > 0;

    // Major axes are positive
    ok = ok && A > 0;
    ok = ok && B > 0;
    ok = ok && C > 0;

    return ok;
  }

  double
  State::logLik() {
    if (isValid()) {
      double logLik = 0;

      gsl_vector *lambda = gsl_vector_alloc(3);

      gsl_vector_set(lambda, 0, A);
      gsl_vector_set(lambda, 1, B);
      gsl_vector_set(lambda, 2, C);
      
      gsl_matrix *Q = gsl_matrix_alloc(3, 3);
      double ca = cos(rot1);
      double sa = sin(rot1);
      double cb = cos(rot2);
      double sb = sin(rot2);
      double cc = cos(rot3);
      double sc = sin(rot3);      
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
      gsl_vector_set(c, 0, cx);
      gsl_vector_set(c, 0, cy);
      gsl_vector_set(c, 0, cz);

      logLik += datasetLogLik(obs, L, lambda, Q, c, sdR, sdZ);
      logLik += priors();
      
      return logLik;
    } else {
      return -INFINITY;
    }
  }

  double
  State::priors() {
    double logLik = 0;

    // Rotations
    const double ROT_SD = M_PI/8;
    logLik += log(gsl_ran_gaussian_pdf(rot1, ROT_SD));
    logLik += log(gsl_ran_gaussian_pdf(abs(rot2), ROT_SD/2));
    logLik += log(gsl_ran_gaussian_pdf(rot3, ROT_SD));

    // Axes
    const double AXES_MU = 1.75;
    const double AXES_SD = 0.5;;
    logLik += log(gsl_ran_gaussian_pdf(A - AXES_MU, AXES_SD));
    logLik += log(gsl_ran_gaussian_pdf(B - AXES_MU, AXES_SD));
    logLik += log(gsl_ran_gaussian_pdf(C - AXES_MU, AXES_SD));

    // Measurement error
    logLik -= log(sdR);
    logLik -= log(sdZ);

    // Shift
    logLik += log(gsl_ran_gaussian_pdf(cx, sqrt(sdCx)));
    logLik += log(gsl_ran_gaussian_pdf(cy, sqrt(sdCy)));
    logLik += log(gsl_ran_gaussian_pdf(cz, sqrt(sdCz)));

    // Shift variance
    logLik -= 10*log(sdCx);
    logLik -= 10*log(sdCy);
    logLik -= 10*log(sdCz);

    // L
    logLik += log(gsl_ran_gaussian_pdf(L - muL, sdL));
    logLik += log(gsl_ran_gaussian_pdf(muL - 30.4, 1.5));
    logLik -= log(sdL);

    return logLik;
  }


  double 
  findT(const double a, const double b, const double c,
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
    double z[2*(N-1)];
    
    gsl_poly_complex_workspace * workspace
      = gsl_poly_complex_workspace_alloc(N);
    
    gsl_poly_complex_solve(coefs, N, workspace, z);
    
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
  computeDelta(const double a, const double b, const double c,
               const double u, const double v, const double w,
               gsl_vector *out) {
    
    double t = findT(a, b, c, u, v, w);
    
    gsl_vector_set(out, 0, u*t/(t + pow(a, 2)));
    gsl_vector_set(out, 1, v*t/(t + pow(b, 2)));
    gsl_vector_set(out, 2, w*t/(t + pow(c, 2)));
  }
  
  double
  ellipsoidLogPdf(const double a, const double b, const double c,
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
  obserationLogLik(const gsl_vector *ob,
                   const gsl_vector *lambda,
                   const gsl_matrix *Q,
                   const gsl_vector *c,
                   const double sdR,
                   const double sdZ) {

    size_t D = ob->size;
    
    gsl_vector *uvw    = gsl_vector_alloc(D);
    gsl_vector *ytmp   = gsl_vector_alloc(D);
    
    double logLik =
      _obserationLogLik(ob, lambda, Q, c, sdR, sdZ, uvw, ytmp);
    
    gsl_vector_free(ytmp);
    gsl_vector_free(uvw);

    return logLik;
  } 
  
  double
  _obserationLogLik(const gsl_vector *ob,
                    const gsl_vector *lambda,
                    const gsl_matrix *Q,
                    const gsl_vector *c,
                    const double sdR,
                    const double sdZ,
                    gsl_vector *uvw,
                    gsl_vector *ytmp) {
    
    // Shift and rotate to the eigenspace
    // Qt (y - c) ---> (u, v, w)
    gsl_vector_memcpy(ytmp, c); // ytmp <- c
    gsl_vector_scale(ytmp, -1); // ytmp <- -o
    gsl_vector_add(ytmp, ob);   // ytmp <- ytmp + (an observation)
                                // [ytmp <- (-c) + (an observation)]
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
  
  double
  datasetLogLik(Observations *obs,
                const double L,
                const gsl_vector *lambda,
                const gsl_matrix *Q,
                const gsl_vector *c,
                const double sdR,
                const double sdZ) {
    size_t N = obs->N;
    size_t D = obs->D;
    
    double logLik = 0;
    gsl_vector *uvw    = gsl_vector_alloc(D);
    gsl_vector *ytmp   = gsl_vector_alloc(D);
    gsl_vector *ob     = gsl_vector_alloc(D);

    for (size_t n = 0; n < N; n++) { 
      obs->writeObservation(n, ob);     
      logLik +=
        _obserationLogLik(ob, lambda, Q, c, sdR, sdZ, uvw, ytmp);
    }
    
    gsl_vector_free(ytmp);
    gsl_vector_free(uvw);
    
    return logLik;
  }
}
