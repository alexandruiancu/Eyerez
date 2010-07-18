#include "sequentially.hpp"

#include <vector>

Sequentially::Sequentially (samp_vec samplers_in) {
  samplers = samplers_in;
  S = samplers.size();
}

Sequentially::~Sequentially () {
  // will call of the sampler's deconstructors as well.
  delete &samplers;
}

Sequentially *
Sequentially::clone () {
  samp_vec samplers_new = samplers;
  Sequentially *samp = new Sequentially(samplers_new);
  return samp;
}

void
Sequentially::jump (State *x) {
  samp_vec::iterator it;  
  for (it = samplers.begin() ; it < samplers.end(); it++) {
    (*it)->jump(x);
  }
}

void
Sequentially::freeze () {
  samp_vec::iterator it;
  for (it = samplers.begin() ; it < samplers.end(); it++) {
    (*it)->freeze();
  }
}

bool
Sequentially::isFrozen () {
  bool out = true;
  samp_vec::iterator it;
  for (it = samplers.begin() ; it < samplers.end(); it++) {
    out &= (*it)->isFrozen();
  }
  return out;
}

void
Sequentially::initState (State *st) {
  samp_vec::iterator it;
  for (it = samplers.begin() ; it < samplers.end(); it++) {
    (*it)->initState(st);
  }
}
