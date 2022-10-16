#pragma once
#include "basedefs.h"

/* From wikipedia: https://en.wikipedia.org/wiki/Xorshift */
/* xorwow is used as default PRNG in CUDA Toolkit. */
/* The state array must be initialized to not be all zero in the first four words */
struct xorwow_state {
    unsigned int x[5];
    unsigned int counter;
};

unsigned int xorwow(xorwow_state * state);

class XorwowRNG : public Noncopyable {
    /* directly copy a RNG can be dangerous */
protected:
    xorwow_state state;
public:
    XorwowRNG();
    XorwowRNG(const int& a, const int& b);
    void init();
    void init(const int& a, const int& b);
    REAL uniform();
    REAL uniform(const REAL& min, const REAL& max);
    REAL gaussian();
    REAL gaussian(const REAL& mu, const REAL& sigma);
    int uniform_int(const int& min, const int& max); /* generate a random number between [min, max) */
};

/* generate normally distributed random variable using xorwow method */
REAL uniform();
REAL uniform(REAL min, REAL max);

/* generate a random number between [min, max) */
int uniform_int(const int& min, const int& max); 

/* generate normally distributed random variable using the Box-Muller method */
REAL normal();
REAL normal(REAL mean, REAL std);

