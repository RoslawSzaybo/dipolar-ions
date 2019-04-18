/*
Customised C-version of a vector. 

Allows to store states which are in a superposition 
of a few basis vectors.

A very good example of a standard vector implementation can be found here:
https://gist.github.com/EmilHernvall/953968/0fef1b1f826a8c3d8cfb74b2915f17d2944ec1d0
*/

#ifndef state_h
#define state_h

#include "find-spectrum.h"


typedef struct {
    int n1, n3, n5, j1, j2;
} basis;

typedef struct {
    float mass, charge, dipole, B, omega_1, omega_3;
} parameters;

typedef struct {
    int n1, n3, n5, j1, m1, j2, m2;
} versor;

typedef struct { 
    int length; 
    int cap;
    versor* kets;
    fcomplex* amplitudes;
} state;

void state_init(state* psi);
void state_free(state* psi);
int state_cnt(state* psi);
void state_add(state* psi, versor ket, fcomplex amplitude);
versor state_get_versor(state* psi, int idx);
fcomplex state_get_amplitude(state* psi, int idx);
void state_times_float(state * input, const float factor);
void state_times_i_float(state * input, const float factor);
void state_add_state(state * sum, const state * compound);

#endif // state_h