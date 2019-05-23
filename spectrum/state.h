/*
Customised C-version of a vector.

In C++ that would be something linke vector<versor>

Allows to store states which are superpositions (linear combinations)
of states that can be represented in the product basis (struct versor).

A very good example of a standard vector implementation can be found here:
https://gist.github.com/EmilHernvall/953968/0fef1b1f826a8c3d8cfb74b2915f17d2944ec1d0
*/

#ifndef state_h
#define state_h

#include "fcomplex.h"
#include "versor.h"

typedef struct { 
    int length; 
    int cap;
    versor* kets;
    fcomplex* amplitudes;
} state;

void state_init(state* psi);
void state_free(state* psi);
int state_cnt(state* psi);
int state_versor_location(state *psi, const versor *ket);
int state_contains_versor(state *psi, const versor *ket);
void state_add(state* psi, versor ket, fcomplex amplitude);
versor state_get_versor(state* psi, int idx);
fcomplex state_get_amplitude(state* psi, int idx);
void state_times_float(state * input, const float factor);
void state_times_i_float(state * input, const float factor);
void state_add_state(state * sum, const state * compound);
void state_clean(state* psi, basis b);

#endif // state_h