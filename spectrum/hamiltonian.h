#ifndef hamiltonian_h
#define hamiltonian_h

#include "input.h" 
#include "versor.h"
#include "find-spectrum.h"
#include "state.h"

/* prints |n1,n3,n5,j1,m1,j2,m2> and state index */
void print_versor(const versor psi, const basis b);
void test_idx_to_versor_translation();

int valid_versor(versor psi, basis b);

/*
Upon acting with state $\ket{\psi}$ from the left 
on the Hamiltonian one gets a superposition of many bras.
Those bras are stored in the output state.

Here are listed operators which are already implemented.
*/

void apply_harmonic_oscillator(state * input, state * output, 
                                const parameters pars);
void apply_rotational_kinetic_energy(state * input, state * output, 
                                        const parameters pars);
void apply_a1_plus_a1dagger(state *input, state *output, 
                                    const parameters pars);
void apply_a3_plus_a3dagger(state *input, state *output,
                                    const parameters pars);
void apply_a5_plus_a5dagger(state *input, state *output,
                                    const parameters pars);
void apply_dz1(state * input, state * output, const parameters pars);
void apply_dz2(state * input, state * output, const parameters pars);
void apply_dy1(state * input, state * output, const parameters pars);
void apply_charge_dipole_zero(state *input, state *output,
                                    const parameters pars);

state bra_H(state* psi, const parameters pars);

void construct_Hamiltonian(dcomplex* a, const basis b, const parameters pars);

#endif // hamiltonian_h