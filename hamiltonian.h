#ifndef hamiltonian_h
#define hamiltonian_h

#include "find-spectrum.h"
#include "state.h"


int jm_jump(int j, int m);
/* 
We work in a product basis
$$
\ket{\psi} = \ket{n1}\ket{n3}\ket{n5}\ket{j1,m1}\ket{j2,m2}.
$$
This basis is sorted, so that to each versor there corresponds 
an index $idx$. The order is available only because the basis 
is truncated i.e. $n1 < N1, \ldots, j1 <= J1, \ldots, m2 \le M2$.
This function returns the index of a given set of quantum number 
and a truncation.
*/
int get_index_from_versor(versor psi, const basis b);
versor get_versor_from_index(int idx, const basis b);

/*
* Returns $j$ or $m$ quantum number of a versor $\ket{j,m}$
* which is in the positin `idx` on a list of all
* $\ket{j,m}$. The list is sorted so that $j$ 
* is non-decreasing, and within the same $j$,
* $m$ is increasing: $\ket{0,0}, \ket{1,-1},
* \ket{1,0}, \ket{1,1}, \ket{2,-2}, \ldots$.
*/
int getj(int idx);
int getm(int idx);


void print_versor(const versor psi, const basis b);
void test_idx_to_versor_translation();

/*
Upon acting with state $\ket{\psi}$ from the left 
on the Hamiltonian one gets a superposition of many bras.
Those bras are stored in the output state.
*/

state bra_H(state* psi);
void test_bra_H();
void construct_Hamiltonian(fcomplex* a, basis b);

#endif // hamiltonian_h