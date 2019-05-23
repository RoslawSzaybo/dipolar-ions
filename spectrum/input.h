#ifndef input_h
#define input_h

#include <stdlib.h>
#include "versor.h"

// this structure holds information on which terms
// of the system Hamiltionan are active 
// in the selected run of the simulation
typedef struct {
    int normal_modes, 
    T_rot,
    Vqd_zeroth,
    Vqd_first, 
    Vdd_zeroth;
} hamiltonian;

typedef struct {
    float mass, charge, dipole, B, omega_1, omega_3;
    hamiltonian active_terms;
} parameters;

void check_argc(int argc, char *argv[]);
void check_physical_parameters(char *argv[]);
void check_basis_truncation(basis b);

int get_basis_size(basis b);
basis get_basis_truncation(char *argv[]);
parameters get_system_parameters(char *argv[]);

void print_input(basis b, parameters p, char *argv[]);
void print_active_terms_of_Hamiltonian(parameters p);

void *my_malloc(size_t size, const char *name);

#endif // input_h