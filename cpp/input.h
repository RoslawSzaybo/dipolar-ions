#ifndef input_h
#define input_h

#include <stdlib.h>
#include "versor.h"

typedef struct {
    float mass, charge, dipole, B, omega_1, omega_3;
} parameters;

void check_argc(int argc, char *argv[]);
void check_physical_parameters(char *argv[]);
void check_basis_truncation(basis b);

int get_basis_size(basis b);
basis get_basis_truncation(char *argv[]);
parameters get_system_parameters(char *argv[]);

void print_input(basis b, parameters p, char *argv[]);

void *my_malloc(size_t size, const char *name);

#endif // input_h