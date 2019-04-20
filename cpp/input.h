#ifndef input_h
#define input_h

typedef struct {
    float mass, charge, dipole, B, omega_1, omega_3;
} parameters;

int test_input(int n1, int n3, int n5, int j1, int j2, 
                float omega_rho, float omega_z);

#endif // input_h