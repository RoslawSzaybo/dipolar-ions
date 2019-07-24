#ifndef dcomplex_h
#define dcomplex_h

/* Complex datatype */
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;

dcomplex dcomplex_multiply(const dcomplex* a, const dcomplex* b);
dcomplex dcomplex_times_double(const dcomplex* a, double b);
dcomplex dcomplex_times_i_double(const dcomplex* a, double b);
double dcomplex_amplitude_sqr(const dcomplex* a);
void test_dcomplex_multiply();

#endif // dcomplex_h