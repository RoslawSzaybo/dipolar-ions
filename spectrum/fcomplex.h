#ifndef fcomplex_h
#define fcomplex_h

/* Complex datatype */
struct _fcomplex { float re, im; };
typedef struct _fcomplex fcomplex;

fcomplex fcomplex_multiply(const fcomplex* a, const fcomplex* b);
fcomplex fcomplex_times_float(const fcomplex* a, float b);
fcomplex fcomplex_times_i_float(const fcomplex* a, float b);
float fcomplex_amplitude_sqr(const fcomplex* a);
void test_fcomplex_multiply();

#endif // fcomplex_h