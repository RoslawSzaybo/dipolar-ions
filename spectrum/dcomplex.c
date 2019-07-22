#include <stdio.h>
#include "dcomplex.h"

/* 
* Multiply two complex numbers stored as dcomplex 
*/
dcomplex dcomplex_multiply(const dcomplex* a, const dcomplex* b)
{
    double re=0., im=0.;

    re = a->re * b->re - a->im * b->im;
    im = a->re * b->im  + a->im * b->re;

    return (dcomplex){re, im};
}

/* 
* Multiply complex number `a` by a real number `b` 
* return the result
* No change in the input variable `a`.
*/
dcomplex dcomplex_times_double(const dcomplex* a, double b)
{
    double re = a->re * b;
    double im = a->im * b;

    return (dcomplex){re, im};
}


/*
* Multiply a complex number `a` by a pruely imaginary number `ib` 
* and return the result. No change in the input variable `a`. 
* $$
* a * ib, \quad a \in \mathbb{C}, b\in \mathbb{R}.
* $$ 
*/
dcomplex dcomplex_times_i_double(const dcomplex* a, double b)
{
    // (re + i im)* i b = -im*b + i re*b
    double re = a->re;
    double im = a->im;

    return (dcomplex){-im*b, re*b};
}

/* 
* Return a square of the amplitude of a complex number `a`
* $$
* |a|^2 = \overline{a}a = re(a)^2 + im(a)^2.
* $$
*/
double dcomplex_amplitude_sqr(const dcomplex* a)
{
    return a->im*a->im + a->re*a->re;
}

void test_dcomplex_multiply()
{
    printf("Test dcomplex_multiply\n");
    dcomplex fa = {2.0, 0.0}, fb = {1.0, 1.0};
    dcomplex fc = dcomplex_multiply(&fa, &fb);
    printf("fa = %3.2f+i%3.2f\n", fa.re, fa.im);
    printf("fb = %3.2f+i%3.2f\n", fb.re, fb.im);
    printf("fc = %3.2f+i%3.2f\n", fc.re, fc.im);
}