#include <stdio.h>
#include "fcomplex.h"

/* Multiply two complex numbers stored as fcomplex */
fcomplex fcomplex_multiply(const fcomplex* a, const fcomplex* b)
{
    float re=0.f, im=0.f;

    re = a->re * b->re - a->im * b->im;
    im = a->re * b->im  + a->im * b->re;

    return (fcomplex){re, im};
}

/* Multiply complex number `a` by a float `b` and return the result.
Doesn't change the input variable `a`.*/
fcomplex fcomplex_times_float(const fcomplex* a, float b)
{
        float re = a->re * b;
        float im = a->im * b;

        return (fcomplex){re, im};
}


/* Multiply complex number `a` by an imaginary float `b` and return the result.
Doesn't change the input variable `a`. It would correspond to 
$$
a * ib, \quad a \in \mathbb{C}, b\in \mathbb{R}.
$$ */
fcomplex fcomplex_times_i_float(const fcomplex* a, float b)
{
        // (re + i im)* i b = -im*b + i re*b
        float re = a->re;
        float im = a->im;

        return (fcomplex){-im*b, re*b};
}

/* Return square of an amplitude of the complex number  `a`
$$
|a|^2 = \overline{a}a = re(a)^2 + im(a)^2.
$$*/
float fcomplex_amplitude_sqr(const fcomplex* a)
{
        return a->im*a->im + a->re*a->re;
}

void test_fcomplex_multiply()
{
        printf("Test fcomplex_multiply\n");
        fcomplex fa = {2.0f, 0.0f}, fb = {1.0f, 1.0f};
        fcomplex fc = fcomplex_multiply(&fa, &fb);
        printf("fa = %3.2f+i%3.2f\n", fa.re, fa.im);
        printf("fb = %3.2f+i%3.2f\n", fb.re, fb.im);
        printf("fc = %3.2f+i%3.2f\n", fc.re, fc.im);
}