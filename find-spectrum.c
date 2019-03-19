#include <stdlib.h>
#include <stdio.h>
#include "find-spectrum.h"

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, fcomplex* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", a[i+j*lda].re, a[i+j*lda].im );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, float* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
                printf( "\n" );
        }
}

fcomplex fcomplex_multiply(const fcomplex* a, const fcomplex* b)
{
    float re=0.f, im=0.f;

    re = a->re * b->re - a->im * b->im;
    im = a->re * b->im  + a->im * b->re;

    return (fcomplex){re, im};
}

fcomplex fcomplex_times_float(const fcomplex* a, float b)
{
        float re = a->re * b;
        float im = a->im * b;

        return (fcomplex){re, im};
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