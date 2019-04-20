#include <stdlib.h>
#include <stdio.h>

#include "find-spectrum.h"


/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, fcomplex* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.3f,%6.3f)", a[i+j*lda].re, a[i+j*lda].im );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, float* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for ( i = 0; i < m; i++ ) {
                for ( j = 0; j < n; j++ ) printf( " %6.3f", a[i+j*lda] );
                printf( "\n" );
        }
}


/* Print the most prominent amplitudes with the corresponding kets, 
for a selected eigenstate. 
Warning: The amplitudes are not sorted.
input:
a - matrix with the diagonalisation result
n - basis size
m - which eigenstate do you want to print 
        (eigenstates are orderd with respect to their eigenvalues, in an 
        ascending order i.e. 0th eigenvector has the smallest eigenvalue)
b - descriptor of the basis which was used in the computations */
void print_eigenvector_summary( fcomplex* a, int n, basis b, int m)
{
        const float lower_limit = 1.0e-7f;
        printf( " |%d> = ", m );
        int i;
        fcomplex amp;
        float abs_amp2;
        versor ket;
        for ( i = 0; i < n; i++ )
        {
                amp = a[m+i*n];
                abs_amp2 = fcomplex_amplitude_sqr(&amp);
                if ( abs_amp2 > lower_limit) 
                {
                        printf( "(%6.3f,%6.3f)", amp.re, amp.im );
                        ket =  get_versor_from_index(i, b);
                        show_versor(ket, b);
                        printf(" + ");
                }
        }
        printf("...\n");
}

void print_lower_spectrum(fcomplex* a, int n, basis b, int m)
{ 
        printf("\n Ground state and the lowest excitations:\n");
        // if the desired number of vectors, m, is greater than
        // the their total number, n, then print all n vectors.
        int bound = (m<n) ? m : n ;
        for (int k=0; k < bound; k++)
                print_eigenvector_summary( a, n, b, k);
}