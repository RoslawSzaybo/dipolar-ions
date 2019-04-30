#include <stdlib.h>
#include <stdio.h>

#include "find-spectrum.h"


/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, fcomplex* a, int lda ) {
        int i, j;
        printf( "\n%s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", a[i+j*lda].re, a[i+j*lda].im );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, float* a, int lda ) {
        int i, j;
        printf( "\n%s\n", desc );
        for ( i = 0; i < m; i++ ) {
                for ( j = 0; j < n; j++ ) printf( " %11.9f", a[i+j*lda] );
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
void sort_print_lower_spectrum(fcomplex* a, int n, basis b, int m, fcomplex* work, int* work_int)
{ 
        printf("\n# Ground state and the lowest excitations:\n");
        // if the desired number of vectors, m, is greater than
        // the their total number, n, then print all n vectors.
        int bound = (m<n) ? m : n ;
        for (int k=0; k < bound; k++)
                sort_print_eigenvector_summary( a, work, work_int, n, b, k);
}

void sort_print_eigenvector_summary( fcomplex* a, fcomplex* work, int* work_int, int n, basis b, int m )
{
        // Preselection, the vectors which are of no importance will be neglected;
        const float lower_limit = 1.0e-7f;
        int i, length=0;
        fcomplex amp;
        float abs_amp2;
        for ( i = 0; i < n; i++ )
        {
                amp = a[ m+i*n ];
                abs_amp2 = fcomplex_amplitude_sqr( &amp );
                if ( abs_amp2 > lower_limit ) 
                {
                        work[length].re = amp.re;
                        work[length].im = amp.im;
                        work_int[length] = i;
                        length++;
                }
        }

        // sorting
        sort_fcomplex(work, work_int, length);

        // printing
        printf( " |%2d> = ", m );
        versor ket;
        for ( i = 0; i < length; i++ )
        {
                if(i!=0)
                        printf("        ");
                amp = work[i];
                printf( "(%6.3f,%6.3f)", amp.re, amp.im );
                ket =  get_versor_from_index( work_int[i], b );
                show_versor( ket );
                printf( " +\n" );
        }
        printf( "         ...\n" );
}

void sort_fcomplex(fcomplex* data, int* indices, int length)
{
        int sorted;
        int i;
        float amp0;
        float amp1;
        fcomplex temp_fc;
        int temp_int;
        while( 1 ) 
        {
                sorted = 1;
                amp0  =  fcomplex_amplitude_sqr( data );
                for (i=1; i<length; i++)
                {
                        // amp0 is the amplitude of data[i-1]
                        // amp1 is the amplitude of data[i]
                        amp1 = fcomplex_amplitude_sqr( data+i );
                        if (amp1 > amp0)
                        {
                                // swap data
                                temp_fc = data[i];
                                data[i] = data[i-1]; 
                                data[i-1] = temp_fc;
                                // swap int
                                temp_int = indices[i];
                                indices[i] = indices[i-1];
                                indices[i-1] = temp_int;
                                // raise the flag
                                sorted = 0;
                                // prepare amp0 for the next run of the loop
                                // amp0 = amp0;
                        }
                        else
                        {
                                // prepare amp0 for the next run of the loop
                                amp0 = amp1;
                        }
                }
                if (sorted == 1)
                {
                        break;
                }
        }
}
