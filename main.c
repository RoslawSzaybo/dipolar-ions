/*
this is a file which used to serve as an example on how to use the
cheev function. 
https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/cheev_ex.c.htm

Now it diagonalises the Hamiltonian of two dipolar ions. 
The Hamiltonian is turend into a matrix by computing it's expectation 
values in a basis of (harmonic oscilatior)x(spherical hamonics)^2.
Both bases are truncated. As a resul a the truncated version of the Hamiltonian
is represented as a hermitian matirx. That matrix can be diagonalised.
The eigenvectros corresponding to the lowest eigenvalues constitute 
an approximation to the ground state of the system.
*/

#include <stdio.h>
#include <stdlib.h>
#include "find-spectrum.h"
#include "hamiltonian.h"

/* Main program */
int main(int argc, char *argv[]) {
        /* Locals */
        int n, info, lwork;
        if(argc != 2)
        {
                printf(
                        "Program requires exactly one argument" 
                        " – the matrix size.\n");
                return 0;
        }
        else
        {
               n = atoi(argv[1]);
        }
        test_idx_to_state_translation();
        fcomplex wkopt;
        fcomplex* work;
        /* Local arrays */
        /* rwork dimension should be at least max(1,3*n-2) */
        float *w = (float*)malloc(sizeof(float)*n),
                *rwork = (float*)malloc(sizeof(float)*(3*n-2));
        if(w == NULL || rwork == NULL)
        {
                printf("Not enough memory to allocate w or rwork.\n");
                return 1;
        }

        fcomplex *a = (fcomplex*)malloc(sizeof(fcomplex)*n*n);
        if(a == NULL)
        {
                printf("Not enough memory to allocate a.\n");
                return 1;
        }
        // generate some particular matrix for testing
        for(int row=0; row<n; row++)
        {
                for(int column=0; column<n; column++)
                {
                        fcomplex entry = {0.0f, 0.0f};
                        if(row <= column)
                        {
                                entry.re = (float)(row*n+column);
                                entry.im = (float)(column*n+row);
                        }
                        a[row*n + column] = entry;
                }
        }
        /* Print matrix to be diagonalised */
        //print_matrix( "Matrix", n, n, a, n );
        /* Executable statements */
        printf( " CHEEV Example Program Results\n" );
        /* Query and allocate the optimal workspace */
        lwork = -1;
        cheev( "Vectors", "Lower", &n, a, &n, w, &wkopt, &lwork, rwork, &info );
        lwork = (int)wkopt.re;
        work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
        /* Solve eigenproblem */
        cheev( "Vectors", "Lower", &n, a, &n, w, work, &lwork, rwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        /* Print eigenvalues */
		printf("\n\nFull success in diagonalising %dx%d matirx\n",n,n);
		printf("Five largerst eigenvalues:\n");
        print_rmatrix( "Eigenvalues", 1, 5, w, 1 );
        //print_rmatrix( "Eigenvalues", 1, n, w, 1 );
        /* Print eigenvectors */
        //print_matrix( "Eigenvectors (stored columnwise)", n, n, a, n );
        /* Free workspace */
        free( (void*)work );
        free( w );
        free( rwork );
        free( a );
        exit( 0 );
} /* End of CHEEV Example */