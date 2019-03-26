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
        int n1, n3, n5, j1, j2;
        int is_input_OK = 0;
        if(argc != 6)
        {
                printf( "This is how you can execute the program:\n" );
                printf( "\t$ %s n1 n3 n5 j1 j2\n", argv[0]);
                return 0;
        }
        else
        {
                n1 = atoi(argv[1]);
                n3 = atoi(argv[2]);
                n5 = atoi(argv[3]);
                j1 = atoi(argv[4]);
                j2 = atoi(argv[5]);

        }
        n = n1*n3*n5*(j1*j1+2*j1+1)*(j2*j2+2*j2+1);
        // input test
        is_input_OK = test_input(n1, n3, n5, j1, j2);
        if( !is_input_OK )
        {
                printf("Input basis |%d,%d,%d,%d,%d> is incorrect. \n",n1,n3,n5,j1,j2);
                exit(0);
        }
        basis b = {n1, n3, n5, j1, j2};
        //test_bra_H();
        fcomplex wkopt;
        fcomplex* work;
        /* Local arrays */
        /* rwork dimension should be at least max(1,3*n-2) */
        float *w = (float*)malloc(sizeof(float)*n);
        if(w == NULL)
        {
                printf("Not enough memory to allocate w.\n");
                return 1;
        }
        float *rwork = (float*)malloc(sizeof(float)*(3*n-2));
        if(rwork == NULL)
        {
                printf("Not enough memory to allocate rwork.\n");
                return 1;
        }

        fcomplex *a = (fcomplex*)malloc(sizeof(fcomplex)*n*n);
        if(a == NULL)
        {
                printf("The size of a matirx is %dx%d=%d\n",n,n,n*n);
                printf("Not enough memory to allocate a.\n");
                return 1;
        }
        // generate some particular matrix for testing
        construct_Hamiltonian(a, b);
        /* Print matrix to be diagonalised */
        print_matrix( "Hamiltonian", n, n, a, n );
        /* Executable statements */
        //printf( " Results of the diagonalisation\n" );
        /* Query and allocate the optimal workspace */
        //lwork = -1;
        //cheev( "Vectors", "Lower", &n, a, &n, w, &wkopt, &lwork, rwork, &info );
        //lwork = (int)wkopt.re;
        //work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
        /* Solve eigenproblem */
        //cheev( "Vectors", "Lower", &n, a, &n, w, work, &lwork, rwork, &info );
        /* Check for convergence */
        /*
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        */
        /* Print eigenvalues */
        //print_rmatrix( "Five smallest eigenvalues:", 1, 5, w, 1 );
        //print_rmatrix( "Eigenvalues", 1, n, w, 1 );
        /* Print eigenvectors */
        //print_matrix( "Eigenvectors (stored columnwise)", n, n, a, n );
        /* Free workspace */
        //free( (void*)work );
        free( w );
        free( rwork );
        free( a );
        exit( 0 );
} /* End of diagonalisation. */