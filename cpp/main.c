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

#include "find-spectrum.h"
#include "hamiltonian.h"
#include "input.h"

#include <stdlib.h>
#include <stdio.h>

/*
MgH+
q = 1e
dipole = 3D (Journal of Physics B: AMO, vol 42, no 15 (2009) M. Aymar et al.)
B = 180 GHz (9K) (New Journal of Physics 11 (2009) 055026)
m = 25.3 u
*/

/* Main program */
int main(int argc, char *argv[]) {
        /* Work with input */
        check_argc(argc, argv);

        check_physical_parameters(argv);
        parameters pars = get_system_parameters(argv);

        basis b = get_basis_truncation(argv);
        check_basis_truncation(b);
        int n = get_basis_size(b);

        print_input(b, pars, argv);

        //test_bra_H();

        /* Executable statements */
        /* allocate memory for the matrix */
        long matrix_size = sizeof(fcomplex);
        matrix_size *= n;
        matrix_size *= n;
        fcomplex *a = my_malloc((size_t)matrix_size, "a");
        /* Define the Hamiltonian matrix */
        construct_Hamiltonian(a, b, pars);
        /* WARTNING: Print the Hamiltonian Matrix */
        /* WARTNING: Usually you don't want to do it! */
        /* WARTNING: It's mostly for testing purposese */
        // print_matrix( "# Hamiltonian", n, n, a, n );


        /* Query and allocate the optimal workspace */
        float *w = my_malloc(sizeof(float)*n, "w");
        fcomplex wkopt;
        /* rwork dimension should be at least max(1,3*n-2) */
        float *rwork = my_malloc(sizeof(float)*(3*n-2), "rwork");
        int lwork = -1, info;
        cheev( "Vectors", "Lower", &n, a, &n, w, &wkopt, &lwork, rwork, &info );
        lwork = (int)wkopt.re;
        fcomplex *work = my_malloc( lwork*sizeof(fcomplex), "work" );

        /* Solve eigenproblem */
        cheev( "Vectors", "Lower", &n, a, &n, w, work, &lwork, rwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        /* Free workspace - part I */
		// It's larger than work_int so it should 'make space' 
		// for work_int, but you know.
        free( rwork ); 

		/* Generate the output */
        printf( "# Results of the diagonalisation\n" );
        /* Print eigenvalues */
        print_rmatrix( "# 100 smallest eigenvalues", 1, (n<100)?n:100 , w, 1 );
        /* Print all eigenvalues */
        // print_rmatrix( "# 100 smallest eigenvalues", 1, n, w, 1 );

        /* Print states more explicitly */
        // Space required in  sorting
        int *work_int = my_malloc(sizeof(int)*n, "work_int");
        sort_print_lower_spectrum(a, n, b, 25, work, work_int);

        /* Free workspace - part II */
        free( a );
        free( w );
        free( work );
        free( work_int );
        exit( 0 );
} /* End of diagonalisation. */
