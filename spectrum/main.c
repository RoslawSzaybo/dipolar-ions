/*
This is a file which used to serve as an example on how to use the
cheev function. 
https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/cheev_ex.c.htm

Now it diagonalises the Hamiltonian matrix of two dipolar ions. 
The Hamiltonian is turend into a matrix by computing it's expectation 
values in a basis of (harmonic oscilatior)x(spherical hamonics)^{x2}.
Both bases are truncated. As a resul the truncated version of the Hamiltonian
is represented as a hermitian matirx. That matrix can be diagonalised.
The eigenvectors corresponding to the lowest eigenvalues constitute 
an approximation to the ground state and the low excitations of the system.

This method of solving the Schr\"odinger equation in literature is often 
reffered to as an exact diagonalisation.
*/

#include <stdlib.h>
#include <stdio.h>

#include "find-spectrum.h"
//#include "hamiltonian.h"
#include "input.h"

/* Main program */
int main(int argc, char *argv[]) {
    /* Work with input */
    check_argc(argc, argv);

    check_physical_parameters(argv);
    parameters pars = get_system_parameters(argv);

    basis b = get_basis_truncation(argv);
    check_basis_truncation(b);
    int n = get_basis_size(b);

    print_system_parameters(pars);
    print_basis_details(b);
    print_active_terms_of_Hamiltonian(pars);

    /* Executable statements */
    /* allocate memory for the matrix */
    long matrix_size = sizeof(dcomplex);
    matrix_size *= n;
    matrix_size *= n;
    dcomplex *a = my_malloc((size_t)matrix_size, "a");
    /* Define the Hamiltonian matrix */
    // construct_Hamiltonian(a, b, pars);
    int p,r;
    for (p=0; p < n; p++)
        for (r=0; r<n; r++)
        {
            a[p*n+r].re = 0.;
            a[p*n+r].im = 0.;
        }
    // upper triangle
    for (p=0; p < n; p++)
        for (r=p; r<n; r++)
        {
            a[p*n+r].re = p/100.;
            a[p*n+r].im = r/100.;
        }
    /* Print the Hamiltonian Matrix */
    /* WARNING: Usually you don't want to do it! */
    /* WARNING: It's mostly for testing purposes */
    print_matrix( "# Hamiltonian", n, n, a, n );

    /* Query and allocate the optimal workspace */
    double *w = my_malloc(sizeof(double)*n, "w");
    dcomplex wkopt;
    /* rwork dimension should be at least max(1,3*n-2) */
    double *rwork = my_malloc(sizeof(double)*(3*n-2), "rwork");
    int lwork = -1, info;
    zheev( "Vectors", "Lower", &n, a, &n, w, &wkopt, &lwork, rwork, &info );
    lwork = (int)wkopt.re;
    dcomplex *work = my_malloc( lwork*sizeof(dcomplex), "work" );

    /* Solve eigenproblem */
    zheev( "Vectors", "Lower", &n, a, &n, w, work, &lwork, rwork, &info );
    /* Check for convergence */
    if( info > 0 ) {
            printf( "The algorithm failed to compute eigenvalues.\n" );
            exit( 1 );
    }
    /* End of diagonalisation. */

    /* Free workspace - part I */
    // It's larger than work_int so it should 'make space' 
    // for work_int, but you know.
    free( rwork ); 

    /* Generate the output */
    /* Print eigenvalues and eigenvectors of the lower spectrum */
    printf( "# Results of the diagonalisation\n" );
    print_rmatrix( "# 100 smallest eigenvalues", 1, (n<100)?n:100 , w, 1 );
    printf("# sanity check\n");
    print_eigenvalues( 0, (n<100)?n:100, w);

    /* Print states more explicitly */
    // Space required for sorting
    int *work_int = my_malloc(sizeof(int)*n, "work_int");
    int how_many_states = 25;
    how_many_states = ( n<how_many_states ) ? n : how_many_states;
    sort_print_lower_spectrum( a, n, b, how_many_states, work, work_int );

    /* Print eigenvalues and eigenvectors of the lower spectrum */
    // this will also print the last withouth j excited, for 
    // sanity check
    int first = b.n1*b.n3*b.n5;
    int last = first + 50;
    first = (n < first) ? n : first;
    last = (n < last) ? n : last;
    printf( " # The upper part of the spectrum\n" );
    printf( " # first: %d\t last: %d\n", first, last );
    print_eigenvalues( first, last, w );
    sort_print_upper_spectrum( a, n, b, first, last, work, work_int );

    /* Free workspace - part II */
    free( a );
    free( w );
    free( work );
    free( work_int );
    exit( 0 );
} 