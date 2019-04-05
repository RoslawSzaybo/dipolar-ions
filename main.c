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
#include <math.h>
#include "find-spectrum.h"
#include "hamiltonian.h"


/*
MgH+
q = 1e
dipole = 3D (Journal of Physics B: AMO, vol 42, no 15 (2009) M. Aymar et al.)
B = 180 GHz (9K) (New Journal of Physics 11 (2009) 055026)
m = 25.3 u
*/


/* Main program */
int main(int argc, char *argv[]) {
        /* Locals */
        int n, info, lwork;
        int n1, n3, n5, j1, j2;
        float mass, charge, dipole, B, omega_rho, omega_z;
        int is_input_OK = 0;
        if(argc != 12)
        {
                printf( "This is how you can execute the program:\n" );
                printf( "\t$ %s n1 n3 n5 j1 j2 "
                "mass[u] charge[e] dipole[D] B[MHz] omega_rho[MHz] omega_z[MHz]\n", argv[0]);
                return 0;
        }
        else
        {
                n1 = atoi(argv[1]);
                n3 = atoi(argv[2]);
                n5 = atoi(argv[3]);
                j1 = atoi(argv[4]);
                j2 = atoi(argv[5]);
                mass = atof(argv[6]);
                charge = atof(argv[7]);
                dipole = atof(argv[8]);
                B = atof(argv[9]);
                omega_rho = atof(argv[10]);
                omega_z = atof(argv[11]);
                printf("# Parameters of the model\n");
                printf("#  mass:\t\t%10.2f u\n", mass);
                printf("#  charge:\t%10.2f e\n", charge);
                printf("#  dipole:\t%10.2f D\n", dipole);
                printf("#  B:\t\t%10.2f MHz\n", B);
                printf("#  omega_rho:\t%10.2f MHz\n", omega_rho);
                printf("#  omega_z:\t%10.2f MHz\n\n", omega_z);
        }
        n = n1*n3*n5*(j1*j1+2*j1+1)*(j2*j2+2*j2+1);
        // input test
        is_input_OK = test_input(n1, n3, n5, j1, j2, omega_rho, omega_z);
        if( !is_input_OK )
        {
                printf(" Basis truncation |%d,%d,%d,%d,%d>, or"
                " trap frequency is incorrect. \n",n1,n3,n5,j1,j2);
                exit(0);
        }
        basis b = {n1, n3, n5, j1, j2};
        float omega_1 = sqrt(3.f)*omega_z;
        float omega_3 = sqrt(omega_rho*omega_rho - omega_z*omega_z);
        parameters pars = {mass, charge, dipole, B, omega_1, omega_3};
        printf("# omega_1:\t%10.2f MHz\n", pars.omega_1);
        printf("# omega_3/5:\t%10.2f MHz\n", pars.omega_3);
        printf("# Basis truncation:\t|%d,%d,%d,%d,%d>\n",n1,n3,n5,j1,j2);
        printf("# Basis size:\t\t%d  \n\n",n);
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
        // Construct the Hamiltonian matrix
        construct_Hamiltonian(a, b, pars);
        /* Print matrix to be diagonalised */
        //print_matrix( "# Hamiltonian", n, n, a, n );
        /* Executable statements */
        printf( "# Results of the diagonalisation\n" );
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
        print_rmatrix( "# 100 smallest eigenvalues", 1, 100, w, 1 );
		/* Print all eigenvalues */
        // print_rmatrix( "# 100 smallest eigenvalues", 1, n, w, 1 );
        /* Print ground state */
        print_matrix( "# Ground state (stored columnwise)", n, 1, a, n );
        /* Print all eigenvectors */
        //print_matrix( "Eigenvectors (stored columnwise)", n, n, a, n );
        /* Free workspace */
        free( (void*)work );
        free( w );
        free( rwork );
        free( a );
        exit( 0 );
} /* End of diagonalisation. */