/*
this is a file which served as an example on how to use the
cheev function. 

I have barely edited it to test how large matrices can be 
diagonalised on my machine.
*/

#include <stdlib.h>
#include <stdio.h>

/* Complex datatype */
struct _fcomplex { float re, im; };
typedef struct _fcomplex fcomplex;

/* CHEEV prototype */
extern void cheev( char* jobz, char* uplo, int* n, fcomplex* a, int* lda,
                float* w, fcomplex* work, int* lwork, float* rwork, int* info );
/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, fcomplex* a, int lda );
extern void print_rmatrix( char* desc, int m, int n, float* a, int lda );

/* Parameters */
#define N 581
#define LDA N

/* Main program */
int main() {
        /* Locals */
        int n = N, lda = LDA, info, lwork;
        fcomplex wkopt;
        fcomplex* work;
        /* Local arrays */
        /* rwork dimension should be at least max(1,3*n-2) */
        float w[N], rwork[3*N-2];
        fcomplex a[LDA*LDA];
        // generate some particular matrix
        for(int row=0; row< N; row++)
        {
                for(int column=0; column<N;column++)
                {
                        fcomplex entry = {0.0f, 0.0f};
                        if(row <= column)
                        {
                                entry.re = (float)(row*N+column);
                                entry.im = (float)(column*N+row);
                        }
                        a[row*N + column] = entry;
                }
        }
        /* Print matrix to be diagonalised */
        //print_matrix( "Matrix", n, n, a, lda );
        /* Executable statements */
        printf( " CHEEV Example Program Results\n" );
        /* Query and allocate the optimal workspace */
        lwork = -1;
        cheev( "Vectors", "Lower", &n, a, &lda, w, &wkopt, &lwork, rwork, &info );
        lwork = (int)wkopt.re;
        work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
        /* Solve eigenproblem */
        cheev( "Vectors", "Lower", &n, a, &lda, w, work, &lwork, rwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        /* Print eigenvalues */
        //print_rmatrix( "Eigenvalues", 1, n, w, 1 );
        /* Print eigenvectors */
        //print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
        /* Free workspace */
        free( (void*)work );
        exit( 0 );
} /* End of CHEEV Example */

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
