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

/* Main program */
int main(int argc, char *argv[]) {
        /* Locals */
        /* idk why the variable lda exists as it is always equal to n,
        but it was organised like it in the library so lets keep it this way.
        At leas at this iteration of the program */
        int n, lda, info, lwork;
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
               lda = n;
        }
        fcomplex wkopt;
        fcomplex* work;
		printf("The size of fcomplex is = %d\n",sizeof(fcomplex));
		printf("The size of float is = %d\n",sizeof(float));
        /* Local arrays */
        /* rwork dimension should be at least max(1,3*n-2) */
        float *w = (float*)malloc(sizeof(float)*n),
                *rwork = (float*)malloc(sizeof(float)*(3*n-2));
        if(w == NULL || rwork == NULL)
        {
                printf("Not enough memory to allocate w or rwork.\n");
                return 1;
        }

        fcomplex *a = (fcomplex*)malloc(sizeof(fcomplex)*lda*lda);
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
		printf("\n\nFull success in diagonalising %dx%d matirx\n",n,n);
		printf("Five largerst eigenvalues:\n");
        print_rmatrix( "Eigenvalues", 1, 5, w, 1 );
        //print_rmatrix( "Eigenvalues", 1, n, w, 1 );
        /* Print eigenvectors */
        //print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
        /* Free workspace */
        free( (void*)work );
        free( w );
        free( rwork );
        free( a );
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
