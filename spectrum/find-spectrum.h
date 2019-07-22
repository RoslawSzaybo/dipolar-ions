#ifndef find_spectrum_h
#define find_spectrum_h

#include "dcomplex.h"
#include "versor.h" // in print_eigenvector_summary

/* CHEEV prototype */
extern void zheev( char* jobz, char* uplo, int* n, dcomplex* a, int* lda,
                double* w, dcomplex* work, int* lwork, double* rwork, int* info );
/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, dcomplex* a, int lda );
extern void print_rmatrix( char* desc, int m, int n, double* a, int lda );
void print_eigenvalues( int first, int last, double* a );
void sort_print_lower_spectrum(dcomplex* a, int n, basis b, int m, 
                                dcomplex* work, int* work_int);
void sort_print_upper_spectrum(dcomplex* a, int n, basis b, int first, 
                                int last, dcomplex* work, int* work_int);
void sort_print_eigenvector_summary( dcomplex* a, dcomplex* work, 
                                    int* work_int, int n, basis b, int m );
void sort_dcomplex(dcomplex* data, int* indices, int length);

#endif // find_spectrum_h