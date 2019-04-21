#ifndef find_spectrum_h
#define find_spectrum_h

#include "fcomplex.h"
#include "versor.h" // in print_eigenvector_summary

/* CHEEV prototype */
extern void cheev( char* jobz, char* uplo, int* n, fcomplex* a, int* lda,
                float* w, fcomplex* work, int* lwork, float* rwork, int* info );
/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, fcomplex* a, int lda );
extern void print_rmatrix( char* desc, int m, int n, float* a, int lda );
void print_eigenvector_summary( fcomplex* a, int n, basis b, int m);
void print_lower_spectrum(fcomplex* a, int n, basis b, int m); 
void sort_print_lower_spectrum(fcomplex* a, int n, basis b, int m, fcomplex* work, int* work_int);
void sort_print_eigenvector_summary( fcomplex* a, fcomplex* work, int* work_int, int n, basis b, int m );
void sort_fcomplex(fcomplex* data, int* indices, int length);

#endif // find_spectrum_h