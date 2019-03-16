#ifndef find_spectrum_h
#define find_spectrum_h

/* Complex datatype */
struct _fcomplex { float re, im; };
typedef struct _fcomplex fcomplex;

/* CHEEV prototype */
extern void cheev( char* jobz, char* uplo, int* n, fcomplex* a, int* lda,
                float* w, fcomplex* work, int* lwork, float* rwork, int* info );
/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, fcomplex* a, int lda );
extern void print_rmatrix( char* desc, int m, int n, float* a, int lda );

#endif // find_spectrum_h