 ### Segmentation fauls:
 24 
 25 For matrices larger than 1021x1021 there was a segfault in the original program. I have changed the variable handling i.e. in the original version, the size of a matrix was stored in `#define N xxx`, which i    s fortran style. I have changed it, so that the matrix size is read as a command line argument. Then I had to change the matrix size definitions from e.g. `fcomplex a[N*N]` to more secure `fcomplex* a = (fco    mplex*)malloc(sizeof(fcomplex)*n*n)` and it turned out that is has resolved the segfault problem.
 26 
 27 I was testing it on the cluster and for the matrix size ~48kx48k the segfault has returned and it is not solved yet.

