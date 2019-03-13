Matrix diagonalisation.

How to use
1. Install Intel(R) MKL

2. Setup local variables, by default it's
``` bash
	source /opt/intel/mkl/bin/mklvars.sh # libraries and stuff
	export PATH="/opt/intel/bin:$PATH" # for icc
```

3. Compilation and linkning
``` bash
	icc -mkl dsyev.c
```



The first problem:
I was trying to figure out how much time does it take the MKL to diagonalise
a Hermitian matirx. So I was running the an edited verion of the cheev on larger and larger matrices.
At the matrix size 1021 the segmentiation fault occured.

To make the code run, I didn't change anything. All I did was just to move the code to the plgrid cluster
and there all was working good.

The problem got resolved automatically as soon as I have changed the variable handling i.e. in the original
version the size of a matrix is stored in `#define N xxx`, which is fortran-alike. I have changed it, so that
the matrix size is read as a command line argument. Then I had to change the matrix size definitions from e.g. 
`fcomplex a[N*N]` to more secure `fcomplex* a = (fcomplex*)malloc(sizeof(fcomplex)*n*n)` 
and it turned out that is has resolved the segfault problem.
