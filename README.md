# Matrix diagonalisation.

### How to use
1. Install Intel(R) MKL

2. Setup local variables, by default it's
``` bash
	source /opt/intel/mkl/bin/mklvars.sh # libraries and stuff
	export PATH="/opt/intel/bin:$PATH" # for icc
```

3. Compilation and linkning
``` bash
	icc -mkl diagonalise.c
```

4. Program execution
``` bash
	./a.out 7000 # the matirx size
```



### The first problem:
I was trying to figure out how much time does it take the MKL to diagonalise
a Hermitian matirx. So I was running an edited version of the `cheev` example on larger and larger matrices.
At the matrix size 1021 the segmentiation fault occured. However it was not a problem for cluster, that was
only present on my pc.

The problem got resolved automatically as soon as I have changed the variable handling i.e. in the original
version the size of a matrix is stored in `#define N xxx`, which is fortran-alike. I have changed it, so that
the matrix size is read as a command line argument. Then I had to change the matrix size definitions from e.g. 
`fcomplex a[N*N]` to more secure `fcomplex* a = (fcomplex*)malloc(sizeof(fcomplex)*n*n)` 
and it turned out that is has resolved the segfault problem.

### Statistics from plgrid:
8000x8000 Hermitian matrix: 500Mb; 7min;
16k x 16k Hermitian matrix: 2Gb; 51min;
