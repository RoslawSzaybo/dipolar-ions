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



### Segmentation fauls:

For matrices larger than 1021x1021 there was a segfault in the original program. I have changed the variable handling i.e. in the original version, the size of a matrix was stored in `#define N xxx`, which is fortran style. I have changed it, so that the matrix size is read as a command line argument. Then I had to change the matrix size definitions from e.g. `fcomplex a[N*N]` to more secure `fcomplex* a = (fcomplex*)malloc(sizeof(fcomplex)*n*n)` and it turned out that is has resolved the segfault problem.

I was testing it on the cluster and for the matrix size ~48kx48k the segfault has returned and it is not solved yet.

### Statistics from plgrid:
- 8000x8000 Hermitian matrix: 500Mb; 7min;
- 16k x 16k Hermitian matrix: 2Gb; 51min;
- 32k x 32k Hermitian matrix: 8Gb; 7h 01min;
- 44k x 44k Hermitian matrix: 15Gb; 20h 25min;
