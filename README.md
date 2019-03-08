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
When I was trying to figure how much time does it take the MKL to diagonalise
a Hermitian matirx the segmentiation fault occured.
A diagonalisation of matrices of the size 1000x1000 complex numbers was immidiate (0.5s),
but for sizes larger thatn 1021x1021 the segfault occured.