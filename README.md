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


Statistics from plgrid:
8000x8000 Hermitian matrix: 500Mb; 7min;
16k x 16k Hermitian matrix: 2Gb; 51min;

