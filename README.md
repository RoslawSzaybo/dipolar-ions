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
	icc -mkl dsyev.cpp 
```
