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
	icc -mkl *.c
```

4. Program execution (this one applies only to the file in tests)
``` bash
	./a.out 7000 # the matirx size
```

### Statistics from plgrid:
- 8000x8000 Hermitian matrix: 500Mb; 7min;
- 16k x 16k Hermitian matrix: 2Gb; 51min;
- 32k x 32k Hermitian matrix: 8Gb; 7h 01min;
- 44k x 44k Hermitian matrix: 15Gb; 20h 25min;
- 64k x 64k Hermitian matrix: 30.7Gb; 60h (2d 10h 30min);

### todo
- implement `void state_clean()`
