# Matrix diagonalisation.

### How to use
1. Install [Intel(R) Math Kernel Library](https://software.intel.com/mkl)

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

### Statistics of diagonalise.c from plgrid:
- 8000x8000 Hermitian matrix: 500Mb; 7min;
- 16k x 16k Hermitian matrix: 2Gb; 51min;
- 32k x 32k Hermitian matrix: 8Gb; 7h 01min;
- 44k x 44k Hermitian matrix: 15Gb; 20h 25min;
- 64k x 64k Hermitian matrix: 30.7Gb; 60h (2d 10h 30min);

### Statistics of the main program from plgrid:
- 42k x 42k Hermitian matrix: 13Gb; 11-22h ;

### Selected moleules
- MgH<sup>+</sup>
m = 25.3 u
q = 1 e
dipole = 3 D (Journal of Physics B: AMO, vol 42, no 15 (2009) M. Aymar et al.)
B = 180 000 MHz (9K) (New Journal of Physics 11 (2009) 055026)

- SrYb<sup>+</sup>
mass="261.0" u
charge="1.0" e
dipole="4.745" D
B="503.7" MHz
