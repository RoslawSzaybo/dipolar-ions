#!/bin/bash

date=$(date +"%m.%d")

gcc -c ../spectrum/input.c -o input.o
gcc -c ../spectrum/versor.c -o versor.o
gcc -c ../spectrum/state.c -o state.o
gcc -c ../spectrum/dcomplex.c -o dcomplex.o
gcc -c ../spectrum/hamiltonian.c -o hamiltonian.o
gcc -c main.c -o main.o
gcc *.o -lm -o quench.$date
rm *.o
