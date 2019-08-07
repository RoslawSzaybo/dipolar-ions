#!/bin/bash

icc -mkl -c ../input.c -o input.o
icc -mkl -c ../versor.c -o versor.o
icc -mkl -c ../state.c -o state.o
icc -mkl -c ../dcomplex.c -o dcomplex.o
icc -mkl -c ../hamiltonian.c -o hamiltonian.o
icc -mkl -c main.c -o main.o
icc -mkl *.o

rm *.o

./a.out
