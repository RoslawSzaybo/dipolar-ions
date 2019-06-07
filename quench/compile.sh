#!/bin/bash

icc -mkl -c ../spectrum/input.c -o input.o
icc -mkl -c ../spectrum/versor.c -o versor.o
icc -mkl -c ../spectrum/state.c -o state.o
icc -mkl -c ../spectrum/fcomplex.c -o fcomplex.o
icc -mkl -c ../spectrum/hamiltonian.c -o hamiltonian.o
icc -mkl -c main.c -o main.o
icc -mkl *.o
rm *.o

./a.out
