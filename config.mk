# Compiler and compilation flags
cc=gcc-14
cxx=g++-14
cflags=-Wall -ansi -pedantic -fopenmp -O3 -std=c++11

# MPI compiler
mpicxx=mpicxx -Wno-long-long

# Flags for linking to PNG library
png_iflags=-DHAS_PNG
png_lflags=-lpng

# Flags for FFTW library
fftw_iflags=
fftw_lflags=-lfftw3

# LAPACK flags for dense linear algebra
lp_lflags=-llapack -lblas
