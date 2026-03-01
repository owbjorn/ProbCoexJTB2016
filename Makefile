# A Makefile for the Runge-Kutta solver RKF45 and
#                ODE problem REPLICATOR.
# Note that the names of compilers and Makefile variables
# (macro names) are not standardized.

OBJS = rkf45.o

# Name of your compiler
FC = gfortran #ifort
CC = gcc

# Compiler options
FFLAGS = -O2#-pthread -O2 # for choleskyopen
#CFLAGS 

replicator: $(OBJS)
#	$(FC) -o $@ $(FFLAGS) choleskylinux.f90 $(OBJS) -llapack
#	$(FC) -o $@ $(FFLAGS) choleskyopen.f90 $(OBJS) -L./Netlib -lDPOTRF -L./OpenBLAS -lopenblas
	$(FC) -o $@ $(FFLAGS) replicator.f90 $(OBJS) #-L. -lrkf45 

rkf45.o: rkf45.f90
	$(FC) -c rkf45.f90

clean:
	rm -f replicator

# Some make-programs do not know about the f90 suffix,
# so we may have to add a rule for how to produce
# an object file (.o) from a Fortran90 file (.f90). 
.SUFFIXES: .f90
.f90.o:
	$(FC) -c $(FFLAGS) $<

