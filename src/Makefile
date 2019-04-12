#local
export FC = gfortran
#export FFLAGS = -O3 
export FFLAGS = -g -fcheck=bounds -Wall -Wextra -Wshadow

#local
#export LINAL = -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llapack -lblas

#HiPerGator
export MKL = /apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
export LINAL = -L$(MKL) -Wl,-R$(MKL) -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

#both
export bindir = ../bin

objects = input.o vho.o ints_HO.o proc.o linal.o fcs.o tdp.o conv.o scat.o nints.o fit.o val.o

all: $(objects)
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	cp vho.o $(bindir)/vho
	cp fcs.o $(bindir)/fcs
	cp tdp.o $(bindir)/tdp	
	cp conv.o $(bindir)/conv
	cp scat.o $(bindir)/scat

input.o : input.f90 
	$(FC) $(FFLAGS) -c input.f90

ints_HO.o : ints_HO.f90 nints.o fit.o val.o
	$(FC) $(FFLAGS) -c ints_HO.f90 nints.o fit.o val.o

proc.o : proc.f90
	$(FC) $(FFLAGS) -c proc.f90

linal.o : linal.f90
	$(FC) $(FFLAGS) $(LINAL) -c linal.f90 

nints.o : nints.f90
	$(FC) $(FFLAGS) -c nints.f90

fit.o : fit.f90 linal.o val.o nints.o
	$(FC) $(FFLAGS) $(LINAL) -c fit.f90 linal.o val.o nints.o

val.o : val.f90
	$(FC) $(FFLAGS) -c val.f90

vho.o : vho.f90 input.o ints_HO.o proc.o linal.o nints.o val.o fit.o
	$(FC) $(FFLAGS) -o vho.o vho.f90 input.o ints_HO.o nints.o proc.o fit.o val.o linal.o $(LINAL)

fcs.o : fcs.f90
	$(FC) $(FFLAGS) -o fcs.o fcs.f90

tdp.o : tdp.f90
	$(FC) $(FFLAGS) -o tdp.o tdp.f90

conv.o : conv.f90
	$(FC) $(FFLAGS) -o conv.o conv.f90

scat.o : scat.f90
	$(FC) $(FFLAGS) -o scat.o scat.f90

clean:
	rm *.o *.mod


