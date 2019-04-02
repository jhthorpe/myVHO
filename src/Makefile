#local
export FC = gfortran
export FFLAGS = -O3 
#export FFLAGS = -g -fcsheck=bounds -Wall -Wextra -Wshadow

#local
#export LINAL = -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llapack -lblas

#HiPerGator
export MKL = /apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
export LINAL = -L$(MKL) -Wl,-R$(MKL) -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

#both
export bindir = ../bin

objects = input.o vho.o ints_HO.o plots.o linal.o fcs.o

all: $(objects)
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	cp vho.o $(bindir)/vho
	cp fcs.o $(bindir)/fcs

input.o : input.f90 
	$(FC) $(FFLAGS) -c input.f90

ints_HO.o : ints_HO.f90 
	$(FC) $(FFLAGS) -c ints_HO.f90

plots.o : plots.f90
	$(FC) $(FFLAGS) -c plots.f90

linal.o : linal.f90
	$(FC) $(FFLAGS) $(LINAL) -c linal.f90 

vho.o : vho.f90 input.o ints_HO.o plots.o linal.o
	$(FC) $(FFLAGS) -o vho.o vho.f90 input.o ints_HO.o plots.o linal.o $(LINAL)

fcs.o : fcs.f90
	$(FC) $(FFLAGS) -o fcs.o fcs.f90

clean:
	rm *.o *.mod


