#local
export FC = gfortran
export FFLAGS = -O3
export LINAL = -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llapack -lblas

#both
export bindir = ../bin

objects = input.o vscf.o ints_HO.o plots.o

all: $(objects)
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	cp vscf.o $(bindir)/vscf

input.o : input.f90 
	$(FC) $(FFLAGS) -c input.f90

ints_HO.o : ints_HO.f90 
	$(FC) $(FFLAGS) -c ints_HO.f90

plots.o : plots.f90
	$(FC) $(FFLAGS) -c plots.f90

vscf.o : vscf.f90 input.o ints_HO.o plots.o
	$(FC) $(FFLAGS) -o vscf.o vscf.f90 input.o ints_HO.o plots.o

clean:
	rm *.o *.mod


