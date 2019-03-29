#local
export FC = gfortran
export FFLAGS = -O3
export LINAL = -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llapack -lblas

#both
export bindir = ../bin

objects = input.o vscf.o 

all: $(objects)
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	cp vscf.o $(bindir)/vscf

input.o : input.f90 
	$(FC) $(FFLAGS) -c input.f90

vscf.o : vscf.f90 input.o
	$(FC) $(FFLAGS) -o vscf.o vscf.f90 input.o

clean:
	rm *.o *.mod


