buildir = ..
include $(buildir)/make.config

objects = xvho input.o

all : $(objects)
	cp xvho $(bindir)/xvho

clean : 
	rm *.o *.mod

input.o : input.f90
	$(FC) $(FFLAGS) -c input.f90

xvho : vho.f90 input.o
	$(FC) $(FFLAGS) -o xvho vho.f90 input.o
	
