buildir = ..
include $(buildir)/make.config

objects = xvho input.o gauss.o H.o V.o fname.o fit.o val.o fcon.o basis.o

all : $(objects)
	cp xvho $(bindir)/xvho

clean : 
	rm *.o *.mod

input.o : input.f90
	$(FC) $(FFLAGS) -c input.f90

gauss.o : gauss.f90
	$(FC) $(FFLAGS) -c gauss.f90

fname.o : fname.f90
	$(FC) $(FFLAGS) -c fname.f90

val.o : val.f90
	$(FC) $(FFLAGS) -c val.f90

fcon.o : fcon.f90
	$(FC) $(FFLAGS) -c fcon.f90

fit.o : fit.f90
	$(FC) $(FFLAGS) -c fit.f90

basis.o : basis.f90 input.o
	$(FC) $(FFLAGS) -c basis.f90 input.o

V.o : V.f90 fname.o input.o fit.o val.o
	$(FC) $(FFLAGS) -c V.f90 fname.o input.o fit.o val.o

H.o : H.f90 V.o fcon.o basis.o
	$(FC) $(FFLAGS) -c H.f90 V.o fcon.o basis.o

xvho : vho.f90 input.o gauss.o H.o V.o fname.o fit.o val.o fcon.o basis.o
	$(FC) $(FFLAGS) -o xvho vho.f90 input.o gauss.o H.o V.o fname.o fit.o val.o fcon.o basis.o
	
