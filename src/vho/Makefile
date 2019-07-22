buildir = ..
include $(buildir)/make.config

objects = xvho input.o gauss.o H.o V.o fname.o fit.o

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

fit.o : fit.f90
	$(FC) $(FFLAGS) -c fit.f90

V.o : V.f90 fname.o input.o fit.o
	$(FC) $(FFLAGS) -c V.f90 fname.o input.o fit.o

H.o : H.f90 V.o
	$(FC) $(FFLAGS) -c H.f90 V.o

xvho : vho.f90 input.o gauss.o H.o V.o fname.o fit.o
	$(FC) $(FFLAGS) -o xvho vho.f90 input.o gauss.o H.o V.o fname.o fit.o
	
