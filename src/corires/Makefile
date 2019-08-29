buildir = ..
include $(buildir)/make.config

objects = xcorires input.o quad.o cubi.o rota.o conv.o mu.o cori.o vstate.o calc.o

all : $(objects)
	cp xcorires $(bindir)/xcorires

input.o : input.f90
	$(FC) $(FFLAGS) -c input.f90

conv.o : conv.f90
	$(FC) $(FFLAGS) -c conv.f90

vstate.o : vstate.f90
	$(FC) $(FFLAGS) -c vstate.f90

quad.o : quad.f90 input.o 
	$(FC) $(FFLAGS) -c quad.f90 input.o

cubi.o : cubi.f90 input.o
	$(FC) $(FFLAGS) -c cubi.f90 input.o

calc.o : calc.f90 conv.o
	$(FC) $(FFLAGS) -c calc.f90 conv.o

rota.o : rota.f90 conv.o input.o
	$(FC) $(FFLAGS) -c rota.f90 input.o conv.o

mu.o : mu.f90 input.o
	$(FC) $(FFLAGS) -c mu.f90 input.o

cori.o : cori.f90 input.o
	$(FC) $(FFLAGS) -c cori.f90 input.o

xcorires : corires.f90 input.o quad.o rota.o conv.o cubi.o mu.o cori.o vstate.o calc.o
	$(FC) $(FFLAGS) -o xcorires corires.f90 input.o quad.o rota.o conv.o cubi.o mu.o cori.o vstate.o calc.o

clean : 
	rm *.o *.mod
