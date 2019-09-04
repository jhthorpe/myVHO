buildir = ..
include $(buildir)/make.config

objects = xcorires input.o quad.o cubi.o quar.o rota.o conv.o mu.o cori.o vstate.o calc.o ints.o term.o linal.o sort.o

all : $(objects)
	cp xcorires $(bindir)/xcorires

input.o : input.f90
	$(FC) $(FFLAGS) -c input.f90

conv.o : conv.f90
	$(FC) $(FFLAGS) -c conv.f90

vstate.o : vstate.f90
	$(FC) $(FFLAGS) -c vstate.f90

sort.o : sort.f90
	$(FC) $(FFLAGS) -c sort.f90

linal.o : linal.f90
	$(FC) $(FFLAGS) $(LINAL) -c linal.f90

ints.o : ints.f90 sort.o
	$(FC) $(FFLAGS) -c ints.f90 sort.o

term.o : term.f90 ints.o
	$(FC) $(FFLAGS) -c term.f90 ints.o

quad.o : quad.f90 input.o 
	$(FC) $(FFLAGS) -c quad.f90 input.o

cubi.o : cubi.f90 input.o
	$(FC) $(FFLAGS) -c cubi.f90 input.o

quar.o : quar.f90 input.o
	$(FC) $(FFLAGS) -c quar.f90 input.o

rota.o : rota.f90 conv.o input.o
	$(FC) $(FFLAGS) -c rota.f90 input.o conv.o

mu.o : mu.f90 input.o
	$(FC) $(FFLAGS) -c mu.f90 input.o

cori.o : cori.f90 input.o
	$(FC) $(FFLAGS) -c cori.f90 input.o

calc.o : calc.f90 conv.o ints.o term.o linal.o
	$(FC) $(FFLAGS) -c calc.f90 conv.o ints.o term.o linal.o

xcorires : corires.f90 input.o quad.o rota.o conv.o cubi.o quar.o mu.o cori.o vstate.o calc.o ints.o term.o linal.o sort.o
	$(FC) $(FFLAGS) $(LINAL) -o xcorires corires.f90 input.o quad.o rota.o conv.o cubi.o quar.o mu.o cori.o vstate.o calc.o ints.o term.o linal.o sort.o

clean : 
	rm *.o *.mod
