buildir = ..
include $(buildir)/make.config

objects = xvho input.o gauss.o H_HO.o V.o fname.o fit.o val.o fcon.o basis.o ints_HO.o memory.o evec.o sort.o key.o rota.o cori.o quad.o cubi.o quar.o mome.o

all : $(objects)
	cp xvho $(bindir)/xvho

clean : 
	rm *.o *.mod

input.o : input.f90
	$(FC) $(FFLAGS) -c input.f90

fname.o : fname.f90
	$(FC) $(FFLAGS) -c fname.f90

val.o : val.f90
	$(FC) $(FFLAGS) -c val.f90

fit.o : fit.f90
	$(FC) $(FFLAGS) -c fit.f90

key.o : key.f90
	$(FC) $(FFLAGS) -c key.f90

sort.o : sort.f90
	$(FC) $(FFLAGS) -c sort.f90

quad.o : quad.f90
	$(FC) $(FFLAGS) -c quad.f90

cubi.o : cubi.f90
	$(FC) $(FFLAGS) -c cubi.f90

quar.o : quar.f90
	$(FC) $(FFLAGS) -c quar.f90

mome.o : mome.f90
	$(FC) $(FFLAGS) -c mome.f90

linal.o : linal.f90
	$(FC) $(FFLAGS) $(LINAL) -c linal.f90

rota.o : rota.f90 input.o
	$(FC) $(FFLAGS) -c rota.f90 input.o

fcon.o : fcon.f90 sort.o
	$(FC) $(FFLAGS) -c fcon.f90 sort.o

evec.o	: evec.f90 key.o
	$(FC) $(FFALGS) -c evec.f90 key.o

basis.o : basis.f90 input.o
	$(FC) $(FFLAGS) -c basis.f90 input.o

cori.o : cori.f90 input.o sort.o
	$(FC) $(FFLAGS) -c cori.f90 input.o sort.o

memory.o : memory.f90 val.o linal.o
	$(FC) $(FFLAGS) -c memory.f90 val.o linal.o

gauss.o : gauss.f90 input.o fname.o key.o
	$(FC) $(FFLAGS) -c gauss.f90 input.o fname.o key.o

V.o : V.f90 fname.o input.o fit.o val.o key.o
	$(FC) $(FFLAGS) -c V.f90 fname.o input.o fit.o val.o key.o

ints_HO.o : ints_HO.f90 val.o gauss.o key.o quad.o cubi.o quar.o mome.o
	$(FC) $(FFLAGS) -c ints_HO.f90 val.o gauss.o key.o quad.o cubi.o quar.o mome.o

H_HO.o : H_HO.f90 V.o fcon.o basis.o ints_HO.o memory.o linal.o evec.o key.o rota.o cori.o
	$(FC) $(FFLAGS) -c H_HO.f90 V.o fcon.o basis.o ints_HO.o memory.o linal.o evec.o key.o rota.o cori.o

xvho : vho.f90 input.o gauss.o H_HO.o V.o fname.o fit.o val.o fcon.o basis.o ints_HO.o memory.o linal.o evec.o sort.o key.o rota.o cori.o quad.o cubi.o quar.o mome.o
	$(FC) $(FFLAGS) $(LINAL) -o xvho vho.f90 input.o gauss.o H_HO.o V.o fname.o fit.o val.o fcon.o basis.o ints_HO.o memory.o linal.o evec.o sort.o key.o rota.o cori.o quad.o cubi.o quar.o mome.o
	
