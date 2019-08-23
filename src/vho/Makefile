buildir = ..
include $(buildir)/make.config

objects = xvho input.o gauss.o H_HO.o V.o fname.o fit.o valu.o fcon.o basis.o ints_HO.o memory.o evec.o sort.o key.o rota.o cori.o quad.o cubi.o quar.o mome.o pseu.o didq.o test

all : $(objects)
	cp xvho $(bindir)/xvho
	cp test $(bindir)/testi

clean : 
	rm *.o *.mod

input.o : input.f90
	$(FC) $(FFLAGS) -c input.f90

fname.o : fname.f90
	$(FC) $(FFLAGS) -c fname.f90

valu.o : valu.f90
	$(FC) $(FFLAGS) -c valu.f90

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

memory.o : memory.f90 valu.o linal.o
	$(FC) $(FFLAGS) -c memory.f90 valu.o linal.o

didq.o : didq.f90 sort.o input.o
	$(FC) $(FFLAGS) -c didq.f90 sort.o input.o

cori.o : cori.f90 input.o sort.o ints_HO.o
	$(FC) $(FFLAGS) -c cori.f90 input.o sort.o ints_HO.o

gauss.o : gauss.f90 input.o fname.o key.o
	$(FC) $(FFLAGS) -c gauss.f90 input.o fname.o key.o

V.o : V.f90 fname.o input.o fit.o valu.o key.o
	$(FC) $(FFLAGS) -c V.f90 fname.o input.o fit.o valu.o key.o

ints_HO.o : ints_HO.f90 valu.o gauss.o key.o quad.o cubi.o quar.o mome.o 
	$(FC) $(FFLAGS) -c ints_HO.f90 valu.o gauss.o key.o quad.o cubi.o quar.o mome.o 

pseu.o : pseu.f90 ints_HO.o sort.o
	$(FC) $(FFLAGS) -c pseu.f90 ints_HO.o sort.o

H_HO.o : H_HO.f90 V.o fcon.o basis.o ints_HO.o memory.o linal.o evec.o key.o rota.o cori.o sort.o pseu.o didq.o
	$(FC) $(FFLAGS) -c H_HO.f90 V.o fcon.o basis.o ints_HO.o memory.o linal.o evec.o key.o rota.o cori.o sort.o pseu.o didq.o

xvho : vho.f90 input.o gauss.o H_HO.o V.o fname.o fit.o valu.o fcon.o basis.o ints_HO.o memory.o linal.o evec.o sort.o key.o rota.o cori.o quad.o cubi.o quar.o mome.o pseu.o didq.o
	$(FC) $(FFLAGS) $(LINAL) -o xvho vho.f90 input.o gauss.o H_HO.o V.o fname.o fit.o valu.o fcon.o basis.o ints_HO.o memory.o linal.o evec.o sort.o key.o rota.o cori.o quad.o cubi.o quar.o mome.o pseu.o didq.o
	
test : test.f90 cori.o ints_HO.o sort.o input.o valu.o gauss.o key.o quad.o cubi.o quar.o mome.o fname.o
	$(FC) $(FFLAGS) $(LINAL) -o test test.f90 cori.o ints_HO.o sort.o input.o valu.o gauss.o key.o quad.o cubi.o quar.o mome.o fname.o
