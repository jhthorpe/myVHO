buildir = ..
include $(buildir)/make.config

objects = xvho


all : $(objects)
	cp xvho $(bindir)/xvho

clean : 
	rm *.o *.mod

xvho : vho.f90
	$(FC) $(FFLAGS) -o xvho vho.f90
	
