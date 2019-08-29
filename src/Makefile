include make.config
export bindir = $(shell pwd)/bin

progdirs:= vho corires

.PHONY: all clean

all: 
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	@set -e; \
	for i in $(progdirs); do \
		if [ -d $$i ]; then\
			if [ -f $$i/Makefile ]; then \
				$(MAKE) -C $$i all; \
			fi; \
		fi; \
	done;

clean:
	@set -e;\
	for i in $(progdirs); do \
		if [ -d $$i ]; then\
			if [ -f $$i/Makefile ]; then \
				$(MAKE) -C $$i clean; \	
			fi; \
		fi; \
	done;

