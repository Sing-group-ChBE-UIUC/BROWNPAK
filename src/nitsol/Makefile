SHELL = /bin/bash
DEBUG = false#

FC = $(gfortran9)

ifeq ($(DEBUG), true)
FCFLAGS = -c -O0 -Wall -fcheck=all -fbacktrace  \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
else
FCFLAGS = -c -O3 -funroll-loops -flto \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
endif

all: libnitsol.a

NITSOL = nitsol.o nitdrv.o nitjv.o nitbt.o nitgm.o nitstb.o nittfq.o nitbd.o

# Dependencies.

nitsol.o:	nitsol.f nitprint.h nitinfo.h nitparam.h nitdflts.h
nitdrv.o:	nitdrv.f nitprint.h nitinfo.h nitparam.h
nitjv.o:	nitjv.f nitprint.h
nitbt.o:	nitbt.f nitprint.h nitparam.h
nitgm.o:	nitgm.f nitprint.h
nitstb.o:	nitstb.f nitprint.h
nittfq.o: 	nittfq.f nitprint.h
nitbd.o:	nitbd.f nitprint.h nitparam.h nitdflts.h

# Default compilation rules.

.f.o:
	$(FC) -c $(FFLAGS) $*.f

# Rules to build library.

libnitsol.a:  $(NITSOL)
	ar rs $@ $?

clean:
	- /bin/rm $(NITSOL) libnitsol.a

remake:
	make clean
	make
