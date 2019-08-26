#Makefile for Taylor with GCC-8.2.0 compilers

SHELL = /bin/bash

DEBUG = false#{true, false}  (NOTE: no whitespace at end)

FC = $(gfortran8)

#*************************************************

ifeq ($(DEBUG), true)
FCFLAGS = -c -O0 -Wall -fcheck=all -fbacktrace -std=f2008ts \
		  -fdefault-double-8 \
		  -fdefault-real-8 \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
else
FCFLAGS = -c -O3 -std=f2008ts \
		  -fdefault-double-8 \
		  -fdefault-real-8 \
		  -funroll-loops \
		  -flto          \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
endif

SRC_DIR = ../../src
BUILD_DIR = .
INSTALL_DIR = .
EXECUTABLE = post

vpath %.f90 $(SRC_DIR)
vpath $(EXECUTABLE) $(INSTALL_DIR)

MKL_INTRFC = $(HOME)/soft/mkl/gcc-8

INCS = \
	-I$(SRC_DIR) \
	-I$(MKL_INTRFC)/include \
	-I$(MKLROOT)/include

LIBS = \
	$(MKL_INTRFC)/lib/libmkl_blas95_lp64.a \
	$(MKL_INTRFC)/lib/libmkl_lapack95_lp64.a \
	-L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib \
	-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread \
	-lz -ldl -lm

OBJ_NAMES = \
		  m_precision        \
		  m_constants_math   \
		  m_strings          \
		  m_logger           \
		  m_trajectory       \
		  m_simbox           \
		  m_globals          \
		  m_config_io        \
		  post  

OBJECTS = $(addsuffix .o,$(OBJ_NAMES))

$(EXECUTABLE) : $(OBJECTS)
	$(FC) $(INCS) $^ $(LIBS) -o $@

$(OBJECTS) : %.o : %.f90
	$(FC) $(INCS) $(FCFLAGS) $^ -o $@

.PHONY : clean clobber
clean:
	-rm -rf *.o *.mod
clobber:
	-rm -rf $(INSTALL_DIR)/$(EXECUTABLE) *.o *.mod
