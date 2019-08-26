#Makefile for Taylor with GCC-8.2.0 compilers

SHELL = /bin/bash

DEBUG = true#{true, false}  (NOTE: no whitespace at end)

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

SRC_DIR = ../src
BUILD_DIR = .
INSTALL_DIR = ../bin
EXECUTABLE = brownpak

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
		  m_utils_math       \
		  m_strings          \
		  m_logger           \
          m_ran_num          \
		  m_qsort            \
		  m_vector           \
		  m_table            \
		  m_trajectory       \
		  m_simbox           \
		  m_globals          \
		  m_control_io       \
		  m_config_io        \
		  m_stats_io         \
		  m_connectivity     \
		  m_verlet	         \
		  m_cell_list        \
		  m_ia_bond          \
		  m_ia_angle         \
		  m_ia_dihedral      \
		  m_ia_vdw           \
		  m_ia_tether        \
		  m_ia_external      \
		  m_interaction      \
		  m_relax            \
		  m_bd_solver        \
		  m_mpcd             \
		  m_setup            \
		  main  

OBJECTS = $(addsuffix .o,$(OBJ_NAMES))

move : $(EXECUTABLE)
	if [ -e "$(EXECUTABLE)" ]; then mv $(EXECUTABLE) $(INSTALL_DIR);fi 

$(EXECUTABLE) : $(OBJECTS)
	$(FC) $(INCS) $^ $(LIBS) -o $@

$(OBJECTS) : %.o : %.f90
	$(FC) $(INCS) $(FCFLAGS) $^ -o $@

docs :
	ford ../docs/brownpak.md 

.PHONY : clean clobber
clean:
	-rm -rf *.o *.mod
clobber:
	-rm -rf $(INSTALL_DIR)/$(EXECUTABLE) *.o *.mod
