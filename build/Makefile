#Makefile for Taylor with GCC-9.2.0 compilers

SHELL = /bin/bash

DEBUG = false#{true, false}  (NOTE: no whitespace at end)

FC = $(gfortran9)

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
		  -flto \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
endif

SRC_DIR = ../src
BUILD_DIR = .
INSTALL_DIR = ../bin
LIB_DIR = ../lib
EXECUTABLE = brownpak

vpath %.f90 $(SRC_DIR)
vpath %.fpp $(SRC_DIR)
vpath $(EXECUTABLE) $(INSTALL_DIR)

MKL_INTRFC = $(HOME)/soft/mkl/gcc-9

INCS = \
	-I$(SRC_DIR) \
	-I$(MKL_INTRFC)/include
	#-I$(MKLROOT)/include

LIBS = \
	$(LIB_DIR)/libnitsol.a \
	$(MKL_INTRFC)/lib/libmkl_blas95_lp64.a \
	$(MKL_INTRFC)/lib/libmkl_lapack95_lp64.a \
	-L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib \
	-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread \
	-lz -ldl -lm

FYPPFLAGS = #-D DEBUG=1

OBJ_NAMES = \
		  constants_m    \
		  utils_math_m   \
		  strings_m      \
		  timestamp_m    \
		  logger_m       \
          random_m       \
		  qsort_m        \
		  vector_m       \
		  table_m        \
		  aabb_m         \
		  aabbtree_m     \
		  aabbtree_sm    \
		  trajectory_m   \
		  simbox_m       \
		  atmcfg_m       \
		  control_m      \
		  config_io_m    \
		  stats_m        \
		  connectivity_m \
		  cell_list_m    \
		  pairtab_m      \
		  ia_bond_m      \
		  ia_angle_m     \
		  ia_dihedral_m  \
		  ia_vdw_m       \
		  ia_tether_m    \
		  ia_external_m  \
		  interaction_m  \
		  brown_m        \
		  bd_solver_m    \
		  setup_m        \
		  main  

OBJECTS = $(addsuffix .o,$(OBJ_NAMES))

move : $(EXECUTABLE)
	if [ -e "$(EXECUTABLE)" ]; then mv $(EXECUTABLE) $(INSTALL_DIR);fi 

$(EXECUTABLE) : $(OBJECTS)
	$(FC) $(INCS) $^ $(LIBS) -o $@

$(OBJECTS) : %.o : %.f90
	$(FC) $(INCS) $(FCFLAGS) $^ -o $@

%.f90 : %.fpp
	fypp $(FYPPFLAGS) $< $@

%.o : %.f90
	$(FC) $(INCS) $(FCFLAGS) $^ -o $@ 

docs :
	ford ../docs/brownpak.md 

.PHONY : clean
clean:
	-rm -rf *.o *.mod *.smod $(BUILD_DIR)/*.f90

.PHONY : clobber
clobber:
	-rm -rf $(INSTALL_DIR)/$(EXECUTABLE) *.o *.mod *.smod $(BUILD_DIR)/*.f90
