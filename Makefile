# =====
# Setup
# =====

# The compiler and flags (uncomment extra for debugging)
FC = ifort
FCFLAGS += -O2 #-g -traceback -warn all -check all

# Install location (note that bin will be appended)
PREFIX = /usr/local/

# Parallel mode we want to use
# Define mode to openmp or mpi for OpenMP or MPI, respectively
PARMODE = serial
# Set this to define the MPI flavor
# So far only OpenMPI and MPICH are inmplemented
FLAVOR = openmpi

# Define the compiler based on MPI or not
ifeq (${PARMODE}, mpi)
  FCFLAGS += -D_MPI
  ifeq (${FLAVOR}, openmpi)
    FCCOMP := OMPI_FC=${FC} mpif90
  else
    ifeq (${FLAVOR}, mpich)
      FCCOMP := mpif90 -f90=${FC}
    endif
  endif
else
  # If OpenMP, select the correct flags
  ifeq (${PARMODE}, openmp)
    ifeq (${FC}, ifort)
      FCFLAGS += -openmp
    else
      ifeq (${FC}, gfortran)
        FCFLAGS += -fopenmp -lgomp
      else
        ifeq (${FC}, pgf90)
          FCFLAGS += -mp
        endif
          # g95 doesn't support OpenMP
      endif
    endif
  endif
  # Reguardless of serial or OpenMP, use compiler as-is
  FCCOMP := ${FC}
endif

# ==============
# Building rules
# ==============

# "make" builds all
all: NPVol

# General rule for building executable from objects. 
# $@ is the name of the target (in this case the executable)
# $^ is all dependencies
NPVol: npvol.o constants.o CollectCoordinates.o GetOptions.o global.o parallel.o quit.o Random.o
	${FCCOMP} ${FCFLAGS} -o $@ $^

# $< is used in order to list only the first dependency (the source file)
# and not the additional prerequisites such as module or include files
npvol.o: npvol.f90 constants.o CollectCoordinates.o GetOptions.o global.o parallel.o quit.o Random.o
	${FCCOMP} ${FCFLAGS} -c $< -o $@

quit.o: quit.f90 constants.o parallel.o global.o
	${FCCOMP} ${FCFLAGS} -c $< -o $@

CollectCoordinates.o: CollectCoordinates.f90 constants.o parallel.o global.o
	${FCCOMP} ${FCFLAGS} -c $< -o $@

GetOptions.o: GetOptions.f90
	${FCCOMP} ${FCFLAGS} -c $< -o $@

Random.o: Random.f90 constants.o
	${FCCOMP} ${FCFLAGS} -c $< -o $@

parallel.o: parallel.F90 constants.o
	${FCCOMP} ${FCFLAGS} -c $< -o $@

global.o: global.f90 constants.o
	${FCCOMP} ${FCFLAGS} -c $< -o $@

constants.o: constants.f90
	${FCCOMP} ${FCFLAGS} -c $< -o $@

# =============
# Special rules
# =============
#
#  Make sure we don't remove this by accident if interrupted at the wrong time.
.PRECIOUS: Makefile

# Utility targets
.PHONY: install clean cleanall

debug:
	@echo ${FLAVOR} ${PARMODE} ${FCCOMP}

# Installs executable to your home binary directory
install:
	install -D ${PROG} ${PREFIX}/${PROG}

# Removes objects
clean:
	@rm -vf *.o *.mod

# Removes executable
cleanall: clean
	@rm -vf ${PROG}
