# =====
# Setup
# =====

# Install location (note that bin will be appended)
PREFIX = /usr/local/

# The compiler and flags (uncomment extra for debugging)
FC = ifort
FCFLAGS += -O2 -Iinclude #-g -traceback -warn all -check all

# Choose the flag that defines the module directory for this compiler
ifeq (${FC}, ifort)
  FCFLAGS += -module include
else
  ifeq (${FC}, gfortran)
    FCFLAGS += -Jinclude
  else
    ifeq (${FC}, pgf90)
      FCFLAGS += -module include
    else
      ifeq (${FC}, g95)
        FCFLAGS += -fmod=include
      endif
    endif
  endif
endif

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

SRC = src
BLD = build

MODS = ${BLD}/constants.o \
       ${BLD}/parallel.o \
       ${BLD}/global.o

SUBOBJ = ${MODS} \
         ${BLD}/CollectCoordinates.o \
         ${BLD}/GetOptions.o \
         ${BLD}/quit.o \
         ${BLD}/Random.o

OBJ = ${BLD}/npvol.o ${SUBOBJ}

# "make" builds all
all: NPVol

# General rule for building executable from objects. 
# $@ is the name of the target (in this case the executable)
# $^ is all dependencies
NPVol: ${OBJ}
	${FCCOMP} ${FCFLAGS} -o $@ $^

# $< is used in order to list only the first dependency (the source file)
# and not the additional prerequisites such as module or include files
${BLD}/npvol.o: ${SRC}/npvol.f90 ${SUBOBJ}
	${FCCOMP} ${FCFLAGS} -c $< -o $@

${BLD}/quit.o: ${SRC}/quit.f90 ${MODS}
	${FCCOMP} ${FCFLAGS} -c $< -o $@

${BLD}/CollectCoordinates.o: ${SRC}/CollectCoordinates.f90 ${MODS}
	${FCCOMP} ${FCFLAGS} -c $< -o $@

${BLD}/GetOptions.o: ${SRC}/GetOptions.f90
	${FCCOMP} ${FCFLAGS} -c $< -o $@

${BLD}/Random.o: ${SRC}/Random.f90 ${BLD}/constants.o
	${FCCOMP} ${FCFLAGS} -c $< -o $@

${BLD}/parallel.o: ${SRC}/parallel.F90 ${BLD}/constants.o
	${FCCOMP} ${FCFLAGS} -c $< -o $@

${BLD}/global.o: ${SRC}/global.f90 ${BLD}/constants.o
	${FCCOMP} ${FCFLAGS} -c $< -o $@

${BLD}/constants.o: ${SRC}/constants.f90 ${BLD}/build-dir-exists
	${FCCOMP} ${FCFLAGS} -c $< -o $@

${BLD}/build-dir-exists:
	mkdir ${BLD}
	touch ${BLD}/build-dir-exists

# =============
# Special rules
# =============
#
#  Make sure we don't remove this by accident if interrupted at the wrong time.
.PRECIOUS: Makefile

# Utility targets
.PHONY: install clean cleanall

# Installs executable to your home binary directory
install:
	install -D ${PROG} ${PREFIX}/${PROG}

# Removes objects
clean:
	@rm -vf build/*.o include/*.mod

# Removes executable
cleanall: clean
	@rm -vf ${PROG}
	@rm -vrf build/
