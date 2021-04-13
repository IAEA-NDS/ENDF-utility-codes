# set default compiler

FC   =  gfortran
#FC  =  ifort
#FC  =  af90
#FC  =  lf95

# FFLAGS are based on compiler setting

ifeq ($(FC),gfortran)

  #-------------------------
  #----GNU FORTRAN compiler
  #-------------------------
  #----flags for production compilation
  FFLAGS = -std=legacy -O3
  #----flags for debuging
  #FFLAGS =  -g --bounds-check -std=legacy

else ifeq ($(FC),ifort)

  #-----------------------
  #----INTEL f95 compiler
  #-----------------------
  FFLAGS =  -O3 -x=host
  #----flags for debuging using ifort compiler
  #FFLAGS =  -C -ftrapuv -g -debug all
  #----flags for automatic parallelization
  #FFLAGS =  -parallel -openmp-report1 -par-threshold2 -openmp -vec-report1

else ifeq ($(FC),af90)

  #---------------------------
  #----Absoft FORTRAN compiler
  #---------------------------
  FFLAGS = -O3 -s
  #FFLAGS = -O0 -g

else ifeq ($(FC),lf90)

  #--------------------------------
  #----Lahey/Fujitsu f95 compiler
  #--------------------------------
  #----flag for Lahey/Fujitsu for production compilation
  FFLAGS = -O3
  #----flag for Lahey/Fujitsu for debugging
  #FFLAGS = -g --chk

endif

OBJF =	psyche.o

psyche:   $(OBJF)
	$(FC) -o $@ $(OBJF) $(STATIC)

clean:
	rm -f *.o
