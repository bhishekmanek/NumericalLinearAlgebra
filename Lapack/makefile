#                Modified by Bhishek Manek, AMS 213A , Winter 2017
#                Adapted from Prof. Dongwook's lecture notes for AMS 209
########################################################################################

FC      = gfortran

FFLAGS_DEBUG = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace

LFLAGS = -lblas -llapack

#FFLAGS_OPT = -ggdb -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -Wuninitialized

EXE_FILE = run_lp1.exe

OBJS  = linal.o \
        Driver_LAPACK1.o

.PHONY: clean
########################################################################################
#COMPILING AND LINKING USING GENERIC SUFFIX RULE FOR F90

$(EXE_FILE) : $(OBJS)
	@$(FC) $(FFLAGS_OPT) $(OBJS) $(LFLAGS) -o $(EXE_FILE)
	@echo "code is now linking..."

#LET'S APPLY GENERIC SUFFIX RULE HERE FOR FORTRAN 90

# method 1 using generic suffix rule
#.SUFFIXES: 
#.SUFFIXES: .F90 .o
#.F90.o:
#       $(FC) $(FFLAGS_OPT) -c $<

# method 2 using inference rule
%.o: %.F90
	$(FC) $(FFLAGS_OPT) -c $<

#######################################################################################
#SOME USEFUL COMMANDS
clean:
	@rm -f *.o *.mod *~ run_lp1.exe

# debug: clean
debug: FFLAGS_OPT = $(FFLAGS_DEBUG)
debug: $(EXE_FILE)

#######################################################################################
#LET'S DEFINE SOME MODULE DEPENDENCIES!
Driver_LAPACK1.o: linal.o

#######################################################################################
