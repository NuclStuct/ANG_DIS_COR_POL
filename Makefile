# Define variables
FC = gfortran                 # Fortran compiler
FFLAGS = -O2 -g               # Compiler flags: -O2 for optimization, -g for debugging
OBJ = ang_dis_wm.o COM.o MATH.o QM_AM.o ANG_DIS.o  # Object files
EXEC = ang_dis  # Name of the final executable

# Default target: compile and link
all: $(EXEC)

# Link the object files into the final executable
$(EXEC): $(OBJ)
	$(FC) $(OBJ) -o $(EXEC)

# Compile the main program (ang_dis_wm.F90)
ang_dis_wm.o: ang_dis_wm.F90 COM.o MATH.o QM_AM.o ANG_DIS.o
	$(FC) $(FFLAGS) -c ang_dis_wm.F90

# Compile the COM module (COM.F90)
COM.o: COM.F90
	$(FC) $(FFLAGS) -c COM.F90

# Compile the MATH module (MATH.F90), which depends on COM.o
MATH.o: MATH.F90 COM.o
	$(FC) $(FFLAGS) -c MATH.F90

# Compile the QM_AM module (QM_AM.F90), which depends on MATH.o
QM_AM.o: QM_AM.F90 MATH.o
	$(FC) $(FFLAGS) -c QM_AM.F90

# Compile the ANG_DIS module (ANG_DIS.F90), which depends on QM_AM.o
ANG_DIS.o: ANG_DIS.F90 QM_AM.o
	$(FC) $(FFLAGS) -c ANG_DIS.F90


# Clean up the compiled files
clean:
	rm -f *.o *.mod $(EXEC)

# Phony targets (non-file targets)
.PHONY: all clean

