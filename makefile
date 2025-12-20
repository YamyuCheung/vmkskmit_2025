# ==========================================
# Makefile for Wave Equation Solver (Polus)
# ==========================================

# [System Config] Use bash for module command
SHELL := /bin/bash

# [Environment] Load MPI module
PRE_CMD = module load SpectrumMPI/10.1.0

# [Compilers]
CXX      = xlc++_r
MPICXX   = mpixlC

# [Flags]
CXXFLAGS = -std=c++11 -O3 -qarch=pwr8 -qtune=pwr8 -qhot -qipa
OMPFLAG  = -qsmp=omp

# --- Output Names ---
TARGET_OMP    = wave_eq_variant7
TARGET_MPI    = wave_mpi
TARGET_HYBRID = wave_hybrid

# --- Source Files ---
SRC_OMP    = wave_eq_variant7_omp.cpp
SRC_MPI    = wave_eq_variant7_mpi.cpp
SRC_HYBRID = wave_eq_variant7_hybrid.cpp

# --- Target Folders (Based on your image) ---
DIR_OMP    = script_polus_omp_version
DIR_MPI    = script_polus_mpi_version
DIR_HYBRID = script_polus_hybrid_version

# ==========================================
# Rules
# ==========================================

.PHONY: all clean omp mpi hybrid

# Build all versions
all: omp mpi hybrid

# 1. OpenMP Implementation -> move to script_polus_omp_version/
omp: $(SRC_OMP)
	$(PRE_CMD) && $(CXX) $(CXXFLAGS) $(OMPFLAG) -o $(TARGET_OMP) $(SRC_OMP)
	@mv -f $(TARGET_OMP) $(DIR_OMP)/
	@echo "Build successful: $(DIR_OMP)/$(TARGET_OMP)"

# 2. MPI Implementation -> move to script_polus_mpi_version/
mpi: $(SRC_MPI)
	$(PRE_CMD) && $(MPICXX) $(CXXFLAGS) -o $(TARGET_MPI) $(SRC_MPI)
	@mv -f $(TARGET_MPI) $(DIR_MPI)/
	@echo "Build successful: $(DIR_MPI)/$(TARGET_MPI)"

# 3. Hybrid Implementation -> move to script_polus_hybrid_version/
hybrid: $(SRC_HYBRID)
	$(PRE_CMD) && $(MPICXX) $(CXXFLAGS) $(OMPFLAG) -o $(TARGET_HYBRID) $(SRC_HYBRID)
	@mv -f $(TARGET_HYBRID) $(DIR_HYBRID)/
	@echo "Build successful: $(DIR_HYBRID)/$(TARGET_HYBRID)"

# Clean up (Delete files inside the folders)
clean:
	rm -f $(DIR_OMP)/$(TARGET_OMP)
	rm -f $(DIR_MPI)/$(TARGET_MPI)
	rm -f $(DIR_HYBRID)/$(TARGET_HYBRID)
	@echo "Cleaned executables from subdirectories."