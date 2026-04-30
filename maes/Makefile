# =======================================================
# CHOIX DU COMPILATEUR
# =======================================================

FC = gfortran

# =======================================================
# ORGANISATION DES DOSSIERS
# =======================================================

DIR_SRC = src
DIR_OBJ = build
DIR_MOD = $(DIR_OBJ)/mod
DIR_LIB = lib
DIR_NUMERICS = $(DIR_SRC)/numerics
DIR_PHYSICS = $(DIR_SRC)/physics
DIR_IO = $(DIR_SRC)/io
DIR_ODEPACK = $(DIR_LIB)/odepack

# =======================================================
# INFORMATIONS GÉNÉRALES SUR LA COMPILATION
# =======================================================

TARGET = main

# =======================================================
# DÉTECTION AUTOMATIQUE DES FLAGS MPI
# =======================================================

SHOW_MPI := $(shell mpif90 --showme:compile) $(shell mpif90 --showme:link)
CFLAGS_MPI := $(shell printf '%s\n' "$(SHOW_MPI)" | tr ' ' '\n' | grep '^-I')
LFLAGS_MPI := $(shell printf '%s\n' "$(SHOW_MPI)" | tr ' ' '\n' | grep '^-L')
LIBS_MPI   := $(shell printf '%s\n' "$(SHOW_MPI)" | tr ' ' '\n' | grep '^-l')
LIBDIR_MPI := $(shell printf '%s\n' "$(SHOW_MPI)" | tr ' ' '\n' | grep '^-L' | head -n1 | sed 's/^-L//')

# =======================================================
# DÉTECTION AUTOMATIQUE DES FLAGS HDF5
# =======================================================

SHOW_HDF5 := $(shell h5fc -show)
CFLAGS_HDF5 := $(shell printf '%s\n' "$(SHOW_HDF5)" | tr ' ' '\n' | grep '^-I')
LFLAGS_HDF5 := $(shell printf '%s\n' "$(SHOW_HDF5)" | tr ' ' '\n' | grep '^-L')
LIBS_HDF5 := $(shell printf '%s\n' "$(SHOW_HDF5)" | tr ' ' '\n' | grep -E '^-l|^/.+\.(a|so|so\.[0-9.]+|dylib)$$')
LIBDIR_HDF5 := $(shell printf '%s\n' "$(SHOW_HDF5)" | tr ' ' '\n' | grep '^-L' | head -n1 | sed 's/^-L//')

# =======================================================
# OPTIONS DE COMPILATION
# =======================================================

FFLAGS = -O3 \
         -fopenmp \
         -J$(DIR_MOD) \
         -I$(DIR_MOD) \
         -I$(DIR_SRC) \
         -I$(DIR_NUMERICS) \
         -I$(DIR_PHYSICS) \
         -I$(DIR_IO) \
         $(CFLAGS_MPI) \
         $(CFLAGS_HDF5)
#         -ffree-line-length-none

FFLAGS_LEGACY = $(FFLAGS) -std=legacy

LDFLAGS = -fopenmp \
          $(LFLAGS_MPI) \
          $(LFLAGS_HDF5)

ifneq ($(LIBDIR_MPI),)
LDFLAGS += -Wl,-rpath,$(LIBDIR_MPI)
endif

ifneq ($(LIBDIR_HDF5),)
LDFLAGS += -Wl,-rpath,$(LIBDIR_HDF5)
endif

# =======================================================
# LISTE DES FICHIERS SOURCE
# =======================================================

SRC_ODEPACK = $(DIR_ODEPACK)/odepack_common.f90 \
              $(DIR_ODEPACK)/odepack_interface.f90 \
              $(DIR_ODEPACK)/odepack_mod.f90 \
              $(DIR_ODEPACK)/odepack_sub1.f \
              $(DIR_ODEPACK)/odepack_sub2.f \
              $(DIR_ODEPACK)/odepack.f \
              $(wildcard $(DIR_ODEPACK)/lapack/*.f)

SRC_MAIN = $(DIR_NUMERICS)/m_numerical_init.f90 \
           $(DIR_PHYSICS)/m_physical_init.f90 \
           $(DIR_IO)/m_io.f90 \
           $(DIR_PHYSICS)/m_mhd.f90 \
           $(DIR_NUMERICS)/m_integrator.f90 \
           $(DIR_PHYSICS)/m_invariants.f90 \
           $(DIR_PHYSICS)/m_mach_numbers.f90 \
           $(DIR_PHYSICS)/m_critical_points.f90 \
           $(DIR_SRC)/m_mhd_solver.f90 \
           $(DIR_SRC)/m_solver.f90 \
           $(DIR_SRC)/main.f90

SRC = $(SRC_ODEPACK) $(SRC_MAIN)

# =======================================================
# LISTE DES FICHIERS OBJETS
# =======================================================

OBJ = $(SRC:%=$(DIR_OBJ)/%.o)

# =======================================================
# CIBLES UTILITAIRES
# =======================================================

.PHONY: all run clean

all: $(TARGET)

run: $(TARGET)
	./$(TARGET)

# =======================================================
# ÉDITION DE LIENS : FABRIQUER L’EXÉCUTABLE FINAL
# =======================================================

$(TARGET): $(OBJ)
	$(FC) $(OBJ) $(LDFLAGS) $(LIBS_HDF5) $(LIBS_MPI) -o $@

# =======================================================
# DÉPENDANCES EXPLICITES ENTRE MODULES FORTRAN
# =======================================================

$(DIR_OBJ)/$(DIR_NUMERICS)/m_integrator.f90.o: $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_mod.f90.o \
                                               $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_common.f90.o \
                                               $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_interface.f90.o

$(DIR_OBJ)/$(DIR_ODEPACK)/odepack_mod.f90.o: $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_common.f90.o \
                                             $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_interface.f90.o

$(DIR_OBJ)/$(DIR_ODEPACK)/odepack_sub1.f.o: $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_common.f90.o \
                                            $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_interface.f90.o \
                                            $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_mod.f90.o

$(DIR_OBJ)/$(DIR_ODEPACK)/odepack_sub2.f.o: $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_common.f90.o \
                                            $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_interface.f90.o \
                                            $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_mod.f90.o

$(DIR_OBJ)/$(DIR_ODEPACK)/odepack.f.o: $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_common.f90.o \
                                       $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_interface.f90.o \
                                       $(DIR_OBJ)/$(DIR_ODEPACK)/odepack_mod.f90.o

# =======================================================
# RÈGLE DE COMPILATION : .f90 --> .o
# =======================================================

$(DIR_OBJ)/%.f90.o: %.f90
	@mkdir -p $(dir $@) $(DIR_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

# =======================================================
# RÈGLE DE COMPILATION : .f --> .o
# =======================================================

$(DIR_OBJ)/%.f.o: %.f
	@mkdir -p $(dir $@) $(DIR_MOD)
	$(FC) $(FFLAGS_LEGACY) -c $< -o $@

# =======================================================
# CIBLE DE NETTOYAGE
# =======================================================

clean:
	rm -rf $(TARGET) $(DIR_OBJ)