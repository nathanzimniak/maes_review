# =======================================================
#  Projet MAES — Compilation du code principal + LSODA
# =======================================================

# Nom de l’exécutable final à produire
TARGET = maes

# Répertoires utilisés pour les sources, modules et objets compilés
SRC_DIR = src
ODEPACK_DIR = lib/odepack
OBJ_DIR = build/obj
MOD_DIR = build/mod
NUMERICS_DIR = $(SRC_DIR)/numerics
PHYSICS_DIR = $(SRC_DIR)/physics
IO_DIR = $(SRC_DIR)/io

# -------------------------------------------------------
# Options de compilation
# -------------------------------------------------------

# Définition du compilateur Fortran
FC = gfortran

# Flags reliés à HDF5
HDF5_CFLAGS := $(shell h5fc -show | tr ' ' '\n' | grep '^-I')
HDF5_LFLAGS := $(shell h5fc -show | tr ' ' '\n' | grep '^-L\|^-l')
HDF5_LIBDIR := $(shell h5fc -show | tr ' ' '\n' | grep '^-L' | head -n1 | sed 's/^-L//')

# Options de compilation en mode release (exécution rapide)
FFLAGS = -O3 \
         -J$(MOD_DIR) -I$(MOD_DIR) \
         -fopenmp \
         -ffree-line-length-none \
         -fno-backtrace -ffpe-summary=none \
         -w \
         $(HDF5_CFLAGS)

## Options de compilation en mode debug (compilation détaillée)
#FFLAGS = -O0 \
#         -g -Wall -Wextra -Wconversion -Wimplicit-interface \
#         -fcheck=all -fbacktrace \
#         -J$(MOD_DIR) -I$(MOD_DIR) \
#         -fopenmp \
#         -ffree-line-length-none \
#         $(HDF5_CFLAGS)

# Options passées à l’édition de liens (utilisées pour assembler tous les fichiers objets en un exécutable)
LDFLAGS = -fopenmp \
          $(HDF5_LFLAGS) \
          -Wl,-rpath,$(HDF5_LIBDIR)

# -------------------------------------------------------
# Liste des fichiers sources à compiler
# -------------------------------------------------------

# Modules du cœur scientifique (ordonnés selon leurs dépendances)
SRC_CORE = \
  $(NUMERICS_DIR)/m_numerical_init.f90 \
  $(PHYSICS_DIR)/m_physical_init.f90 \
  $(IO_DIR)/m_io.f90 \
  $(PHYSICS_DIR)/m_mhd.f90 \
  $(NUMERICS_DIR)/m_integrator.f90 \
  $(PHYSICS_DIR)/m_invariants.f90 \
  $(PHYSICS_DIR)/m_mach_numbers.f90 \
  $(PHYSICS_DIR)/m_critical_points.f90 \
  $(SRC_DIR)/m_mhd_solver.f90 \
  $(SRC_DIR)/m_solver.f90 \

# Programme principal
SRC_MAIN = $(SRC_DIR)/main.f90

# Fichiers ODEPACK en Fortran 90 (interface moderne thread-safe de LSODA)
SRC_ODEPACK_F90 = \
  $(ODEPACK_DIR)/odepack_common.f90 \
  $(ODEPACK_DIR)/odepack_interface.f90 \
  $(ODEPACK_DIR)/odepack_mod.f90

# Fichiers ODEPACK en Fortran 77 (routines originales LSODA)
SRC_ODEPACK_F = \
  $(ODEPACK_DIR)/odepack_sub1.f \
  $(ODEPACK_DIR)/odepack_sub2.f \
  $(ODEPACK_DIR)/odepack.f

# Fichiers LAPACK (toutes les routines mathématiques du sous-dossier lapack/)
SRC_LAPACK = $(wildcard $(ODEPACK_DIR)/lapack/*.f)

# Regroupement de tous les fichiers sources
SRC = $(SRC_ODEPACK_F90) $(SRC_ODEPACK_F) $(SRC_LAPACK) $(SRC_CORE) $(SRC_MAIN)

# -------------------------------------------------------
# Liste des fichiers objets correspondant à chaque source
# -------------------------------------------------------

OBJ_CORE        = $(SRC_CORE:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)
OBJ_MAIN        = $(SRC_MAIN:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)
OBJ_ODEPACK_F90 = $(SRC_ODEPACK_F90:$(ODEPACK_DIR)/%.f90=$(OBJ_DIR)/%.o)
OBJ_ODEPACK_F   = $(SRC_ODEPACK_F:$(ODEPACK_DIR)/%.f=$(OBJ_DIR)/%.o)
OBJ_LAPACK      = $(SRC_LAPACK:$(ODEPACK_DIR)/lapack/%.f=$(OBJ_DIR)/lapack/%.o)
OBJ             = $(OBJ_CORE) $(OBJ_MAIN) $(OBJ_ODEPACK_F90) $(OBJ_ODEPACK_F) $(OBJ_LAPACK)

# -------------------------------------------------------
# Édition de liens (assemblage des objets en exécutable)
# -------------------------------------------------------

# Pour construire l’exécutable $(TARGET), il faut que tous les fichiers objets listés dans $(OBJ) existent et soient à jour.
$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LDFLAGS) -o $@

# -------------------------------------------------------
# Dépendances explicites entre modules Fortran
# -------------------------------------------------------

# Ces lignes assurent que les modules ODEPACK (.mod) sont compilés avant les fichiers qui en dépendent.
$(OBJ_DIR)/numerics/m_integrator.o: $(OBJ_DIR)/odepack_mod.o $(OBJ_DIR)/odepack_common.o $(OBJ_DIR)/odepack_interface.o
$(OBJ_DIR)/odepack_mod.o:  $(OBJ_DIR)/odepack_common.o $(OBJ_DIR)/odepack_interface.o
$(OBJ_DIR)/odepack_sub1.o: $(OBJ_DIR)/odepack_common.o $(OBJ_DIR)/odepack_interface.o $(OBJ_DIR)/odepack_mod.o
$(OBJ_DIR)/odepack_sub2.o: $(OBJ_DIR)/odepack_common.o $(OBJ_DIR)/odepack_interface.o $(OBJ_DIR)/odepack_mod.o
$(OBJ_DIR)/odepack.o:      $(OBJ_DIR)/odepack_common.o $(OBJ_DIR)/odepack_interface.o $(OBJ_DIR)/odepack_mod.o

# -------------------------------------------------------
# Règles de compilation : comment produire un .o à partir d’un .f90 ou .f
# -------------------------------------------------------

# Pour compiler un fichier du répertoire src/ (.f90)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Pour compiler un fichier du répertoire odepack/ (.f90)
$(OBJ_DIR)/%.o: $(ODEPACK_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Pour compiler un fichier du répertoire odepack/ (.f)
$(OBJ_DIR)/%.o: $(ODEPACK_DIR)/%.f
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Pour compiler un fichier du sous-dossier odepack/lapack/ (.f)
$(OBJ_DIR)/lapack/%.o: $(ODEPACK_DIR)/lapack/%.f
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# -------------------------------------------------------
# Cibles utilitaires
# -------------------------------------------------------

# Cible principale : compile tout le projet (make)
.PHONY: all
all: $(TARGET)

# Cible d’exécution : lance le programme une fois compilé (make run)
.PHONY: run
run: $(TARGET)
	./$(TARGET)

# Cible de nettoyage : supprime les fichiers compilés et l’exécutable (make clean)
.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(TARGET)