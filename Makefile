# configure generated Makefile
#
# Makefile for Packmol: Read the comments if you have some
#                       problem while compiling.
#
# You may use the ./configure script to search automatically for
# some fortran compiler.
#
# This make file will try to compile packmol with the default
# fortran compiler, defined by the FC directive. For doing this,
# just type
#
#          make 
#
# If you want to compile with some specific fortran compiler,
# you must change the line below to the path of your fortran compiler.
#
FORTRAN = /usr/bin/gfortran
#
# Change "AUTO" to the fortran command you want. After changing
# this line, compile the package.
#
# Change the flags of the compilation if you want:
#
FLAGS= -O3 -ffast-math 
 
###################################################################
#                                                                 #
# Generally no modifications are required after this.             #
#                                                                 #
###################################################################
#
# Get the default fortran compiler
#
 
ifeq ($(FORTRAN),AUTO)
FORTRAN = $(FC)
endif 
#
# Files required
#
oall = cenmass.o \
       gencan.o \
       pgencan.o \
       initial.o \
       title.o \
       io.o \
	 fgcommon.o \
       packmol.o \
       polartocart.o \
       heuristics.o \
       random.o \
       sizes.o \
       usegencan.o \
       molpa.o \
       feasy.o \
       geasy.o
#
# Linking 
#
all : $(oall)
	@echo " ------------------------------------------------------ " 
	@echo " Compiling packmol with $(FORTRAN) " 
	@echo " Flags: $(FLAGS) " 
	@echo " ------------------------------------------------------ " 
	@$(FORTRAN) -o packmol $(oall) $(FLAGS) 
	@\rm -f *.mod *.o
	@echo " ------------------------------------------------------ " 
	@echo " Packmol succesfully built." 
	@echo " ------------------------------------------------------ " 
#
# Compiling with flags for development
#
devel : $(oall) 
	@echo " ------------------------------------------------------ " 
	@echo " Compiling packmol with $(FORTRAN) " 
	@echo " Flags: -Wunused -fcheck=bounds"
	@echo " ------------------------------------------------------ "
	@$(FORTRAN) -o packmol $(oall) -Wunused -fcheck=bounds 
	@echo " ------------------------------------------------------ " 
	@echo " Packmol succesfully built. " 
	@echo " ------------------------------------------------------ " 
#
# Modules
#
modules = sizes.o molpa.o usegencan.o
sizes.o : sizes.f90 
	@$(FORTRAN) $(FLAGS) -c sizes.f90
molpa.o : molpa.f90 sizes.o
	@$(FORTRAN) $(FLAGS) -c molpa.f90
usegencan.o : usegencan.f90 
	@$(FORTRAN) $(FLAGS) -c usegencan.f90
#
# Code compiled only for all versions
#
cenmass.o : cenmass.f90 $(modules)
	@$(FORTRAN) $(FLAGS) -c cenmass.f90
initial.o : initial.f90 $(modules)
	@$(FORTRAN) $(FLAGS) -c initial.f90
title.o : title.f90 
	@$(FORTRAN) $(FLAGS) -c title.f90
io.o : io.f90  $(modules)
	@$(FORTRAN) $(FLAGS) -c io.f90
fgcommon.o : fgcommon.f90 $(modules)
	@$(FORTRAN) $(FLAGS) -c fgcommon.f90
packmol.o : packmol.f90 $(modules)
	@$(FORTRAN) $(FLAGS) -c packmol.f90
polartocart.o : polartocart.f90 $(modules)   
	@$(FORTRAN) $(FLAGS) -c polartocart.f90
heuristics.o : heuristics.f90 $(modules)   
	@$(FORTRAN) $(FLAGS) -c heuristics.f90
random.o : random.f90 $(modules)   
	@$(FORTRAN) $(FLAGS) -c random.f90
pgencan.o : pgencan.f90 $(modules)
	@$(FORTRAN) $(FLAGS) -c pgencan.f90
gencan.o : gencan.f
	@$(FORTRAN) $(FLAGS) -c gencan.f 
feasy.o : feasy.f90 $(modules)   
	@$(FORTRAN) $(FLAGS) -c feasy.f90
geasy.o : geasy.f90 $(modules)   
	@$(FORTRAN) $(FLAGS) -c geasy.f90
#
# Clean build files
#
clean: 
	@\rm -f ./*.o ./*.mod 
