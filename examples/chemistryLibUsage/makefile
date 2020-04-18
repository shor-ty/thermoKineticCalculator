#
# Tobias Holzmann
# Februar 2017
#
# Description
#     This makefile compiles the libraries and the application for generating 
#     steady-state flamelets that can be used with the flameletHolzmannFoam
#
#------------------------------------------------------------------------------ 

include ../../src/.compilerFlags

PROGRAM=chemistryTestProg
COMPILER=g++
SRC_PATH=../../src/gcc/lnInclude
LIB_PATH=../../src/gcc/libs/


#------------------------------------------------------------------------------ 

build :
	$(COMPILER) $(CPPFLAGS) -I$(SRC_PATH) -L$(LIB_PATH) chemistryTestProg.cpp -lthermoKinetics -o $(PROGRAM)


#------------------------------------------------------------------------------ 
