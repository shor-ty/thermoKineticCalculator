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

PROGRAM=homogeneousReactor
COMPILER=g++
MAKE_DIR=mkdir -p
SRC_PATH=../../src/gcc/lnInclude
LIB_PATH=../../platforms/libs/
DIR_APP=../../platforms/bin/

#------------------------------------------------------------------------------

build : pre
	$(shell echo $(APP_PATH))
	$(COMPILER) $(CPPFLAGS) -I$(SRC_PATH) -L$(LIB_PATH) homogeneousReactor.cpp -lthermoKinetics -o $(addprefix $(DIR_APP), $(PROGRAM))


pre:
	$(shell $(MAKE_DIR) $(DIR_APP))


#------------------------------------------------------------------------------
