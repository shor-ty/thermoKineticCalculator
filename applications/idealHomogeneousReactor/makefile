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

PROGRAM=idealHomogeneousReactor
COMPILER=g++
MAKE_DIR=mkdir -p
RM_DIR=rm -rf
SRC_PATH=../../src/gcc/lnInclude
LIB_PATH=../../platforms/libs/
DIR_APP=../../platforms/bin/

#------------------------------------------------------------------------------

build: pre
	$(shell echo $(APP_PATH))
	$(COMPILER) $(CPPFLAGS) -I$(SRC_PATH) -L$(LIB_PATH) $(addsuffix .cpp, $(PROGRAM)) -lthermoKinetics -o $(addprefix $(DIR_APP), $(PROGRAM))


pre:
	$(shell $(MAKE_DIR) $(DIR_APP))


rebuild: clean build

clean:
	$(shell $(RM_DIR) $(DIR_APP))


#------------------------------------------------------------------------------
