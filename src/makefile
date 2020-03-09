#------------------------------------------------------------------------------
# AFC -- MAKEFILE
#------------------------------------------------------------------------------
#
# Tobias Holzmann
# March 2020
#
# Description
#     This makefile compiles all source files and builds the libraries that 
#     can be used for your personal program for individual usage
#
#------------------------------------------------------------------------------

include .makefile

#------------------------------------------------------------------------------

PROGRAM=automaticFlameletCreator


#------------------------------------------------------------------------------
# Linux commands
#------------------------------------------------------------------------------

MAKE_DIR=mkdir -p
LINK=ln -sf
CLEAR=clear
CLEAN_DIR=rm -rfv
WORK_DIR=$(PWD)
DIR_OBJ=gcc/obj
DIR_LIB=gcc/libs
DIR_INCLUDE=gcc/lnInclude


#------------------------------------------------------------------------------
# Compiler and compiling mode: DEBUG | OPT
#------------------------------------------------------------------------------

COMPILER=g++
COMPILE_MODE=OPT


#------------------------------------------------------------------------------
# Compiling flags for DEBUG 
#------------------------------------------------------------------------------

CPPFLAGS_DEBUG=-Wall -Wextra -std=c++17 -g -pedantic-errors \
	-ggdb -Wno-unused-parameter -Wno-unused-variable -D_GLIBCXX_DEBUG


#------------------------------------------------------------------------------
# Compiling flags for OPT
#------------------------------------------------------------------------------

CPPFLAGS_OPT=-Wall -Wextra -std=c++17 -pedantic-errors -Wno-unused-parameter \
	-Wno-unused-variable -O3


#------------------------------------------------------------------------------

ifeq ($(COMPILE_MODE), DEBUG)
CPPFLAGS=$(CPPFLAGS_DEBUG)
endif

ifeq ($(COMPILE_MODE), OPT)
CPPFLAGS=$(CPPFLAGS_OPT)
endif


#------------------------------------------------------------------------------
# Stuff for symbolic links (all header files) and libraries (folder generation)
#------------------------------------------------------------------------------
HEADER_SYM_LINK=$(shell find . -type f -iname '*.hpp')
DIR_OBJECTS=$(shell find . -type d -not -path "./gcc*" | sed "s|\.||" | sed "s|^/||")



#------------------------------------------------------------------------------
# Stuff for the thermoKinetics library
#------------------------------------------------------------------------------

DIR_THERMOKINETICS=thermoKinetics

#SOURCE_FILES_THERMO=$(shell find $(DIR_THERMOKINETICS) -name '*.cpp')
SOURCE_FILES_THERMO=\
	thermoKinetics/thermo/thermo.cpp \
	thermoKinetics/thermo/thermoCalc.cpp \
	thermoKinetics/thermo/thermoData.cpp \
	thermoKinetics/thermo/thermoReader.cpp \
	manipulation/stringManipulator/stringManipulator.cpp \
	definitions/definitions.cpp

HEADER_FILES_THERMOKINETICS=$(shell find $(DIR_THERMOKINETICS) -name '*.hpp') \
	definitions/definitions.hpp \
	manipulation/stringManipulator/stringManipulator.hpp

OBJECT_FILES_THERMOKINETICS=$(addprefix gcc/obj/, $(addsuffix .o, $(basename $(SOURCE_FILES_THERMOKINETICS))))

OBJECT_FILES_THERMO=$(addprefix gcc/obj/, $(addsuffix .o, $(basename $(SOURCE_FILES_THERMO))))


#------------------------------------------------------------------------------
# Make Commands | Probably could be done much nicer - I am not a makefile guru
#------------------------------------------------------------------------------

buildThermoKinetics: pre linkHeaders libthermo


#------------------------------------------------------------------------------
pre: clear
	$(info $(PRE_MSG))
	$(info Creating directories...)
	$(shell $(MAKE_DIR) $(DIR_OBJ) $(DIR_INCLUDE) $(DIR_LIB) $(addprefix gcc/obj/, $(DIR_OBJECTS)))
	$(shell sleep 0.5)


#------------------------------------------------------------------------------
clear: 
	$(CLEAR)


#------------------------------------------------------------------------------
linkHeaders:
	$(info Linking header files...)
	$(foreach file, $(HEADER_SYM_LINK), $(info $(file) ))
	$(foreach file, $(HEADER_SYM_LINK), $(shell $(LINK) $(WORK_DIR)/$(file) $(WORK_DIR)/$(DIR_INCLUDE)))
	$(info Linking done...)
	$(info )
	$(info )


#------------------------------------------------------------------------------
libthermo: $(OBJECT_FILES_THERMO)
	$(info Create thermo library...)
	ar rvs $(DIR_LIB)/$(addsuffix .so, $(basename $@)) $(OBJECT_FILES_THERMO)
	


#------------------------------------------------------------------------------
gcc/obj/%.o : %.cpp
	$(COMPILER) $(CPPFLAGS) -I$(DIR_INCLUDE) -c $< -o $@


#------------------------------------------------------------------------------
clean: clear
	$(info $(CLEAN_MSG))
	$(info )
	$(CLEAN_DIR) gcc $(PROGRAM)


#------------------------------------------------------------------------------
rebuild: clean buildThermoKinetics


#------------------------------------------------------------------------------
post:
	$(POST_MSG)


#------------------------------------------------------------------------------
help:
	@echo "Targets:"
	@echo "build        - build all libraries"
	@echo "clean        - clean all libraries"
	@echo "pre          - build directory tree for libraries"
	@echo "rebuild      - execute clean, pre and build"


#------------------------------------------------------------------------------
