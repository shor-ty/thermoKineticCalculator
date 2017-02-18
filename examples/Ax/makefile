#
# Tobias Holzmann
# Februar 2017
#
# Description
#     This makefile compiles the libraries and the application for generating 
#     steady-state flamelets that can be used with the flameletHolzmannFoam
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

PROGRAM=solveMatrixSystem
OBJECTDIR=obj
MAKE_DIR=mkdir -p
CLEAN_DIR=rm -rfv
WORK_DIR=$(PWD)
COMPILER=g++
SRC_PATH=../../src/

# Compiling mode: DEBUG | OPT
COMPILE_MODE=DEBUG

CLASSES=numerics/LUDecompose \
		tensors/matrix \
		tensors/vector \
		tensors/tensor \
        typedef


CPPFLAGS=-Wall \
		 -Wextra \
         -std=c++11 \
		 -D_GLIBCXX_DEBUG \
         -g \
         -pedantic-errors \
         -ggdb \
         -Wno-unused-parameter \
         -Wno-unused-variable 


SOURCE_FILES=$(shell find $(CLASSES) -name '*.cpp')

OBJECT_FILES=$(addprefix obj/, $(addsuffix .o, $(basename $(SOURCE_FILES))))

INCLUDE_DIR=$(addprefix -I,$(CLASSES))


# Make  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

pre:
	$(MAKE_DIR) $(CLASSES)
	$(MAKE_DIR) $(addprefix obj/,$(CLASSES))
	./link.sh N

clean:
	$(CLEAN_DIR) obj $(PROGRAM) $(CLASSES) tensors numerics

build : $(OBJECT_FILES) 
	$(COMPILER) $(CPPFLAGS) $(INCLUDE_DIR) -o $(PROGRAM) $(addsuffix .cpp, $(PROGRAM)) $^

obj/%.o: %.cpp | obj
	$(COMPILER) $(CPPFLAGS) $(INCLUDE_DIR) -c $< -o $@
			
help: 
	@echo "Targets:"
	@echo "build        - build all libraries"
	@echo "clean        - clean all libraries"
	@echo "pre          - build directory tree for libraries"


