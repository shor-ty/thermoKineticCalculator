#
# Tobias Holzmann
# Februar 2017
#
# Description
#     This makefile compiles the libraries and the application for generating 
#     steady-state flamelets that can be used with the flameletHolzmannFoam
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

PROGRAM=SolveAx
OBJECTDIR=obj
MAKE_DIR=mkdir -p
CLEAN_DIR=rm -rf
WORK_DIR=$(PWD)
COMPILER=g++
PATH=../../src/

# Compiling mode: DEBUG | OPT
COMPILE_MODE=DEBUG

DIRECTORIES=numerics/ODE \
            numerics/LUDecompose \
            typedef \
            tensors/tensor \
            tensors/matrix \
            tensors/vector


CPPFLAGS=-Wall \
		 -Wextra \
         -std=c++11 \
		 -D_GLIBCXX_DEBUG \
         -g \
         -pedantic-errors \
         -ggdb \
         -Wno-unused-parameter \
         -Wno-unused-variable 


DIRECTORIES_LIBS=$(addprefix $(PATH), $(DIRECTORIES))

SOURCE_FILES=$(shell find $(DIRECTORIES_LIBS) -name '*.cpp')

OBJECT_FILES=$(addprefix obj/, $(addsuffix .o, $(basename $(SOURCE_FILES))))

INCLUDE_DIR=$(addprefix -I,$(DIRECTORIES))


# Make  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

show:
	@echo "Foo"
	@echo "me"
	@echo $(MAKE_DIR)

pre:
	$(MAKE_DIR) $(addprefix obj/,$(DIRECTORIES))

clean:
	$(CLEAN_DIR) obj

build : $(OBJECT_FILES) 
	$(COMPILER) $(CPPFLAGS) $(INCLUDE_DIR) -o $(PROGRAM) $(addsuffix .cpp, $(PROGRAM)) $^

obj/%.o: %.cpp | obj
	$(COMPILER) $(CPPFLAGS) $(INCLUDE_DIR) -c $< -o $@
			
help: 
	@echo "Targets:"
	@echo "build        - build all libraries"
	@echo "clean        - clean all libraries"
	@echo "pre          - build directory tree for libraries"


