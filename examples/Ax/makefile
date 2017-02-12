WORK_DIR=$(PWD)/../../src

$(OUT_O_DIR)/%.o: %.cpp
	        @mkdir -p $(@D)
	        @$(CC) -c $< -o $@

INC_DIR=-I$(WORK_DIR)/numerics/LUDecompose \
		-I$(WORK_DIR)/typedef \
		-I$(WORK_DIR)/tensors/tensor \
		-I$(WORK_DIR)/tensors/matrix \
		-I$(WORK_DIR)/tensors/vector \
		-I$(WORK_DIR)/tensors/scalar



SOURCE_FILES=solveMatrixSystem.cpp \
			 typedef/typedef.cpp \
			 numerics/LUDecompose/LUDecompose.cpp \
			 tensors/tensor/tensor.cpp \
			 tensors/matrix/matrix.cpp \
			 tensors/vector/vector.cpp \
			 tensors/scalar/scalar.cpp \
			 tensors/dimensions/dimensions.cpp


CPPFLAGS=-g -Wall -Wextra -pedantic-errors -std=c++11 -Wno-unused-parameter -Wno-unused-variable -D_GLIBCXX_DEBUG

OBJ := $(addprefix obj/,$(addsuffix .o,$(basename $(SOURCE_FILES))))


obj/%.o: %.cpp | obj
	    g++ $(CPPFLAGS) $(INC_DIR) -c $< -o $@

solveMatrixSystem : $(OBJ)
	g++ $(CPPFLAGS) $(INC_DIR) -o $@ $^
