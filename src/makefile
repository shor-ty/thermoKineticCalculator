WORK_DIR=$(PWD)

$(OUT_O_DIR)/%.o: %.cpp
	        @mkdir -p $(@D)
	        @$(CC) -c $< -o $@

INC_DIR=-I$(WORK_DIR)/chemistry \
		-I$(WORK_DIR)/thermo \
		-I$(WORK_DIR)/transport \
		-I$(WORK_DIR)/properties \
		-I$(WORK_DIR)/mixtureFraction \
		-I$(WORK_DIR)/constants \
		-I$(WORK_DIR)/numerics \
		-I$(WORK_DIR)/numerics/ODE \
		-I$(WORK_DIR)/numerics/jacobian \
		-I$(WORK_DIR)/numerics/LUDecompose \
		-I$(WORK_DIR)/numerics/seulex \
		-I$(WORK_DIR)/numerics/stepStatus \
		-I$(WORK_DIR)/typedef \
		-I$(WORK_DIR)/interpreter \
		-I$(WORK_DIR)/tensors/tensor \
		-I$(WORK_DIR)/tensors/matrix \
		-I$(WORK_DIR)/tensors/vector \
		-I$(WORK_DIR)/tensors/scalar \
		-I$(WORK_DIR)/tensors/dimensions \
		-I$(WORK_DIR)/stringManipulator



SOURCE_FILES=automaticFlameletCreator.cpp \
			 typedef/typedef.cpp \
			 numerics/LUDecompose/LUDecompose.cpp \
			 tensors/tensor/tensor.cpp \
			 tensors/matrix/matrix.cpp \
			 tensors/vector/vector.cpp \
			 tensors/scalar/scalar.cpp \
			 tensors/dimensions/dimensions.cpp
#			 numerics/seulex/seulex.cpp \
			 numerics/numerics.cpp \
			 numerics/stepStatus/stepStatus.cpp \
			 numerics/jacobian/jacobian.cpp \
			 numerics/ODE/ODE.cpp \
			 chemistry/chemistry.cpp \
			 chemistry/chemistryReader.cpp \
			 chemistry/chemistryData.cpp \
			 chemistry/chemistryCalc.cpp \
			 thermo/thermo.cpp \
			 thermo/thermoReader.cpp \
			 thermo/thermoData.cpp \
			 thermo/thermoCalc.cpp \
			 transport/transport.cpp \
			 transport/transportReader.cpp \
			 transport/transportData.cpp \
			 transport/transportCalc.cpp \
			 properties/properties.cpp \
			 properties/propertiesReader.cpp \
			 mixtureFraction/mixtureFraction.cpp \
			 interpreter/interpreter.cpp \
			 stringManipulator/stringManipulator.cpp


#SRC_FILES=mixtureFraction/numerics.cpp

CPPFLAGS=-g -Wall -Wextra -pedantic-errors -std=c++11 -Wno-unused-parameter -Wno-unused-variable -D_GLIBCXX_DEBUG

OBJ := $(addprefix obj/,$(addsuffix .o,$(basename $(SOURCE_FILES))))


obj/%.o: %.cpp | obj
	    g++ $(CPPFLAGS) $(INC_DIR) -c $< -o $@

automaticFlameletCreator : $(OBJ)
	g++ $(CPPFLAGS) $(INC_DIR) -o $@ $^


#automaticFlameltCreator : automaticFlameletCreator.cpp
#	g++ $(CPPFLAGS) $(INC_DIR) $(OBJ) -o automaticFlameletCreator 

#CC=g++
#CPPFLAGS=-g -Wall -Wextra -pedantic-errors -std=c++11 $(INC_DIR)
# LDFLAGS=
# OBJECTS=$(SOURCES:.cpp=.o)
#	EXECUTABLE=automaticFlameletCreator
#
#all: $(SOURCE_FILES) $(EXECUTABLE)
#	    
#$(EXECUTABLE): $(OBJECTS) 
#	    $(CC) $(LDFLAGS) $(OBJECTS) -o $@
#
#.cpp.o:
#	    $(CC) $(CFLAGS) $< -o $@
