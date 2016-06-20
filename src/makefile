WORK_DIR=$(PWD)

INC_DIR=-I$(WORK_DIR)/chemistry \
		-I$(WORK_DIR)/thermo \
		-I$(WORK_DIR)/transport \
		-I$(WORK_DIR)/properties \
		-I$(WORK_DIR)/mixtureFraction \
		-I$(WORK_DIR)/constants \
		-I$(WORK_DIR)/typedef \
		-I$(WORK_DIR)/stringManipulator \
		-I$(WORK_DIR)/thirdParty



SOURCE_FILES=automaticFlameletCreator.cpp \
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
			 typedef/typedef.cpp\
			 stringManipulator/stringManipulator.cpp


SRC_FILES=mixtureFraction/numerics.cpp \
		  mixtureFraction/interprete

CPPFLAGS=-g -Wall -Wextra -pedantic-errors -std=c++11 -Wno-unused-parameter -Wno-unused-variable

OBJ := $(addprefix obj/,$(addsuffix .o,$(basename $(SOURCE_FILES))))


obj/%.o: %.cpp | obj
	    g++ $(CPPFLAGS) $(INC_DIR) -c $< -o $@

automaticFlameletCreator : $(OBJ) | mixtureFraction/numerics.cpp mixtureFraction/interprete.cpp
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
