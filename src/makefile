WORK_DIR=home/shorty/OpenFOAM/development/flameletcreator/src

INC_DIR=-I/$(WORK_DIR)/chemistry \
		-I/$(WORK_DIR)/typedef \
		-I/$(WORK_DIR)/stringManipulator

SOURCE_FILES=automaticFlameletCreator.cpp \
			 chemistry/chemistry.cpp \
			 chemistry/chemistryReader.cpp \
			 chemistry/chemistryData.cpp \
			 typedef/typedef.cpp\
			 stringManipulator/stringManipulator.cpp

CPPFLAGS=-g -Wall -Wextra -pedantic-errors -std=c++11

automaticFlameltCreator : automaticFlameletCreator.cpp
	g++ $(CPPFLAGS) $(INC_DIR) -o automaticFlameletCreator $(SOURCE_FILES)
