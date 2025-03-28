###########################################################

## USER SPECIFIC DIRECTORIES ##

##########################################################

## CC COMPILER OPTIONS ##

# CC compiler options:
CC=/usr/bin/g++
CC_FLAGS= -O3
CC_LIBS= -lstdc++fs

##########################################################

LFLAGS= -lm

##########################################################

## Project file structure ##

# Source file directory:
SRC_DIR = src

# Object file directory:
OBJ_DIR = bin

# Include header file diretory:
INC_DIR = include

##########################################################

## Make variables ##

# Target executable name:
# make samples
EXE = compress

# run dynamics
#EXE = runNVE
#EXE = runNVT

# testing
#EXE = testInteraction

# Object files:
OBJS = $(OBJ_DIR)/$(EXE).o $(OBJ_DIR)/simSoft.o $(OBJ_DIR)/FIRE.o

##########################################################

## Compile ##

# Ensure bin directory exists
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	
# Link C++ compiled object files to target executable:
$(EXE) : $(OBJS)
	$(CC) $(OBJS) -o $@ $(CC_LIBS)

# Compile main .cpp file to object files:
$(OBJ_DIR)/$(EXE).o : $(EXE).cpp
	$(CC) $(CC_FLAGS) -c $< -o $@

# Compile C++ source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CC_FLAGS) -c $< -o $@

# Clean objects in object directory.
clean:
	$(RM) bin/* *.o $(EXE)

