#=================================================================
#
#	Author: Sebastian Daberdaku
#	Project: InterfaceDescriptors	
#
#	This makefile searches SRC_DIR and its subdirectories
#	recursively for source files and builds them.
#
#	Object and dependency files can be placed in separate 
#	directories.
#
#=================================================================
CC := g++ -m64
MKDIR_P := @mkdir -p 

GPPOPT := -c -std=c++11 -MMD -O3 -march=native -funroll-loops -fopenmp
LDOPT := -static -std=c++11 -O3 -march=native -fopenmp 
LIBS := -lboost_program_options -lboost_regex -lboost_system -lboost_filesystem

BIN_DIR := bin
OBJ_DIR := bin/obj
DEP_DIR := bin/dep
OUT_DIR := bin/output
SRC_DIR := src

# Defines: just add the defines to this variable
DEFS := -D _GLIBCXX_PARALLEL -D NDEBUG -D STATIC_GM #-D PRINT_TEST_G
#-D PROGRESS_BAR -D INT_STAT# -D PRINT_TEST_G -D PRINT_TEST_Z#-D TIGHT_PACKING #-D NO_OUTPUT_TEST
#-D  #-D WORD32 -D RANGECHECK -D NO_OUTPUT_TEST -D PRINT_TEST_G -D PRINT_TEST_Z 

# Make does not offer a recursive wildcard function, so here's one:
rwildcard = $(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))

# Recursively find all source files in SRC_DIR
SOURCES := $(call rwildcard,$(SRC_DIR)/,*.cpp)
OBJECTS := $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(SOURCES:$(SRC_DIR)/%.cpp=$(DEP_DIR)/%.d)
EXECUTABLE := $(BIN_DIR)/test

.PHONY: all clean	

all: $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDOPT) $^ -o $@ $(LIBS)
	$(MKDIR_P) $(OUT_DIR)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(MKDIR_P) $(@D) $(@D:$(OBJ_DIR)%=$(DEP_DIR)%)
	$(CC) $(GPPOPT) $< -o $@ -MF $(@:$(OBJ_DIR)%.o=$(DEP_DIR)%.d) $(DEFS)
	
clean: 
	rm -rf $(EXECUTABLE) $(OBJ_DIR) $(DEP_DIR)

-include $(DEPENDENCIES)
	
