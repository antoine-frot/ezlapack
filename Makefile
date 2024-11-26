# Compiler and Flags
FC = gfortran
FLAGS = -O3

# Directories
SRC_DIR = src
TEST_DIR = test
BIN_DIR = bin

# Source Files
SRC_FILES = $(SRC_DIR)/dgemm_wrapper.f90
TEST_FILES = $(TEST_DIR)/test_matrice_mult_d.f90

# Object Files
OBJ_FILES = $(BIN_DIR)/dgemm_wrapper.o

# Executable
EXEC = $(BIN_DIR)/test_matrices

# Default target
.PHONY: all clean

all: $(EXEC)

# Rule for creating the executable
$(EXEC): $(OBJ_FILES) $(TEST_FILES)
	$(FC) $(FLAGS) -o $@ $(TEST_FILES) $(OBJ_FILES) -lblas
	rm $(OBJ_FILES)

# Rule for compiling dgemm_wrapper.f90
$(BIN_DIR)/dgemm_wrapper.o: $(SRC_FILES) | $(BIN_DIR)
	$(FC) $(FLAGS) -c $< -o $@

# Ensure bin directory exists
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Clean up build files
clean:
	rm -rf $(BIN_DIR)

