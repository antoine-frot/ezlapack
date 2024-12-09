# Compiler and Flags
FC = gfortran
FLAGS = -O3 -g -Wall
LIB = -lblas -llapack

# Directories
SRC_DIR = src
TEST_DIR = test
BIN_DIR = bin

# Files
SRC_FILES = $(wildcard $(SRC_DIR)/*.f90)
TEST_FILES = $(wildcard $(TEST_DIR)/*.f90)
OBJ_FILES = $(SRC_FILES:$(SRC_DIR)/%.f90=$(BIN_DIR)/%.o)
TEST_OBJ_FILES = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%.o)
EXEC = $(BIN_DIR)/test_ez_matmul

# Rules
all: $(EXEC)

$(EXEC): $(OBJ_FILES) $(TEST_OBJ_FILES)
	$(FC) $(FLAGS) -o $@ $? $(LIB)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -c $< -o $@

$(BIN_DIR)/%.o: $(TEST_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -c $< -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

.PHONY: clean
clean:
	rm -f $(BIN_DIR)/*.o $(BIN_DIR)/*.mod $(EXEC)

