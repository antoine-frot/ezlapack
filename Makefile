# Compiler and Flags
FC = gfortran
FLAGS = -O3 -g -Wall
LIB = -lblas -llapack

# Directories
SRC_DIR = src
TEST_DIR = test
BIN_DIR = bin

# Files
MOD_FILES = $(wildcard $(SRC_DIR)/*.f90)
TEST_FILES = $(wildcard $(TEST_DIR)/*.f90)
OBJ_FILES = $(MOD_FILES:$(SRC_DIR)/%.f90=$(BIN_DIR)/%.o)
TEST_OBJ_FILES = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%.o)
EXEC = $(BIN_DIR)/main_prog

# Rules
all: $(EXEC)

$(EXEC): $(OBJ_FILES)
	$(FC) $(FLAGS) -o $@ $? $(LIB)
	@rm $(MOD_FILES:$(MOD_DIR)/%.f90=%.mod)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -c $< -o $@

$(BIN_DIR)/%.o: $(TEST_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -c $< -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

.PHONY: clean
clean:
	rm -f $(BIN_DIR)/*.o $(BIN_DIR)/*.mod $(EXEC)
	rmdir $(BIN_DIR)
