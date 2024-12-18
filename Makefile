# Compiler and Flags
FC = gfortran
FLAGS = -O3 -g -Wall
LIB = -lblas -llapack

# Directories
SRC_DIR = src
TEST_DIR = test
BIN_DIR = bin

# Files
MOD_FILES = $(wildcard $(SRC_DIR)/mod*.f90)
MOD_OBJ_FILES = $(MOD_FILES:$(SRC_DIR)/%.f90=$(BIN_DIR)/%.o)

EXEC = $(BIN_DIR)/main_program
MAIN_FILE = $(EXEC:$(BIN_DIR)/%=$(SRC_DIR)/%.f90)
MAIN_OBJ_FILE = $(EXEC:$(BIN_DIR)/%=$(BIN_DIR)/%.o)

TEST_FILES = $(wildcard $(TEST_DIR)/*.f90)
TEST_OBJ_FILES = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%.o)

# Rules

all: $(EXEC)

$(EXEC): $(MOD_OBJ_FILES) $(MAIN_OBJ_FILE)
	$(FC) $(FLAGS) -o $@ $^ $(LIB)

# Explicitly compile the module file first
$(BIN_DIR)/mod%.o: $(SRC_DIR)/mod%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@
	@mv *.mod $(BIN_DIR) 2>/dev/null || true

# Compile the main program after modules
$(MAIN_OBJ_FILE): $(MAIN_FILE) $(MOD_MOD_FILES)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@

# Compile test files, ensuring module files are available in $(BIN_DIR)
$(BIN_DIR)/%.o: $(TEST_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@

# Create binary directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

.PHONY: clean
clean:
	rm -rf $(BIN_DIR)/

