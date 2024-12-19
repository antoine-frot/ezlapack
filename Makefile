# Compiler and Flags
FC = gfortran
FLAGS = -O3 -g -Wall
LIB = -lblas -llapack

# Directories
SRC_DIR = src
TEST_DIR = test
BIN_DIR = bin
MOD_DIR = module

# Files
MOD_FILES = $(wildcard $(MOD_DIR)/mod*.f90)
MOD_OBJ_FILES = $(MOD_FILES:$(MOD_DIR)/%.f90=$(BIN_DIR)/%.o)

EXEC = $(BIN_DIR)/main_program
MAIN_FILE = $(EXEC:$(BIN_DIR)/%=$(SRC_DIR)/%.f90)
MAIN_OBJ_FILE = $(EXEC:$(BIN_DIR)/%=$(BIN_DIR)/%.o)

TEST_FILES = $(wildcard $(TEST_DIR)/test*.f90)
EXEC_TEST = $(TEST_FILES:$(TEST_DIR)/test%.f90=$(BIN_DIR)/test%)
TEST_OBJ_FILES = $(TEST_FILES:$(TEST_DIR)/test%.f90=$(BIN_DIR)/test%.o)

# Rules

all: $(EXEC)

test: $(EXEC_TEST)
	$(EXEC_TEST)

$(EXEC): $(MOD_OBJ_FILES) $(MAIN_OBJ_FILE)
	$(FC) $(FLAGS) -o $@ $^ $(LIB)

$(EXEC_TEST): $(MOD_OBJ_FILES) $(TEST_OBJ_FILE)
	$(FC) $(FLAGS) -o $@ $^ $(LIB)

# Explicitly compile the module file first
$(BIN_DIR)/mod%.o: $(MOD_DIR)/mod%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@
	@mv *.mod $(BIN_DIR) 2>/dev/null || true

# Compile the main program after modules
$(MAIN_OBJ_FILE): $(MAIN_FILE) $(MOD_MOD_FILES)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@

# Compile test files, ensuring module files are available in $(BIN_DIR)
$(TEST_OBJ_FILE): $(TEST_FILES) $(MOD_MOD_FILES) | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@

# Create binary directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

.PHONY: clean
clean:
	rm -rf $(BIN_DIR)/
	
