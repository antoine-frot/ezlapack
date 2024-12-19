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

MAIN_FILE = $(SRC_DIR)/main_program.f90
MAIN_OBJ_FILE = $(BIN_DIR)/main_program.o
EXEC = $(BIN_DIR)/main_program

TEST_FILES = $(wildcard $(TEST_DIR)/test*.f90)
TEST_OBJ_FILES = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%.o)
EXEC_TEST = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%)

# Default target: build the main program
all: $(EXEC)

# Build and run the tests
test: $(EXEC_TEST)
	@echo "Running tests..."
	@echo ""
	@for exec in $(EXEC_TEST:$(BIN_DIR)/%=%); do \
		echo "Running $$exec"; \
		echo ""; \
		$(BIN_DIR)/$$exec; \
	done

# Main program build rule
$(EXEC): $(MOD_OBJ_FILES) $(MAIN_OBJ_FILE)
	$(FC) $(FLAGS) -o $@ $^ $(LIB)

# Test program build rule
$(BIN_DIR)/%: $(BIN_DIR)/%.o $(MOD_OBJ_FILES)
	$(FC) $(FLAGS) -o $@ $^ $(LIB)

# Compile test files
$(BIN_DIR)/%.o: $(TEST_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@

# Compile main program
$(MAIN_OBJ_FILE): $(MAIN_FILE) $(MOD_OBJ_FILES)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@

# Compile module files
$(BIN_DIR)/mod%.o: $(MOD_DIR)/mod%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@
	@mv *.mod $(BIN_DIR) 2>/dev/null || true

# Create binary directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Clean rule
.PHONY: clean
clean:
	rm -rf $(BIN_DIR)

# Phony targets
.PHONY: all test

