# Compiler and Flags
FC = gfortran
FLAGS = -O3 -g -Wall
LIB = -lblas -llapack

# Directories
SRC_DIR = src
TEST_DIR = test
BIN_DIR = bin

# Files
MOD_FILES = $(wildcard $(SRC_DIR)/mod_*.f90)
MOD_OBJ_FILES = $(MOD_FILES:$(SRC_DIR)/%.f90=$(BIN_DIR)/%.o)
MOD_MOD_FILES = $(MOD_FILES:$(SRC_DIR)/%.f90=$(BIN_DIR)/%.mod)
EXEC = $(BIN_DIR)/main_program
MAIN_FILE = $(EXEC:$(BIN_DIR)/%=$(SRC_DIR)/%.f90)
MAIN_OBJ_FILE = $(EXEC:$(BIN_DIR)/%=$(BIN_DIR)/%.o)

TEST_FILES = $(wildcard $(TEST_DIR)/*.f90)
TEST_OBJ_FILES = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%.o)

# Rules
all: $(EXEC)

$(EXEC): $(MOD_OBJ_FILES) $(MAIN_OBJ_FILE)
	$(FC) $(FLAGS) -o $@ $^ $(LIB)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -c $< -o $@
	@mv *.mod $(BIN_DIR) 2>/dev/null || true

$(MAIN_OBJ_FILE): $(MAIN_FILE) $(MOD_MOD_FILES)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@

$(BIN_DIR)/%.o: $(TEST_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

.PHONY: clean
clean:
	rm -f $(BIN_DIR)/*.o $(BIN_DIR)/*.mod $(EXEC)
	rmdir $(BIN_DIR)

