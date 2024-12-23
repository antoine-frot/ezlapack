# Compiler and Flags
FC = gfortran
FLAGS = -O3 -g -Wall
LIB = -lezlapack -lblas -llapack
LIB_NAME = ezlapack
PATH_LIBRARY = /usr/local/lib
PATH_MOD = /usr/local/include

# Directories
SRC_DIR = src
TEST_DIR = test
BIN_DIR = bin
MOD_DIR = module

# Source Files
MOD_FILES = $(wildcard $(MOD_DIR)/module_*.f90)
MOD_OBJ_FILES = $(MOD_FILES:$(MOD_DIR)/%.f90=$(BIN_DIR)/%.o)
MAIN_FILE = $(SRC_DIR)/main_program.f90
MAIN_OBJ_FILE = $(MAIN_FILE:$(SRC_DIR)/%.f90=$(BIN_DIR)/%.o)
TEST_FILES = $(wildcard $(TEST_DIR)/test*.f90)
TEST_OBJ_FILES = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%.o)
EXEC = $(BIN_DIR)/main_program
EXEC_TEST = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%)

# Default target
all: $(EXEC)

# Build the static library
$(BIN_DIR)/lib$(LIB_NAME).a: $(BIN_DIR)/$(LIB_NAME).o $(MOD_OBJ_FILES)
	@echo "Archiving library..."
	ar rcs $@ $^

$(BIN_DIR)/$(LIB_NAME).o: $(MOD_OBJ_FILES)
	$(FC) $(FLAGS) -J$(BIN_DIR) -o $@ -c $(MOD_DIR)/$(LIB_NAME).f90 $(LIB)

# Install the library
install: $(BIN_DIR)/lib$(LIB_NAME).a
	@echo "Installing library..."
	mkdir -p $(PATH_LIBRARY)/lib $(PATH_LIBRARY)/include
	cp $(BIN_DIR)/lib$(LIB_NAME).a $(PATH_LIBRARY)
	cp $(BIN_DIR)/*.mod $(PATH_MOD)
	ldconfig

# Build the main program
$(EXEC): $(MAIN_OBJ_FILE)
	$(FC) $(FLAGS) -o $@ $^ $(LIB)

# Build test programs
test: $(EXEC_TEST)
	@echo "Running tests..."
	@for exec in $(EXEC_TEST:$(BIN_DIR)/%=%); do \
		echo "Running $$exec..."; \
		$(BIN_DIR)/$$exec; \
	done

$(BIN_DIR)/%: $(BIN_DIR)/%.o $(MOD_OBJ_FILES)
	$(FC) $(FLAGS) -o $@ $^ $(LIB)

# Compile object files
$(BIN_DIR)/%.o: $(MOD_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(BIN_DIR) -J$(BIN_DIR)  -c $< -o $@

$(BIN_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(PATH_MOD) -c $< -o $@

$(BIN_DIR)/%.o: $(TEST_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(BIN_DIR) -c $< -o $@

# Create binary directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Clean rule
.PHONY: clean
clean:
	sudo rm -rf $(BIN_DIR)

# Phony targets
.PHONY: all test install

