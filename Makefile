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
LIB_DIR = lib

# Source Files
# MAIN PROGRAM
MAIN_FILE = $(SRC_DIR)/speed_test.f90
EXEC = $(MAIN_FILE:$(SRC_DIR)/%.f90=$(BIN_DIR)/%)
# INSTALLATION
MOD_FILES = $(wildcard $(MOD_DIR)/module_*.f90)
MOD_OBJ_FILES = $(MOD_FILES:$(MOD_DIR)/%.f90=$(BIN_DIR)/%.o)
# TESTS
TEST_FILES = $(wildcard $(TEST_DIR)/test*.f90)
EXEC_TEST = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%)

# MAIN PROGRAM
all: $(EXEC)

run: $(EXEC)
	$(EXEC)

$(EXEC): $(MAIN_FILE) | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(LIB_DIR) -o $@ $< $(LIB) 

# INSTALLATION
install_global:
	$(MAKE) install
	sudo mkdir -p $(PATH_LIBRARY) $(PATH_MOD)
	sudo cp $(LIB_DIR)/lib$(LIB_NAME).a $(PATH_LIBRARY)
	sudo cp $(LIB_DIR)/*.mod $(PATH_MOD)
	sudo ldconfig

install: $(LIB_DIR)/lib$(LIB_NAME).a | $(LIB_DIR)
	@echo "Installing library..."

$(LIB_DIR)/lib$(LIB_NAME).a: $(BIN_DIR)/$(LIB_NAME).o $(MOD_OBJ_FILES)
	@echo "Archiving library..."
	ar rcs $@ $^

$(BIN_DIR)/$(LIB_NAME).o: $(MOD_OBJ_FILES)
	$(FC) $(FLAGS) -J$(LIB_DIR) -o $@ -c $(MOD_DIR)/$(LIB_NAME).f90 $(LIB)

$(BIN_DIR)/%.o: $(MOD_DIR)/%.f90 $(LIB_DIR) | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(LIB_DIR) -J$(LIB_DIR)  -c $< -o $@

# TESTS
test: $(EXEC_TEST)
	@echo "Running tests..."
	@for exec in $(EXEC_TEST:$(BIN_DIR)/%=%); do \
		echo "Running $$exec..."; \
		$(BIN_DIR)/$$exec; \
		sleep 1; \
		echo ""; \
	done

$(BIN_DIR)/%: $(TEST_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(LIB_MOD) -o $@ $< $(LIB)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

clean:
	rm -rf $(BIN_DIR) $(LIB_DIR)

debug:
	@echo "MOD_FILES: $(MOD_FILES)"
	@echo "MOD_OBJ_FILES: $(MOD_OBJ_FILES)"

.PHONY: all clean test install install_global run

