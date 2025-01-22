# Compiler and Flags
FC = gfortran
PYTHON = python3
FLAGS = -O3 -g -Wall -Wno-missing-include-dirs
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
EXEC_PYTHON = $(MAIN_FILE:$(SRC_DIR)/%.f90=$(SRC_DIR)/%.py)
# INSTALLATION
MOD_FILES = $(wildcard $(MOD_DIR)/*.f90)
MOD_OBJ_FILES = $(MOD_FILES:$(MOD_DIR)/%.f90=$(BIN_DIR)/%.o)
MOD_MOD_FILES = $(MOD_FILES:$(MOD_DIR)/%.f90=$(LIB_DIR)/%.mod)
# TESTS
TEST_FILES = $(wildcard $(TEST_DIR)/test*.f90)
EXEC_TEST = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%)

# MAIN PROGRAM
all: $(EXEC)

run: $(EXEC)
	@echo "Processing... This may take a while (Estimated time: 2 minutes)"
	$(EXEC) > $(SRC_DIR)/fortran_results.txt
	$(PYTHON) $(EXEC_PYTHON)

$(EXEC): $(MAIN_FILE) | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(LIB_DIR) -I$(PATH_MOD) -o $@ $< -L$(LIB_DIR) $(LIB) 

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
	$(FC) $(FLAGS) -J$(LIB_DIR) -o $@ -c $(MOD_DIR)/$(LIB_NAME).f90 -L$(LIB_DIR) $(LIB)

$(BIN_DIR)/%.o: $(MOD_DIR)/%.f90 $(LIB_DIR) | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(LIB_DIR) -I$(PATH_MOD) -J$(LIB_DIR)  -c $< -o $@

# UNINSTALLATION
uninstall:
	@echo "Uninstalling library and module files..."
	@for modfile in $(MOD_MOD_FILES); do \
		if [ -f $(PATH_MOD)/$$(basename $$modfile) ]; then \
			echo "Removing module file: $(PATH_MOD)/$$(basename $$modfile)"; \
			sudo rm $(PATH_MOD)/$$(basename $$modfile); \
		fi; \
	done
	@if [ -f $(PATH_LIBRARY)/lib$(LIB_NAME).a ]; then \
		echo "Removing library file: $(PATH_LIBRARY)/lib$(LIB_NAME).a"; \
		sudo rm $(PATH_LIBRARY)/lib$(LIB_NAME).a; \
	fi
	@$(MAKE) --no-print-directory clean

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
	$(FC) $(FLAGS) -I$(LIB_DIR) -I$(PATH_MOD) -o $@ $< -L$(LIB_DIR) $(LIB)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

clean:
	rm -rf $(BIN_DIR) $(LIB_DIR)

.PHONY: all clean test install install_global uninstall run

