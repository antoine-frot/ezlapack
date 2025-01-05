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
# MAIN PROGRAM
EXEC = $(BIN_DIR)/main_program
MAIN_FILE = $(SRC_DIR)/main_program.f90
#INSTALLATION
MOD_FILES = $(wildcard $(MOD_DIR)/module_*.f90)
MOD_OBJ_FILES = $(MOD_FILES:$(MOD_DIR)/%.f90=$(BIN_DIR)/%.o)
#TESTS
EXEC_TEST = $(TEST_FILES:$(TEST_DIR)/%.f90=$(BIN_DIR)/%)
TEST_FILES = $(wildcard $(TEST_DIR)/test*.f90)

# MAIN PROGRAM
all: $(EXEC)
	$(EXEC)

$(EXEC): $(MAIN_FILE) | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(PATH_MOD) -o $@ $< $(LIB) 

# INSTALLATION
install: $(BIN_DIR)/lib$(LIB_NAME).a
	@echo "Installing library..."
	sudo mkdir -p $(PATH_LIBRARY) $(PATH_MOD)
	sudo cp $(BIN_DIR)/lib$(LIB_NAME).a $(PATH_LIBRARY)
	sudo cp $(BIN_DIR)/*.mod $(PATH_MOD)
	sudo ldconfig

$(BIN_DIR)/lib$(LIB_NAME).a: $(BIN_DIR)/$(LIB_NAME).o $(MOD_OBJ_FILES)
	@echo "Archiving library..."
	ar rcs $@ $^

$(BIN_DIR)/$(LIB_NAME).o: $(MOD_OBJ_FILES)
	$(FC) $(FLAGS) -J$(BIN_DIR) -o $@ -c $(MOD_DIR)/$(LIB_NAME).f90 $(LIB)

$(BIN_DIR)/%.o: $(MOD_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(BIN_DIR) -J$(BIN_DIR)  -c $< -o $@

# TESTS
test: $(EXEC_TEST)
	@echo "Running tests..."
	@for exec in $(EXEC_TEST:$(BIN_DIR)/%=%); do \
		echo "Running $$exec..."; \
		$(BIN_DIR)/$$exec; \
		sleep 1; \
	done

$(BIN_DIR)/%: $(TEST_DIR)/%.f90 | $(BIN_DIR)
	$(FC) $(FLAGS) -I$(PATH_MOD) -o $@ $< $(LIB)


# Create binary directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	sudo rm -rf $(BIN_DIR)

.PHONY: all clean test install

