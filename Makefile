SRC_DIR = src
TEST_DIR = test
BIN_DIR = bin

SRC_FILES = $(SRC_DIR)/dgemm_wrapper.f90
TEST_FILES = $(TEST_DIR)/test_dgemm.f90

OBJ_FILES = $(BIN_DIR)/dgemm_wrapper.o $(BIN_DIR)/test_dgemm.o
EXE = $(BIN_DIR)/test_dgemm

FC = gfortran
FCFLAGS = -O2 -g
LIBS = -lblas -llapack

all: $(EXE)

# Create bin directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Compile the dgemm wrapper to an object file
$(BIN_DIR)/dgemm_wrapper.o: $(SRC_FILES) | $(BIN_DIR)
	$(FC) $(FCFLAGS) -c $(SRC_FILES) -o $(BIN_DIR)/dgemm_wrapper.o

# Compile the test program to an object file
$(BIN_DIR)/test_dgemm.o: $(TEST_FILES) | $(BIN_DIR)
	$(FC) $(FCFLAGS) -c $(TEST_FILES) -o $(BIN_DIR)/test_dgemm.o

# Link the test program with the dgemm wrapper
$(EXE): $(OBJ_FILES)
	$(FC) $(FCFLAGS) $(OBJ_FILES) -o $(EXE) $(LIBS)
	rm -f $(OBJ_FILES)

# Clean up object files and executables
clean:
	rm -f $(OBJ_FILES) $(EXE)

distclean: clean
	rm -rf $(BIN_DIR)

.PHONY: all clean distclean

