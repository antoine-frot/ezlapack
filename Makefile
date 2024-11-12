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

$(EXE): $(OBJ_FILES)
	$(FC) $(OBJ_FILES) -o $(EXE) $(LIBS)

$(BIN_DIR)/dgemm_wrapper.o: $(SRC_DIR)/dgemm_wrapper.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/dgemm_wrapper.f90 -o $(BIN_DIR)/dgemm_wrapper.o

$(BIN_DIR)/test_dgemm.o: $(TEST_DIR)/test_dgemm.f90
	$(FC) $(FCFLAGS) -c $(TEST_DIR)/test_dgemm.f90 -o $(BIN_DIR)/test_dgemm.o

clean:
	rm -f $(OBJ_FILES) $(EXE)

distclean: clean
	rm -f $(BIN_DIR)/$(EXE)

.PHONY: all clean distclean

