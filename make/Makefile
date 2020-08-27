FC = gfortran
FFLAGS = -Wall

EXE_DIR = ..

UPDATE = update

SRC_FILES := \
	$(wildcard ../src/*.f90) 

OBJ_FILES := $(notdir $(SRC_FILES:.f90=.o))
SRC_DIR := $(dir $(SRC_FILES))
VPATH := $(SRC_DIR)

$(UPDATE): $(OBJ_FILES)
	$(FC) $(FFLAGS) -o $(EXE_DIR)/a.out $(OBJ_FILES) && \
	touch $(UPDATE)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf $(OBJ_FILES) $(EXE_DIR)/a.out *.mod $(UPDATE)

-include deps
