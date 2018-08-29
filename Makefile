EXECUTABLE=run
SOURCES=metric.c prob.c algogas.c output.c fvm.c
HEADER=defs.h structs.h prototypes.h



CFLAGS=-c -fopenmp -Wall -O3 -g

HDF5LIB=-lhdf5
MATHLIB=-lm

INCLIB=
LDLIB=


BIN=bin/
SRC=src/
CC=gcc

UNAME := $(shell echo $(USER))


LDFLAGS= $(MATHLIB) $(HDF5LIB)


#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))




all: $(CSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS)
	$(CC)  $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER)
	$(CC) $(INCLIB) $(CFLAGS) $< -o $@


clean:
	rm $(COBJECTS) $(EXECUTABLE)
