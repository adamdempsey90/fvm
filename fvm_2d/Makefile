EXECUTABLE=fvm2d
SOURCES=prob.c allocate.c algogas.c output.c fvm.c reconstruct.c riemann.c update.c read_pars.c exact.c hllc.c anrs.c plm.c boundary.c
HEADER=defs.h structs.h prototypes.h



CFLAGS=-c -Wall -O3 -g

HDF5LIB=-lhdf5 -lgomp
MATHLIB=-lm

INCLIB=-I/usr/lib/x86_64-linux-gnu/hdf5/serial/include/
LDLIB=-L/usr/lib/x86_64-linux-gnu/hdf5/serial/


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
