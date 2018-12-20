EXECUTABLE=fvmcuda
CUSOURCES=prob.cu algogas.cu reconstruct.cu riemann.cu update.cu exact.cu hllc.cu anrs.cu boundary.cu calc_dt.cu source.cu conduction.cu viscosity.cu config.cu
SOURCES=allocate.c output.c read_pars.c fvm.c 
HEADER=prob.h defs.h structs.h prototypes.h cuda_defs.h 




BIN=bin/
CUBIN=bin/cuda/
SRC=src/
OUT=out/

$(shell mkdir -p $(BIN))
$(shell mkdir -p $(CUBIN))
$(shell mkdir -p $(OUT))

NVCC=nvcc

CFLAGS=-O3 --compiler-options -Wall
INCLIB=-I/usr/lib/x86_64-linux-gnu/hdf5/serial/include/
LDLIB=-L/usr/lib/x86_64-linux-gnu/hdf5/serial/

OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))

CUOBJECTS=$(CUSOURCES:.cu=.o)
CUDASOURCES=$(addprefix $(SRC),$(CUSOURCES))
CUDAOBJECTS=$(addprefix $(CUBIN),$(CUOBJECTS))


OBJ=$(COBJECTS) $(CUDAOBJECTS)

all: $(CSOURCES) $(CUDASOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS) $(CUDAOBJECTS)
	$(NVCC) -arch=sm_60 $(LDLIB) $(OBJ) -lm -lhdf5 -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER)
	$(NVCC) $(CFLAGS) -x cu -arch=sm_60 $(INCLIB) -dc $< -o $@
$(CUBIN)%.o: $(SRC)%.cu $(CHEADER)
	$(NVCC) $(CFLAGS) -x cu -arch=sm_60 $(INCLIB) -dc $< -o $@
clean:
	rm $(OBJ) $(EXECUTABLE)

