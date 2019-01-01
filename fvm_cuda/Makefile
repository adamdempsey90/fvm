EXECUTABLE=fvmcuda
CUSOURCES=prob.cu algogas.cu reconstruct.cu riemann.cu update.cu exact.cu hllc.cu anrs.cu boundary.cu calc_dt.cu source.cu conduction.cu viscosity.cu config.cu driver.cu
SOURCES=allocate.c output.c read_pars.c fvm.c 
HEADER=prob.h defs.h structs.h prototypes.h cuda_defs.h 



$(info SYSTYPE: "$(SYSTYPE)")

BIN=bin/
CUBIN=bin/cuda/
SRC=src/
OUT=out/

$(shell mkdir -p $(BIN))
$(shell mkdir -p $(CUBIN))
$(shell mkdir -p $(OUT))

NVCC=nvcc

ifeq ("$(SYSTYPE)","quest")
CFLAGS=-O3 --compiler-options -Wall
NVCCFLAGS=-arch=sm_37
INCLIB=-I/software/hdf5/1.8.12-serial/include/
LDLIB=-L/software/hdf5/1.8.12-serial/lib/
endif
ifeq ("$(SYSTYPE)","Ubuntu")
CFLAGS=-O3 --compiler-options -Wall
NVCCFLAGS=-arch=sm_60
INCLIB=-I/usr/lib/x86_64-linux-gnu/hdf5/serial/include/
LDLIB=-L/usr/lib/x86_64-linux-gnu/hdf5/serial/
endif

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
	$(NVCC)  $(NVCCFLAGS) $(LDLIB) $(OBJ) -lm -lhdf5 -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER)
	$(NVCC) $(CFLAGS) -x cu $(NVCCFLAGS) $(INCLIB) -dc $< -o $@
$(CUBIN)%.o: $(SRC)%.cu $(CHEADER)
	$(NVCC) $(CFLAGS) -x cu $(NVCCFLAGS) $(INCLIB) -dc $< -o $@
clean:
	rm $(OBJ) $(EXECUTABLE)

