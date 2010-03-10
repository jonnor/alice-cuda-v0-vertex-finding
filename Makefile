# Name of your binary
EXECUTABLE	:= test_cuda
# Device sources (.cu)
CUFILES		:= \
	test_cuda.cu \
	aliexternaltrackparam.cu
# Headers
CUDEPS		:= \
	common.h \
	aliexternaltrackparam.h


### COMMON ###
# Directories
CUDASDK := /home/lbratrud/cudaWork
CUDAPATH := /usr/local/cuda

# Compilers
NVCC := $(CUDAPATH)/bin/nvcc

# Includes
INCLUDES  += -I$(CUDAPATH)/include -I$(CUDASDK)/common/inc

# Common flags
COMMONFLAGS += $(INCLUDES)

# Flags
NVCCFLAGS += $(COMMONFLAGS)

# Set up object files
OBJDIR := obj
OBJS := $(patsubst %.cu, $(OBJDIR)/%.o,$(notdir $(CUFILES)))

# Set up binary directory
BINDIR := bin

# Rules
$(BINDIR)/$(EXECUTABLE): $(OBJS)
	$(NVCC) $(NVCCFLAGS) -o $@ $^

$(OBJDIR)/%.o: %.cu $(CUDEPS) makedirectories
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

# Other things
makedirectories:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)
	@echo bin/ created
	@echo obj/ created

clean:
	@rm -rf $(OBJDIR)
	@rm -rf $(BINDIR)
	@echo All clean

