# Name of your binary
EXECUTABLE	:= test_cuda
# Device sources (.cu)
CUFILES		:= \
	test_cuda.cu \

# Headers
CUDEPS		:= \
	common.h \
	aliexternaltrackparam.cu \
	aliv0vertexer.cu


### COMMON ###
# Paths
CUDASDK := /home/lbratrud/cudaWork
CUDAPATH := /usr/local/cuda
NVCC := $(CUDAPATH)/bin/nvcc

# Flags
INCLUDES  += -I$(CUDAPATH)/include -I$(CUDASDK)/common/inc
COMMONFLAGS += $(INCLUDES)
NVCCFLAGS += $(COMMONFLAGS)

# Build dirs
BINDIR := bin
LIBDIR := lib
OBJDIR := obj

# Object files
OBJS := $(patsubst %.cu, $(OBJDIR)/%.o,$(notdir $(CUFILES)))

# Rules
$(BINDIR)/$(EXECUTABLE): $(OBJS)
	$(NVCC) $(NVCCFLAGS) -o $@ $^

$(OBJDIR)/%.o: %.cu $(CUDEPS) makedirectories
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

# Other things
makedirectories:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)
	@mkdir -p $(LIBDIR)

clean:
	@rm -rf $(OBJDIR)
	@rm -rf $(BINDIR)
	@rm -rf $(LIBDIR)

