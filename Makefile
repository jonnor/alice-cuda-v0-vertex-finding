# Name of your binary
EXECUTABLE	:= alicudav0vertexer
SOLIB := libcudav0vertexer.so
# Device sources (.cu)
CUFILES		:= \
	AliCudaV0vertexer.cu \

# Headers
CUDEPS		:= \
	AliCudaDefs.h \
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
NVCCFLAGS += -Xcompiler -fPIC

# Build dirs
BINDIR := bin
LIBDIR := lib
OBJDIR := obj

# Object files
OBJS := $(patsubst %.cu, $(OBJDIR)/%.o,$(notdir $(CUFILES)))

# Rules
$(LIBDIR)/$(SOLIB): $(OBJS)
	$(NVCC) $(NVCCFLAGS) -shared -o $@ $^

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

