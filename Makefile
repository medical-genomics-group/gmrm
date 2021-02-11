.PHONY: all clean help info

SOURCEDIR = ./src
BINDIR    = ./bin

SOURCES  := $(wildcard $(SOURCEDIR)/*.cpp)

#SRC_EXCL  =  $(SOURCEDIR)/BayesRRm_mt.cpp
#SRC_EXCL +=  $(SOURCEDIR)/mk_lut.cpp

SOURCES  := $(filter-out $(SRC_EXCL),$(SOURCES))

CXXFLAGS  = -Ofast
CXXFLAGS += -g
#CXXFLAGS += -H
CXXFLAGS += -std=c++17

INCLUDE   = -I$(SOURCEDIR)
INCLUDE  += -I$(BOOST_ROOT)/include

#CXXFLAGS += -DMANVECT

ifeq ($(CXX),g++)

EXEC     ?= ardyh_g
CXX       = mpic++
BUILDDIR  = build_g
CXXFLAGS += -fopenmp
#CXXFLAGS += -march=native -mfma
#CXXFLAGS += -mavx2 -mfma
CXXFLAGS += -march=skylake-avx512 -mfma
CXXFLAGS += -fopt-info-vec-missed=gcc_vec_missed.txt
#CXXFLAGS += -fopt-info-vec=gcc_vec_ok.txt

else ifeq ($(CXX),icpc)

EXEC     ?= ardyh_i
CXX       = mpiicpc
BUILDDIR  = build_i
CXXFLAGS += -qopenmp
CXXFLAGS += -xCORE-AVX512 -qopt-zmm-usage=high
#CXXFLAGS += -xCORE-AVX2, -axCORE-AVX512 -qopt-zmm-usage=high
CXXFLAGS += -qopt-report=2 -qopt-report-phase=vec
else
	@echo "Neither GCC nor Intel compiler available." 1>&2 && false

endif

ifeq (, $(shell which $(CXX)))
$(error "no $(CXX) in $(PATH), please load relevant modules.")
endif

OBJS	:= $(patsubst $(SOURCEDIR)/%.cpp, $(BUILDDIR)/%.o, $(SOURCES))
DEPS	:= $(patsubst $(SOURCEDIR)/%.cpp, $(BUILDDIR)/%.d, $(SOURCES))

LIBS      = -lz

all: create_path $(BINDIR)/$(EXEC)

$(BINDIR)/$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(BUILDDIR)/%.o: $(SOURCEDIR)/%.cpp Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) -MMD -MP -c $< -o $@

-include $(DEPS)

create_path:
	mkdir -p $(BUILDDIR)
	mkdir -p $(BINDIR)

clean:
	rm -vf $(BUILDDIR)/*.o $(BUILDDIR)/*.d $(BINDIR)/$(EXEC)

help:
	@echo "Usage: make [ all | clean | help ]"
