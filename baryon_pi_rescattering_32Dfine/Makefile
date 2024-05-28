SHELL=/bin/bash

CXX=/thfs1/software/mpich/mpi-x-gcc10.2.0/bin/mpicxx
qlat_prefix=$(HOME)/Yusheng_Gao/software/Qlat_zyw

QLAT_INCLUDE=$(qlat_prefix)/include
QLAT_LIB=$(qlat_prefix)/lib

# QLAT_INCLUDE_HEADERS := $(shell find $(QLAT_INCLUDE) -name "*.h")$(foreach dir,$(QLAT_INCLUDE_HEADERS),-I$(dir))

QLAT_CFLAGS= -std=c++14 -fopenmp -march=native -O2 -Wall -funroll-loops -Wno-unused-variable -DNO_ALIGNED_ALLO -fno-strict-aliasing -Wfatal-errors
QLAT_CFLAGS+= -I$(QLAT_INCLUDE) -I$(qlat_prefix)/eigen3
#   -I$(QLAT_INCLUDE)/qlat -I$(QLAT_INCLUDE)/qlat-utils

QLAT_CXXFLAGS=$(QLAT_CFLAGS)

QLAT_LDFLAGS= -lgsl -lgslcblas -lm -lfftw3 -lfftw3f -lmpi

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
DEPS=$(SRCS:.cpp=.d)

all: qlat.x

qlat.x: $(OBJS)
	$(CXX) -o $@ $^ $(QLAT_LDFLAGS)

%.o: %.cpp
	$(CXX) $(QLAT_CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f *.o *.d qlat.x

-include $(DEPS)
