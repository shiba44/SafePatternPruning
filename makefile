
CXX ?= g++;
CFLAGS ?= -std=gnu++17 -O3 -Wall -Wextra 
#-Wcast-align -Wcast-qual -Wdelete-non-virtual-dtor -Wdisabled-optimization -Winit-self -Wmissing-include-dirs -Wold-style-cast -Wredundant-decls -Wswitch-enum -Wswitch-default -Wdouble-promotion

REMOVE = rm -f

SRCS := $(wildcard *.cc)
OBJS := $(SRCS:%.cc=%.o)
DEPS := $(SRCS:%.cc=%.d)

default: train

train: train.cc solver.o miner.o datautil.o testutil.o type.o
	$(CXX) $(CFLAGS) -o $@ $^

%.o : %.cc %.h
	$(CXX) $(CFLAGS) -c -MMD -MP $< -o $@

.PHONY: all
all: clean train

.PHONY: clean
clean:
	$(REMOVE) train $(wildcard *.o) $(wildcard *.d) 

-include $(DEPS)
