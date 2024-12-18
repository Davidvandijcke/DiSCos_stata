# Makefile for Mata QuadProg++ plugin (Mac version)
CXX = g++
CXXFLAGS = -O2 -fPIC -Wall -std=c++11 -fms-extensions
LDFLAGS = -fms-extensions

# Adjust these paths for your system
STATA_DIR = /Applications/Stata
QUADPROG_DIR = /Users/davidvandijcke/University\ of\ Michigan\ Dropbox/David\ Van\ Dijcke/Flo_GSRA/stata_repo/quadprog

# Mac-specific settings
PLUGIN_FLAGS = -dynamiclib -undefined dynamic_lookup
SHARED_LIB_EXT = .dylib

PLUGIN_NAME = quadprog_mata$(SHARED_LIB_EXT)

# List object files needed
OBJECTS = QuadProg++.o quadprog_mata.o

all: $(PLUGIN_NAME)

$(PLUGIN_NAME): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(PLUGIN_FLAGS) -o $@ $(OBJECTS)

quadprog_mata.o: quadprog_mata.cpp
	$(CXX) $(CXXFLAGS) -I. -I$(STATA_DIR)/include -I$(QUADPROG_DIR) -c quadprog_mata.cpp

QuadProg++.o: QuadProg++.cc
	$(CXX) $(CXXFLAGS) -I. -I$(STATA_DIR)/include -I$(QUADPROG_DIR) -c QuadProg++.cc

clean:
	rm -f $(OBJECTS) $(PLUGIN_NAME)