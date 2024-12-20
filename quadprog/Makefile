# Makefile for Mata QuadProg++ plugin (Windows/Mac/Linux compatible)

# Compiler and flags (adjust for your Windows compiler if needed)
CXX = g++
CXXFLAGS = -O2 -fPIC -Wall -std=c++11 -fms-extensions
LDFLAGS = -fms-extensions -shared 

# Stata directory (adjust for each platform)
ifeq ($(OS),Windows_NT)  # Windows
    STATA_DIR = C:/Program Files/Stata17
else
    STATA_DIR = /Applications/Stata
endif

# QuadProg++ directory (adjust for each platform)
QUADPROG_DIR = /Users/davidvandijcke/University\ of\ Michigan\ \
Dropbox/David\ Van\ Dijcke/Flo_GSRA/stata_repo/quadprog

# Platform-specific settings
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)  # Mac
    PLUGIN_FLAGS = -dynamiclib
    SHARED_LIB_EXT = .dylib
endif
ifeq ($(UNAME), Linux)  # Linux
    PLUGIN_FLAGS = -shared -Wl,-soname,
    SHARED_LIB_EXT = .so
endif
ifeq ($(OS),Windows_NT)  # Windows
    PLUGIN_FLAGS = -shared -Wl,-soname,
    SHARED_LIB_EXT = .dll
endif

# Plugin name
PLUGIN_NAME = quadprog_mata$(SHARED_LIB_EXT)

# Object files
OBJECTS = QuadProg++.o quadprog_mata.o

# Default target
all: $(PLUGIN_NAME)

# Linking the plugin
$(PLUGIN_NAME): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(PLUGIN_FLAGS)$@ -o $@ $(OBJECTS)

# Compiling the C++ files
quadprog_mata.o: quadprog_mata.cpp
	$(CXX) $(CXXFLAGS) -I. -I$(STATA_DIR)/include -I$(QUADPROG_DIR) -c $<

QuadProg++.o: QuadProg++.cc
	$(CXX) $(CXXFLAGS) -I. -I$(STATA_DIR)/include -I$(QUADPROG_DIR) -c $<

# Clean target
clean:
	rm -f $(OBJECTS) $(PLUGIN_NAME)

# Install target (adjust for each platform)
install: $(PLUGIN_NAME)
ifeq ($(UNAME), Darwin)  # Mac
	cp quadprog_mata.dylib "/Users/davidvandijcke/Library/Application\ Support/Stata/ado/plus/quadprog_mata.plugin"
endif
ifeq ($(UNAME), Linux)  # Linux
	cp quadprog_mata.so "/home/youruser/stata/ado/plus/quadprog_mata.plugin"  # Adjust path
endif
ifeq ($(OS),Windows_NT)  # Windows
	copy quadprog_mata.dll "C:\Users\youruser\ado\plus\personal\quadprog_mata.plugin"  # Adjust path
endif

# cp quadprog_mata.dylib "/Users/davidvandijcke/Library/Application Support/Stata/ado/plus/quadprog_mata.plugin"