# Project Name (executable)
# PROJECT = cflobdd
# Detect the operating system
ifeq ($(OS),Windows_NT)
    # Windows (MinGW)
    PROJECT      = cflobdd.exe
else
    # Linux/Unix-like
    PROJECT      = cflobdd
endif

# Compiler
# CC = g++-11
CC = g++

# Run Options       
COMMANDLINE_OPTIONS = #/dev/ttyS0

# Compiler options during compilation
COMPILE_OPTIONS = -g -O3 -std=c++2a -w -DCFLOBDD_C_EXPORTS
# -ansi -pedantic -Wall 

#Header include directories
HEADERS = -I. -I../../boost_1_81_0/ -I.Solver/uwr/bit_vector/ -I.Solver/uwr/assert/ -I.Solver/uwr/matrix/ -I.Solver/uwr/parsing/
#Libraries for linking
LIBS =

# Dependency options
DEPENDENCY_OPTIONS = -MM

#-- Do not edit below this line --

# Subdirs to search for additional source files
SOURCE_FILES += $(shell ls *.cpp)
SOURCE_FILES += $(shell ls Solver/uwr/bit_vector/*.cpp)
SOURCE_FILES += $(shell ls Solver/uwr/parsing/*.cpp)
# SOURCE_FILES += $(shell find . -maxdepth 1 -mindepth 1 -name \*.cpp -a -not -name main.cpp)

# Create an object file of every cpp file
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCE_FILES))

# Dependencies
DEPENDENCIES = $(patsubst %.cpp, %.d, $(SOURCE_FILES))

# Create .d files
%.d: %.cpp
	$(CC) $(DEPENDENCY_OPTIONS) $< -MT "$*.o $*.d" -MF $*.d

# Make $(PROJECT) the default target
all: $(PROJECT)
#$(DEPENDENCIES) -shared -dynamiclib 
$(PROJECT): $(OBJECTS)
	$(CC) -o $(PROJECT) $(OBJECTS) $(LIBS)

# Include dependencies (if there are any)
# ifneq "$(strip $(DEPENDENCIES))" ""
#   include $(DEPENDENCIES)
# endif

# Compile every cpp file to an object
# %.cpp 
%.o: %.cpp
	$(CC) -c $(COMPILE_OPTIONS) -o $@ $^ $(HEADERS)

# Build & Run Project
run: $(PROJECT)
	./$(PROJECT) $(COMMANDLINE_OPTIONS)

# Clean & Debug
.PHONY: makefile-debug
makefile-debug:

.PHONY: clean
clean:
	rm -f $(PROJECT) $(OBJECTS)

.PHONY: depclean
depclean:
	rm -f $(DEPENDENCIES)

clean-all: clean depclean

