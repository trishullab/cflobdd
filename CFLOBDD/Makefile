# Project Name (executable)
PROJECT = cflobdd
# Compiler
CC = g++-11

# Run Options       
COMMANDLINE_OPTIONS = #/dev/ttyS0

# Compiler options during compilation
COMPILE_OPTIONS = -g -std=c++2a -w -shared -Wall -Wextra -DHAVE_CONFIG_H -Werror -Wunused-but-set-variable 
# -ansi -pedantic -Wall 

#Header include directories
HEADERS = -I. -I../../../../boost_1_81_0/ -I. -I.Solver/uwr/bit_vector/ -I.Solver/uwr/assert/ -I.Solver/uwr/matrix/ -I.Solver/uwr/parsing/
#Libraries for linking
LIBS = 

# Dependency options
DEPENDENCY_OPTIONS = -MM


# Subdirs to search for additional source files
SOURCE_FILES := $(shell ls *.cpp)
SOURCE_FILES += $(shell ls Solver/uwr/bit_vector/*.cpp)
SOURCE_FILES += $(shell ls Solver/uwr/parsing/*.cpp)

# Create an object file of every cpp file
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCE_FILES))

# Dependencies
DEPENDENCIES = $(patsubst %.cpp, %.d, $(SOURCE_FILES))

dep:
	cd cflobdd/cudd-complex-big && make all && cd ../../

# Make $(PROJECT) the default target
all: $(PROJECT)
#$(DEPENDENCIES) -shared 
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
