CXX = g++
CXXFLAGS = -O3 -std=c++11 -w -I. -Isequitur

SRCS = main.cpp sequitur/symbols.cpp sequitur/symbolwrapper.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = seq.exe

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
