CXX = g++
CXXFLAGS = -Wall -O3
HDRS = ShallowWater.h
LIBS = -lblas
OBJS = main.o ShallowWater.o
TARGET = main

default: $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBS)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)


.PHONY: clean
	
clean: 
	-rm -f .*o myprog
