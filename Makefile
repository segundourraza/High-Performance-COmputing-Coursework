CXX = g++
CXXFLAGS = -Wall -O0 -g
HDRS = ShallowWater.h
LIBS = -lblas -lboost_program_options
OBJS = main.o ShallowWater.o
TARGET = main

default: $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBS)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)

test1: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1

test2: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2

test3: $(TARGET) 
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3

test4: $(TARGET)
	./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4

validation1: $(TARGET)
	./$(TARGET) --dt 0.1 --T 20 --Nx 100 --Ny 100 --ic 1

validation2: $(TARGET)
	./$(TARGET) --dt 0.1 --T 20 --Nx 100 --Ny 100 --ic 2

validation3: $(TARGET)
	./$(TARGET) --dt 0.1 --T 20 --Nx 100 --Ny 100 --ic 3

validation4: $(TARGET)
	./$(TARGET) --dt 0.1 --T 20 --Nx 100 --Ny 100 --ic 4
.PHONY: clean
	
clean: 
	-rm -f *.o $(TARGET)
