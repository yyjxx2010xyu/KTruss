OBJ = main.o graph.o util.o edge.o
TARGET = main
INCL = -Iinclude
CPP = g++

$(TARGET): $(OBJ)
	$(CPP) $(OBJ) -fopenmp $(INCL) -o $(TARGET) -mpopcnt

%.o: %.cpp
	$(CPP) -O3 -std=c++11 -fopenmp $(INCL) -c $< -o $@ -mpopcnt

.PHONY: clean  
clean:
	rm -rf $(OBJ) $(TARGET) 