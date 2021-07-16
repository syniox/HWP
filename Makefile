FLAGS=-g -lm -lcairo -Wall -std=c++11 -fsanitize=undefined
CC=g++

main: main.cpp structs.o
	$(CC) main.cpp structs.o -o main $(FLAGS)
structs.o: structs.cpp
	$(CC) -c structs.cpp $(FLAGS)
draw.o: draw.cpp 

.PHONY: clean
clean:
	rm -f main *.o test.png
