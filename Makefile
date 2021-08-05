FLAGS=-std=c++11 -g -Wall -fsanitize=undefined

all: main data

###

main: main.o types.o
	g++ -o main main.o types.o $(FLAGS) -lcairo

main.o: main.cpp
	g++ -c main.cpp $(FLAGS)

types.o: types.cpp
	g++ -c types.cpp $(FLAGS)

###

data: data.cpp
	g++ data.cpp -o data $(FLAGS)

###

.PHONY: clean
clean:
	rm -f *.o main data test.png
