FLAGS=-std=c++11 -lm -g -Wall -fsanitize=undefined

all: main data

main: main.cpp
	g++ main.cpp -o main $(FLAGS) -lcairo

data:data.cpp
	g++ data.cpp -o data $(FLAGS)

.PHONY: clean
clean:
	rm -f main data test.png
