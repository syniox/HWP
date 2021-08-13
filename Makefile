CXXFLAGS+=-std=c++11 -g -Wall -Wextra -fsanitize=undefined

all: main data

###

main: main.o types.o draw.o
	g++ $(CXXFLAGS) -o main main.o types.o draw.o -lcairo

###

data: data.cpp

###

.PHONY: clean
clean:
	rm -f *.o main data test.png
