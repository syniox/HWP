CXXFLAGS+=-std=c++11 -g -Wall -Wextra -fsanitize=undefined

all: main data

###

main: main.o types.o draw.o utils.o
	g++ $(CXXFLAGS) -o main $^ -lcairo

main.o: types.h draw.h utils.h
types.o: types.h
draw.o: draw.h types.h utils.h

###

data: data.cpp

###

.PHONY: clean
clean:
	rm -f *.o main data *.png
