main: main.cpp
	g++ main.cpp -o main -g -lm -lcairo -Wall -std=c++11 -fsanitize=undefined

.PHONY: clean
clean:
	rm -f main test.png
