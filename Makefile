all:	
	g++ main.cpp -o main -g -lm -lcairo -Wall -std=c++11 -fsanitize=undefined

.PHONY: clean
clean:
	rm main test.png
