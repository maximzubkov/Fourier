all:
	g++ -std=c++11 main.cpp  -o fourier
run:
	valgrind ./fourier.out