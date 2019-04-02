all:
	g++ -std=c++11 main.cpp fourier.cpp -o fourier
run:
	valgrind ./fourier.out