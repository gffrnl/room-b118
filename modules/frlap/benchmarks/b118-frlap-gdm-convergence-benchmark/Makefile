all: Example

calc.a: main.cpp ConvergenceBenchmark.hpp
	g++ -g -std=c++17 -Wall -Wpedantic -fPIC -c main.cpp
	ar rvs calc.a main.o

Example.o: Example.cpp ConvergenceBenchmark.hpp
	g++ -g -std=c++17 -Wall -Wpedantic -c Example.cpp

Example: Example.o calc.a
	g++ -g -std=c++17 -Wall -Wpedantic Example.o calc.a -o Example -lfftw3f -lfftw3 -lfftw3l -lfftw3q -lgsl -lopenblas

run: Example
	./Example

clean:
	rm -f *.o
	rm -f *.a
	rm -f Example

purge: clean
	rm -f *~

.PHONY: clean purge run
