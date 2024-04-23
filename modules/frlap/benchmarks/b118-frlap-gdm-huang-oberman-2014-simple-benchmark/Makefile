all: bench-701

#calc.a: main.cpp ConvergenceBenchmark.hpp
#	g++ -std=c++17 -Wall -Wpedantic -fPIC -c main.cpp
#	ar rvs calc.a main.o

bench-701.o: main.cpp simple_calculator.hpp
	g++ -std=c++17 -Wall -Wpedantic -c main.cpp -o bench-701.o

bench-701: bench-701.o
	g++ -std=c++17 -Wall -Wpedantic bench-701.o -o bench-701 -lfftw3f -lfftw3 -lfftw3l -lfftw3q -lopenblas -lgsl

run: bench-701
	./bench-701

clean:
	rm -f *.o
	rm -f *.a
	rm -f bench-701

purge: clean
	rm -f *~

.PHONY: clean purge run
