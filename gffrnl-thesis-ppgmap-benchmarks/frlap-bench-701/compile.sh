rm -f ./bench-701
g++ -g -std=c++17 bench-701-params.cpp main.cpp -o bench-701 -lb118-frlap -lopenblas -lfftw3 -O3 -march=native -DNDEBUG
