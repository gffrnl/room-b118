rm -f ./bench-701
rm -f bench-701.dat
g++ -g -std=c++17 bench-701-params.cpp main.cpp -o bench-701 -lb118-linalg -lb118-frlap -lopenblas -lfftw3 -lgsl -O3 -march=native -DNDEBUG
