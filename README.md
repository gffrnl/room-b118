# Room-B118 Math Library.

### Authors:
* [Guilherme F. Fornel](https://github.com/gffrnl)             &emsp;&emsp;&emsp;&emsp;                    (<gffrnl@gmail.com>)
* [Fabio Souto de Azevedo](https://github.com/fazedo)          &emsp;&emsp;                                (<fabio.azevedo@ufrgs.br>)

---
### Requirements:
  - C++17 compatible compiler.
  - [CMake](https://cmake.org/) version 3.23 or higher.
  - [FFTW](https://www.fftw.org/) version 3.0 or higher.

### Compilation:

**1. On Linux platforms:**

  ``you@place:/path/to/repo$`` ``mkdir build/ && cd build/``<br>
  ``you@place:/path/to/repo/build$`` ``cmake ..``<br>
  ``you@place:/path/to/repo/build$`` ``make``<br>

  On success, ``you@place:/path/to/repo/build$`` ``ls libgcikpi.so`` returns ``0`` (``echo $?`` prints the last return code).
