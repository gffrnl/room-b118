// Copyright 2024 Guilherme F. Fornel

#pragma once

#include <cstddef>

void fast_symm_toeplitz_prod(std::size_t const                      n,
                             double      const * const __restrict__ A1,
                             double      const * const __restrict__ x,
                             double            * const __restrict__ b);
