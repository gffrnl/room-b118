#pragma once

#include <b118/vendor/cppreference/equal_within_ulps.hpp>

namespace b118 {
    template<typename Real>
    using equal_within_ulps = cppreference::equal_within_ulps;
}
