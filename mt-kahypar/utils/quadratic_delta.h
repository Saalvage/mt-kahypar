#pragma once

#include "../macros.h"

// Calculates (a+delta)^2 - a^2 while avoiding cancellation.
template <typename F>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
static F quadratic_delta(F a, F delta) {
    return 2 * a * delta + delta * delta;
}
