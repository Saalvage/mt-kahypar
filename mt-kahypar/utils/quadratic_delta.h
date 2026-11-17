#pragma once

#include "../macros.h"

#include <utility>

template <typename T, std::size_t... pack>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
T pow_impl(T value, std::index_sequence<pack...>) {
	if constexpr (sizeof...(pack) == 0) {
		return T{1};
	} else {
		return ((pack, value) * ...);
	}
}


template <std::size_t exponent, typename T>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
T pow(T value) {
	return pow_impl(value, std::make_index_sequence<exponent>());
}

template <int P>
struct quadratic_delta_t {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		if constexpr (P == 0) {
			return 0;
		} else if constexpr (P == 1) {
			return delta;
		} else {
			return quadratic_delta_t<P - 1>::compute(a, delta) * a + delta * pow<P - 1>(a + delta);
		}
	}
};

// Calculates (a+delta)^2 - a^2 while avoiding cancellation.
template <typename F>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
static F quadratic_delta(F a, F delta) {
	return quadratic_delta_t<MT_KAHYPAR_LX_NORM>::compute(a, delta);
}
