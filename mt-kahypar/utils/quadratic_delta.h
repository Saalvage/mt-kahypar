#pragma once

#include "../macros.h"

template <int P>
struct quadratic_delta_t;

template <>
struct quadratic_delta_t<0> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return 0;
	}
};

template <>
struct quadratic_delta_t<1> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return delta;
	}
};

template <>
struct quadratic_delta_t<2> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return (2 * a + delta) * delta;
	}
};

template <>
struct quadratic_delta_t<3> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return a * (3 * a * delta + 3 * delta * delta) + delta * delta * delta;
	}
};

template <>
struct quadratic_delta_t<4> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return a * (a * (4 * a * delta + 6 * delta * delta) + 4 * delta * delta * delta) + delta * delta * delta * delta;
	}
};

template <>
struct quadratic_delta_t<5> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return a * (a * (a * (5 * a * delta + 10 * delta * delta) + 10 * delta * delta * delta) + 5 * delta * delta * delta * delta) + delta * delta * delta * delta * delta;
	}
};

// Calculates (a+delta)^2 - a^2 while avoiding cancellation.
template <typename F>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
static F quadratic_delta(F a, F delta) {
	return quadratic_delta_t<2>::compute(a, delta);
}
