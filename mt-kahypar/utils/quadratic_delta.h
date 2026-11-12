#pragma once

#include "../macros.h"

#include <utility>

template<typename T, std::size_t... pack>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
T pow_impl(T value, std::index_sequence<pack...>) {
	if constexpr (sizeof...(pack) == 0) {
		return T{1};
	} else {
		return ((pack, value) * ...);
	}
}


template<std::size_t exponent, typename T>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
T pow(T value) {
	return pow_impl( value, std::make_index_sequence<exponent>());
}

template <int P>
struct quadratic_delta_t;

template <>
struct quadratic_delta_t<0> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F, F) {
		return 0;
	}
};

template <>
struct quadratic_delta_t<1> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F, F delta) {
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
		return a * (3 * a * delta + 3 * pow<2>(delta)) + pow<3>(delta);
	}
};

template <>
struct quadratic_delta_t<4> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return a * (a * (4 * a * delta + 6 * pow<2>(delta)) + 4 * pow<3>(delta)) + pow<4>(delta);
	}
};

template <>
struct quadratic_delta_t<5> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return a * (a * (a * (5 * a * delta + 10 * pow<2>(delta)) + 10 * pow<3>(delta)) + 5 * pow<4>(delta)) + pow<5>(delta);
	}
};

template <>
struct quadratic_delta_t<6> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return a * (a * (a * (a * (6 * a * delta + 15 * pow<2>(delta)) + 20 * pow<3>(delta)) + 15 * pow<4>(delta)) + 6 * pow<5>(delta)) + pow<6>(delta);
	}
};

template <>
struct quadratic_delta_t<7> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return a * (a * (a * (a * (a * (7 * a * delta + 21 * pow<2>(delta)) + 35 * pow<3>(delta)) + 35 * pow<4>(delta)) + 21 * pow<5>(delta)) + 7 * pow<6>(delta)) + pow<7>(delta);
	}
};

template <>
struct quadratic_delta_t<8> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return a * (a * (a * (a * (a * (a * (8 * a * delta + 28 * pow<2>(delta)) + 56 * pow<3>(delta)) + 70 * pow<4>(delta)) + 56 * pow<5>(delta)) + 28 * pow<6>(delta)) + 8 * pow<7>(delta)) + pow<8>(delta);
	}
};

template <>
struct quadratic_delta_t<9> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return delta * (9 * pow<8>(a) + 36 * pow<7>(a) * delta + 84 * pow<6>(a) * pow<2>(delta) + 126 * pow<5>(a) * pow<3>(delta) + 126 * pow<4>(a) * pow<4>(delta) + 84 * pow<3>(a) * pow<5>(delta) + 36 * pow<2>(a) * pow<6>(delta) + 9 * a * pow<7>(delta) + pow<8>(delta));
	}
};

template <>
struct quadratic_delta_t<10> {
	template <typename F>
	MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
	static F compute(F a, F delta) {
		return delta * (10 * pow<9>(a) + 45 * pow<8>(a) * delta + 120 * pow<7>(a) * pow<2>(delta) + 210 * pow<6>(a) * pow<3>(delta) + 252 * pow<5>(a) * pow<4>(delta) + 210 * pow<4>(a) * pow<5>(delta) + 120 * pow<3>(a) * pow<6>(delta) + 45 * pow<2>(a) * pow<7>(delta) + 10 * a * pow<8>(delta) + pow<9>(delta));
	}
};

// Calculates (a+delta)^2 - a^2 while avoiding cancellation.
template <typename F>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
static F quadratic_delta(F a, F delta) {
	return quadratic_delta_t<2>::compute(a, delta);
}
