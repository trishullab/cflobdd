#ifndef BIT_VECTOR_OPS_H
#define BIT_VECTOR_OPS_H

// Arithmetic operations on bit-vectors that use their abstract data type,
// rather than being part of their abstract data type.
#include "bit_vector_1.h"
#include "../assert/uw_assert.h"
// This isn't actually needed, yet, but might be useful.
// If you want what's in the BV_Type_Map header, unwind code so that 
// it winds up here.
// #include "tsl/analysis_components/src/reinterps/ara/ks/BV_Type_Map.hpp"

namespace BitVector
{
	template <typename T> T mult_pow2(T a, size_t w); // compute (a << w)
	template <typename T> T div_pow2(T a, size_t w); // compute (a >> w)
	template <typename T> T multiplicative_inverse(T a);
	template <typename T> size_t compute_rank(T a);
	template <typename T> size_t highest_power();
	
	template <typename T>
	inline size_t highest_power()
	{
		UWAssert::UWAssert::shouldNeverHappen();
	}

	template <>
	inline size_t highest_power<BitVector1>()
	{
		return 1;
	}

	// Compute a * 2^w
	template <typename T>
	inline T mult_pow2(T a, size_t w) {
		return (w == highest_power<T>()) ? T(0) : (a << w);
	}

	// Compute a / 2^w
	// Note: while using this to find the invertible part of a number, 
	// --> 0 = 2^32  * 1 <--
	// and this macro should not be used
	template <typename T>
	inline T div_pow2(T a, size_t w) {
		return (w == highest_power<T>()) ? T(0) : (a >> w);
	}
	
	//typedef unsigned int rank_t;
	
	//
	// compute_rank
	//
	// Version of Fig. 5-14 from Hacker's Delight that uses a loop
	// so that it can be parameterized on highest_power
	/*template <typename T>
	inline size_t compute_rank(BitVectorTemplate<T> x) {
		typedef BitVectorTemplate<T> BV;
		BV y;
		size_t n;

		if (x == BV(0)) return highest_power<BV>();

		n = highest_power<BV>() - 1;
		for (size_t shift = highest_power<BV>() / 2; shift != 0; shift /= 2) {
			y = x << shift;
			if (!(y == BV(0))) {
				n = n - shift;
				x = y;
			}
		}
		return n;
	}*/

	template <typename T>
	inline size_t compute_rank(T) {
		UWAssert::UWAssert::shouldNeverHappen();
	}
	
	template <>
	inline size_t compute_rank<BV1>(BV1 x) {
		if (x == (BV1)0) return 1u;
		else  return 0u;
	}
	
	// find the multiplicative inverse
	template <typename T>
	T multiplicative_inverse(T d) {
		T xn, t;
		if (d % 2 != 1)
			UWAssert::UWAssert::shouldNeverHappen();

		xn = d;
		while (true) {
			t = d*xn;
			if (t == 1) return xn;
			xn = xn*((T)2 - t);
		}
	}

	// A template specialization of multiplicative_inverse for BV1
	//
	template <>
	inline BV1 multiplicative_inverse(BV1 d) {
		if (d != 1)
			UWAssert::UWAssert::shouldNeverHappen();
		return d;
	}

	// Return the min pop count of val and -val. This is useful to figure out how complex the mutiplication with val would be
	template <typename T>
	size_t GetMinPosOrNegPopCount(T val)
	{
		T minus_val = -val;
		return (std::min) (val.pop_count(), minus_val.pop_count());
	}
}

#endif
