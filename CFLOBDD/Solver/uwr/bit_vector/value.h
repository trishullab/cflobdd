#ifndef BIT_VECTOR_VALUE_H
#define BIT_VECTOR_VALUE_H

#include "bit_vector_1.h"
#include "../assert/uw_assert.h"
#include <sys/types.h>

namespace  CFL_OBDD
{
	namespace BitVector
	{
		typedef unsigned long LONG_UINT;
		typedef long LONG_INT;
		enum Bitsize { sixty_four = 64, thirty_two = 32, sixteen = 16, eight = 8, one = 1 };

		Bitsize max_bitsize(Bitsize b1, Bitsize b2);
		
		LONG_UINT GetBitWidth(Bitsize b);

		Bitsize GetBitsize(LONG_UINT w);

		class Value {
		public:
			// Default empty constructor (TODO: Delete as this will create bitsize issues in future)
			Value(const Bitsize size = thirty_two);
			//Value(const LONG_UINT val);
			Value(const Bitsize size, const LONG_UINT val);
			/*Value(const BV64 val);
			Value(const BV32 val);
			Value(const BV16 val);
			Value(const BV8 val);*/
			Value(const BV1 val);
			~Value();

			const Value& operator+() const;
			const Value operator-() const;
			const Value operator~() const;

			bool is_zero() const;
			const Value operator+(const Value& v) const;
			const Value operator-(const Value& v) const;
			const Value operator*(const Value& v) const;
			const Value operator%(const Value& v) const;

			// Signed modulo operation which interprets this and v as signed
			const Value smod(const Value& v) const;
			const Value operator/(const Value& v) const;

			// Signed division operation which interprets this and v as signed
			const Value sdiv(const Value& v) const;
			const Value operator&(const Value& v) const;
			const Value operator|(const Value& v) const;
			const Value operator^(const Value& v) const;
			const Value operator>>(const Value& v) const;
			const Value operator<<(const Value& v) const;
			bool operator< (const Value& v) const;
			bool operator> (const Value& v) const;
			bool operator<= (const Value& v) const;
			bool operator>= (const Value& v) const;
			bool operator== (const Value& v) const;
			bool operator!= (const Value& v) const;
			Value& operator=(const Value& v);
			Value ZeroExtend(Bitsize new_bitsize) const;
			Value SignExtend(Bitsize new_bitsize) const;
			Value Trunc(Bitsize new_bitsize) const;
			Value SignExtendOrTrunc(Bitsize new_bitsize) const;

			Bitsize GetBitsize() const;

			Bitsize size_;
			/*BV64 bv64_;
			BV32 bv32_;
			BV16 bv16_;
			BV8 bv8_;*/
			BV1 bv1_;
		};

		size_t compute_value_rank(Value v);
		Value multiplicative_inverse_value(Value a);

		size_t GetMinPosOrNegPopCount(Value v);

		std::ostream& operator<<(std::ostream& o, const Value v);
	}


} 

#endif
