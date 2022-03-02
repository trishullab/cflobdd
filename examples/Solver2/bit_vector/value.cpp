#include <algorithm>
#include "value.h"
#include "bit_vector_ops.h"


namespace BitVector
{
	Bitsize max_bitsize(Bitsize b1, Bitsize b2) {
		LONG_UINT max_width = std::max(GetBitWidth(b1), GetBitWidth(b2));
		switch (max_width)
		{
			case 64:
				return sixty_four;
			case 32:
				return thirty_two;
			case 16:
				return sixteen;
			case 8:
				return eight;
			case 1:
				return one;
			default:
				UWAssert::UWAssert::shouldNeverHappen();
				return one;
		}
		return one;
	}

	LONG_UINT GetBitWidth(Bitsize b) {
		switch (b)
		{
			case sixty_four:
				return 64ull;
			case thirty_two:
				return 32ull;
			case sixteen:
				return 16ull;
			case eight:
				return 8ull;
			case one:
				return 1ull;
		}

		UWAssert::UWAssert::shouldNeverHappen();
		return 0ull;
	}

	Bitsize GetBitsize(LONG_UINT w) {
		switch (w)
		{
			case 64:
				return sixty_four;
			case 32:
				return thirty_two;
			case 16:
				return sixteen;
			case 8:
				return eight;
			case 1:
				return one;
			default:
				UWAssert::UWAssert::shouldNeverHappen();
			return one;
		}
		return one;
	}

	// Default empty constructor
	Value::Value(const Bitsize size)
		: size_(size)
	{
	}

	/*Value::Value(const LONG_UINT val)
		: size_(thirty_two)
	{
		bv32_ = (UINT)val;
	}*/

	Value::Value(const Bitsize size, const LONG_UINT val)
		: size_(size)
	{
		switch (size_)
		{
			/*case sixty_four:
				bv64_ = val;
				break;
			case thirty_two:
				bv32_ = (UINT)val;
				break;
			case sixteen:
				bv16_ = (USHORT)val;
				break;
			case eight:
				bv8_ = (UCHAR)val;
				break;*/
			case one:
				bv1_ = (bool) val % 2;
		}
	}


	/*Value::Value(const BV64 val)
		: size_(sixty_four), bv64_(val)
	{}

	Value::Value(const BV32 val)
		: size_(thirty_two), bv32_(val)
	{}

	Value::Value(const BV16 val)
		: size_(sixteen), bv16_(val)
	{}

	Value::Value(const BV8 val)
		: size_(eight), bv8_(val)
	{}
	*/
	Value::Value(const BV1 val)
		: size_(one), bv1_(val)
	{}

	Value::~Value() {}

	const Value& Value::operator+() const { return *this; }

	const Value Value::operator-() const
	{
		switch (size_)
		{
			/*case sixty_four:
				return Value(-bv64_);
			case thirty_two:
				return Value(-bv32_);
			case sixteen:
				return Value(-bv16_);
			case eight:
				return Value(-bv8_);*/
			case one:
				return Value(-bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return *this;
	}

	const Value Value::operator~() const
	{
		switch (size_)
		{
			/*case sixty_four:
				return Value(~bv64_);
			case thirty_two:
				return Value(~bv32_);
			case sixteen:
				return Value(~bv16_);
			case eight:
				return Value(~bv8_);*/
			case one:
				return Value(~bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return *this;
	}

	bool Value::is_zero() const
	{
		switch (size_)
		{
			/*case sixty_four:
				return (bv64_ == BV64(0));
			case thirty_two:
				return (bv32_ == BV32(0));
			case sixteen:
				return (bv16_ == BV16(0));
			case eight:
				return (bv8_ == BV8(0));*/
			case one:
				return (bv1_ == BV1(0));
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return false;
	}

	const Value Value::operator+(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
			/*case sixty_four:
				return Value(bv64_ + v.bv64_);
			case thirty_two:
				return Value(bv32_ + v.bv32_);
			case sixteen:
				return Value(bv16_ + v.bv16_);
			case eight:
				return Value(bv8_ + v.bv8_);*/
			case one:
				return Value(bv1_ + v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	const Value Value::operator-(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
			/*case sixty_four:
				return Value(bv64_ - v.bv64_);
			case thirty_two:
				return Value(bv32_ - v.bv32_);
			case sixteen:
				return Value(bv16_ - v.bv16_);
			case eight:
				return Value(bv8_ - v.bv8_);*/
			case one:
				return Value(bv1_ - v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	const Value Value::operator*(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return Value(bv64_*v.bv64_);
		case thirty_two:
			return Value(bv32_*v.bv32_);
		case sixteen:
			return Value(bv16_*v.bv16_);
		case eight:
			return Value(bv8_*v.bv8_);*/
		case one:
			return Value(bv1_*v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	const Value Value::operator%(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
			/*case sixty_four:
				return Value(bv64_%v.bv64_);
			case thirty_two:
				return Value(bv32_%v.bv32_);
			case sixteen:
				return Value(bv16_%v.bv16_);
			case eight:
				return Value(bv8_%v.bv8_);*/
			case one:
				return Value(bv1_%v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	// Signed modulo operation which interprets this and v as signed
	const Value Value::smod(const Value& v) const {
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_) {
		/*case sixty_four: {
			LONG_INT num = bv64_.get_value();
			LONG_INT den = v.bv64_.get_value();
			LONG_INT smod = num%den;
			return Value(BV64(LONG_UINT(smod)));
		}
		case thirty_two: {
			INT num = bv32_.get_value();
			INT den = v.bv32_.get_value();
			INT smod = num%den;
			return Value(BV32(UINT(smod)));
		}
		case sixteen: {
			SHORT num = bv16_.get_value();
			SHORT den = v.bv16_.get_value();
			SHORT smod = num%den;
			return Value(BV16(USHORT(smod)));
		}
		case eight: {
			CHAR num = bv16_.get_value();
			CHAR den = v.bv16_.get_value();
			CHAR smod = num%den;
			return Value(BV8(UCHAR(smod)));
		}*/
		case one:
			return Value(bv1_%v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	const Value Value::operator/(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return Value(bv64_ / v.bv64_);
		case thirty_two:
			return Value(bv32_ / v.bv32_);
		case sixteen:
			return Value(bv16_ / v.bv16_);
		case eight:
			return Value(bv8_ / v.bv8_);*/
		case one:
			return Value(bv1_ / v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	// Signed division operation which interprets this and v as signed
	const Value Value::sdiv(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four: {
			LONG_INT num = bv64_.get_value();
			LONG_INT den = v.bv64_.get_value();
			LONG_INT sdiv = num / den;
			return Value(BV64(LONG_UINT(sdiv)));
		}
		case thirty_two: {
			INT num = bv32_.get_value();
			INT den = v.bv32_.get_value();
			INT sdiv = num / den;
			return Value(BV32(UINT(sdiv)));
		}
		case sixteen: {
			SHORT num = bv16_.get_value();
			SHORT den = v.bv16_.get_value();
			SHORT sdiv = num / den;
			return Value(BV16(USHORT(sdiv)));
		}
		case eight: {
			CHAR num = bv16_.get_value();
			CHAR den = v.bv16_.get_value();
			CHAR sdiv = num / den;
			return Value(BV8(UCHAR(sdiv)));
		}*/
		case one:
			return Value(bv1_ / v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	const Value Value::operator&(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return Value(bv64_&v.bv64_);
		case thirty_two:
			return Value(bv32_&v.bv32_);
		case sixteen:
			return Value(bv16_&v.bv16_);
		case eight:
			return Value(bv8_&v.bv8_);*/
		case one:
			return Value(bv1_&v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	const Value Value::operator|(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return Value(bv64_ | v.bv64_);
		case thirty_two:
			return Value(bv32_ | v.bv32_);
		case sixteen:
			return Value(bv16_ | v.bv16_);
		case eight:
			return Value(bv8_ | v.bv8_);*/
		case one:
			return Value(bv1_ | v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	const Value Value::operator^(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return Value(bv64_^v.bv64_);
		case thirty_two:
			return Value(bv32_^v.bv32_);
		case sixteen:
			return Value(bv16_^v.bv16_);
		case eight:
			return Value(bv8_^v.bv8_);*/
		case one:
			return Value(bv1_^v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	const Value Value::operator>>(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return Value(bv64_ >> v.bv64_);
		case thirty_two:
			return Value(bv32_ >> v.bv32_);
		case sixteen:
			return Value(bv16_ >> v.bv16_);
		case eight:
			return Value(bv8_ >> v.bv8_);*/
		case one:
			if (v.bv1_ == BV1(true))
				return Value(BV1(false));
			else
				return Value(bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	const Value Value::operator<<(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return Value(bv64_ << v.bv64_);
		case thirty_two:
			return Value(bv32_ << v.bv32_);
		case sixteen:
			return Value(bv16_ << v.bv16_);
		case eight:
			return Value(bv8_ << v.bv8_);*/
		case one:
			if (v.bv1_ == BV1(true))
				return Value(BV1(false));
			else
				return Value(bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return v;
	}

	bool Value::operator< (const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return bv64_<v.bv64_;
		case thirty_two:
			return bv32_<v.bv32_;
		case sixteen:
			return bv16_<v.bv16_;
		case eight:
			return bv8_<v.bv8_;*/
		case one:
			return bv1_<v.bv1_;
		}

		UWAssert::UWAssert::shouldNeverHappen();
		return false;
	}

	bool Value::operator>(const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return bv64_>v.bv64_;
		case thirty_two:
			return bv32_>v.bv32_;
		case sixteen:
			return bv16_>v.bv16_;
		case eight:
			return bv8_>v.bv8_;*/
		case one:
			return bv1_>v.bv1_;
		}

		UWAssert::UWAssert::shouldNeverHappen();
		return false;
	}

	bool Value::operator<= (const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return bv64_ <= v.bv64_;
		case thirty_two:
			return bv32_ <= v.bv32_;
		case sixteen:
			return bv16_ <= v.bv16_;
		case eight:
			return bv8_ <= v.bv8_;*/
		case one:
			return bv1_ <= v.bv1_;
		}

		UWAssert::UWAssert::shouldNeverHappen();
		return false;
	}

	bool Value::operator>= (const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return bv64_ >= v.bv64_;
		case thirty_two:
			return bv32_ >= v.bv32_;
		case sixteen:
			return bv16_ >= v.bv16_;
		case eight:
			return bv8_ >= v.bv8_;*/
		case one:
			return bv1_ >= v.bv1_;
		}

		UWAssert::UWAssert::shouldNeverHappen();
		return false;
	}

	bool Value::operator== (const Value& v) const
	{
		UWAssert::UWAssert::shouldNeverHappen(size_ != v.size_);
		switch (size_)
		{
		/*case sixty_four:
			return bv64_ == v.bv64_;
		case thirty_two:
			return bv32_ == v.bv32_;
		case sixteen:
			return bv16_ == v.bv16_;
		case eight:
			return bv8_ == v.bv8_;*/
		case one:
			return bv1_ == v.bv1_;
		}

		UWAssert::UWAssert::shouldNeverHappen();
		return false;
	}

	bool Value::operator!= (const Value& v) const
	{
		return !(*this == v);
	}

	Value& Value::operator=(const Value& v)
	{
		if (this != &v) {
			size_ = v.size_;
			/*bv64_ = v.bv64_;
			bv32_ = v.bv32_;
			bv16_ = v.bv16_;
			bv8_ = v.bv8_;*/
			bv1_ = v.bv1_;
		}
		return *this;
	}

	Value Value::ZeroExtend(Bitsize new_bitsize) const
	{
		LONG_UINT new_bitwidth = BitVector::GetBitWidth(new_bitsize);
		LONG_UINT old_bitwidth = BitVector::GetBitWidth(size_);
		UWAssert::UWAssert::shouldNeverHappen(old_bitwidth > new_bitwidth);
		switch (size_)
		{
		/*case sixty_four:
			return Value(new_bitsize, bv64_.get_value());
		case thirty_two:
			return Value(new_bitsize, LONG_UINT(bv32_.get_value()));
		case sixteen:
			return Value(new_bitsize, LONG_UINT(bv16_.get_value()));
		case eight:
			return Value(new_bitsize, LONG_UINT(bv8_.get_value()));*/
		case one:
			return Value(new_bitsize, LONG_UINT(bv1_.get_value()));
		}

		UWAssert::UWAssert::shouldNeverHappen();

		return Value(new_bitsize);
	}

	Value Value::SignExtend(Bitsize new_bitsize) const
	{
		LONG_UINT new_bitwidth = BitVector::GetBitWidth(new_bitsize);
		LONG_UINT old_bitwidth = BitVector::GetBitWidth(size_);
		UWAssert::UWAssert::shouldNeverHappen(old_bitwidth > new_bitwidth);
		switch (size_)
		{
		/*case sixty_four:
			return Value(new_bitsize, bv64_.get_value());
		case thirty_two:
			return Value(new_bitsize, LONG_UINT(LONG_INT(INT(bv32_.get_value()))));
		case sixteen:
			return Value(new_bitsize, LONG_UINT(LONG_INT(SHORT(bv16_.get_value()))));
		case eight:
			return Value(new_bitsize, LONG_UINT(LONG_INT(CHAR(bv8_.get_value()))));*/
		case one:
			return Value(new_bitsize, LONG_UINT(LONG_INT(bool(bv1_.get_value()))));
		}

		UWAssert::UWAssert::shouldNeverHappen();

		return Value(new_bitsize);
	}

	Value Value::Trunc(Bitsize new_bitsize) const
	{
		LONG_UINT new_bitwidth = BitVector::GetBitWidth(new_bitsize);
		LONG_UINT old_bitwidth = BitVector::GetBitWidth(size_);
		UWAssert::UWAssert::shouldNeverHappen(old_bitwidth < new_bitwidth);
		switch (size_)
		{
		/*case sixty_four:
			return Value(new_bitsize, bv64_.get_value());
		case thirty_two:
			return Value(new_bitsize, LONG_UINT(bv32_.get_value()));
		case sixteen:
			return Value(new_bitsize, LONG_UINT(bv16_.get_value()));
		case eight:
			return Value(new_bitsize, LONG_UINT(bv8_.get_value()));*/
		case one:
			return Value(new_bitsize, LONG_UINT(bv1_.get_value()));
		}

		UWAssert::UWAssert::shouldNeverHappen();

		return Value(new_bitsize);
	}

	Value Value::SignExtendOrTrunc(Bitsize new_bitsize) const
	{
		LONG_UINT new_bitwidth = BitVector::GetBitWidth(new_bitsize);
		LONG_UINT old_bitwidth = BitVector::GetBitWidth(size_);
		if (old_bitwidth == new_bitwidth)
			return Value(*this);

		if (old_bitwidth < new_bitwidth)
			return SignExtend(new_bitsize);

		return Trunc(new_bitsize);
	}

	Bitsize Value::GetBitsize() const
	{
		return size_;
	}

	size_t compute_value_rank(Value v)
	{
		switch (v.size_)
		{
		/*case sixty_four:
			return compute_rank(v.bv64_);
		case thirty_two:
			return compute_rank(v.bv32_);
		case sixteen:
			return compute_rank(v.bv16_);
		case eight:
			return compute_rank(v.bv8_);*/
		case one:
			return compute_rank(v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return 0;
	}

	Value multiplicative_inverse_value(Value v)
	{
		switch (v.size_)
		{
		/*case sixty_four:
			return multiplicative_inverse(v.bv64_);
		case thirty_two:
			return multiplicative_inverse(v.bv32_);
		case sixteen:
			return multiplicative_inverse(v.bv16_);
		case eight:
			return multiplicative_inverse(v.bv8_);*/
		case one:
			return multiplicative_inverse(v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return BV1(0);
		//return BV32(0);
	}

	std::ostream& operator<<(std::ostream& o, const Value v)
	{
		switch (v.size_)
		{
		/*case sixty_four:
			o << v.bv64_;
			break;
		case thirty_two:
			o << v.bv32_;
			break;
		case sixteen:
			o << v.bv16_;
			break;
		case eight:
			o << v.bv8_;
			break;*/
		case one:
			o << v.bv1_;
		}
		return o;
	}

	size_t GetMinPosOrNegPopCount(Value v)
	{
		switch (v.size_)
		{
		/*case sixty_four:
			return GetMinPosOrNegPopCount<BV64>(v.bv64_);
		case thirty_two:
			return GetMinPosOrNegPopCount<BV32>(v.bv32_);
		case sixteen:
			return GetMinPosOrNegPopCount<BV16>(v.bv16_);
		case eight:
			return GetMinPosOrNegPopCount<BV8>(v.bv8_);*/
		case one:
			return GetMinPosOrNegPopCount<BV1>(v.bv1_);
		}
		UWAssert::UWAssert::shouldNeverHappen();
		return 0;
	}

}
