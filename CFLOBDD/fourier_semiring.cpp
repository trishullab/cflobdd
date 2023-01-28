//
//    Copyright (c) 1999 Thomas W. Reps
//    All Rights Reserved.
//
//    This software is furnished under a license and may be used and
//    copied only in accordance with the terms of such license and the
//    inclusion of the above copyright notice.  This software or any
//    other copies thereof or any derivative works may not be provided
//    or otherwise made available to any other person.  Title to and
//    ownership of the software and any derivative works is retained
//    by Thomas W. Reps.
//
//    THIS IMPLEMENTATION MAY HAVE BUGS, SOME OF WHICH MAY HAVE SERIOUS
//    CONSEQUENCES.  THOMAS W. REPS PROVIDES THIS SOFTWARE IN ITS "AS IS"
//    CONDITION, AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
//    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//    AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
//    THOMAS W. REPS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "fourier_semiring.h"

#define THRESHOLD 1e-5

fourierSemiring::fourierSemiring()
	: val(0), ringSize(1), isComplexValueSet(false)
{
	complex_value = BIG_COMPLEX(0);
}

fourierSemiring::fourierSemiring(const BIG_INT i1, const BIG_INT i2)
{
	isComplexValueSet = false;
	complex_value = BIG_COMPLEX(0);
	// ringSize = i2;
	// val = i1;
	// if (i1 < 0)
	// 	val = i2 + i1;
	// if (ringSize != 1)
	// 	val = val % ringSize;
	// if (val == 0 && ringSize != 1)
	// {
	// 	val = 1;
	// 	ringSize = 1;
	// }
	SetValAndRingSize(i1, i2);
}

void fourierSemiring::SetValAndRingSize(BIG_INT v, BIG_INT r) 
{ 
	ringSize = r;
	val = v;
	if (v < 0) 
		val = r + v; 
	if (ringSize != 1)
		val = val % ringSize;
	if (val == 0 && ringSize != 1)
	{
		val = 1;
		ringSize = 1;
	}
	// BIG_INT g = boost::math::gcd(val, ringSize);
	// val = val / g;
	// ringSize = ringSize / g;
}

std::ostream& operator<< (std::ostream & out, const fourierSemiring &p)
{
	if (!p.isComplexValueSet)
		out << "(" << p.GetVal() << ", " << p.GetRingSize() << ") NC";
	else
		out << p.complex_value << " C";
	return(out);
}

BIG_COMPLEX truncateToZeroIfNecessary(BIG_COMPLEX c)
{
	auto real = c.real();
	auto imag = c.imag();
	if (abs(real) < THRESHOLD)
		real = 0;
	if (abs(imag) < THRESHOLD)
		imag = 0;
	return BIG_COMPLEX(real, imag);
}

bool operator==(const fourierSemiring& lhs, const fourierSemiring& rhs)
{
	if (lhs.isComplexValueSet == rhs.isComplexValueSet && !lhs.isComplexValueSet)
		return (lhs.GetVal() == rhs.GetVal()) && (lhs.GetRingSize() == rhs.GetRingSize());
	else if (lhs.isComplexValueSet == rhs.isComplexValueSet && lhs.isComplexValueSet)
		return lhs.complex_value == rhs.complex_value;
	else if (lhs.isComplexValueSet)
	{
		double v = (double)(2 * rhs.GetVal()) / (double) rhs.GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX rhs_complex(cos_v, sin_v);
		rhs_complex = truncateToZeroIfNecessary(rhs_complex);
		if (rhs.GetVal() == 0 && rhs.GetRingSize() == 1)
			rhs_complex = BIG_COMPLEX(0);
		return lhs.complex_value == rhs_complex;
	}
	else
	{
		double v = (double)(2 * lhs.GetVal()) / (double) lhs.GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX lhs_complex(cos_v, sin_v);
		lhs_complex = truncateToZeroIfNecessary(lhs_complex);
		if (lhs.GetVal() == 0 && lhs.GetRingSize() == 1)
			lhs_complex = BIG_COMPLEX(0);
		return rhs.complex_value == lhs_complex;
	}
	abort();
}

fourierSemiring& fourierSemiring::operator= (const fourierSemiring& i)
{
	if (this != &i)      // don't assign to self!
	{
		val = i.GetVal();
		ringSize = i.GetRingSize();
		isComplexValueSet = i.isComplexValueSet;
		complex_value = i.complex_value;
	}
	return *this;
}

fourierSemiring fourierSemiring::operator* (const fourierSemiring& p){
	fourierSemiring ans;
	if (isComplexValueSet == p.isComplexValueSet && !isComplexValueSet)
	{
		if (p.GetRingSize() == 1){
			if (p.GetVal() == 0) { return ans; }
			else if (p.GetVal() == 1) { 
				ans.SetValAndRingSize(val, ringSize);
				return ans;
			}
			else {
				std::cout << *(this) << std::endl;
				std::cout << p << std::endl;
				abort();
			}
		}
		else if (ringSize == 1){
			if (val == 0) return ans;
			else if (val == 1){ 
				ans.SetValAndRingSize(p.GetVal(), p.GetRingSize());
				return ans;
			}
			else { 
				std::cout << *(this) << std::endl;
				std::cout << p << std::endl;
				abort();
			}
		}
		else if (p.GetRingSize() > ringSize){
			BIG_INT factor = p.GetRingSize() / ringSize;
			ans.SetValAndRingSize( (factor * val + p.GetVal()), p.GetRingSize());
		}
		else if (p.GetRingSize() < ringSize){
			BIG_INT factor = ringSize / p.GetRingSize();
			ans.SetValAndRingSize((factor * p.GetVal() + val), ringSize);
		}
		else if (p.GetRingSize() == ringSize){
			ans.SetValAndRingSize((val + p.GetVal()), ringSize);
		}
		else{
			std::cout << *(this) << std::endl;
			std::cout << p << std::endl;
			abort();
		}
	}
	else if (isComplexValueSet == p.isComplexValueSet && isComplexValueSet)
	{
		ans.complex_value = complex_value * p.complex_value;
		ans.complex_value = truncateToZeroIfNecessary(ans.complex_value);
		if (ans.complex_value == BIG_COMPLEX(0))
		{
			ans = fourierSemiring(0, 1);
			ans.isComplexValueSet = false;
		}
		else
			ans.isComplexValueSet = true;
	}
	else if (isComplexValueSet)
	{
		double v = (double)(2 * p.GetVal()) / (double) p.GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX rhs_complex(cos_v, sin_v);
		if (p.GetVal() == 0 && p.GetRingSize() == 1)
			rhs_complex = BIG_COMPLEX(0);
		ans.complex_value = complex_value * rhs_complex;
		ans.complex_value = truncateToZeroIfNecessary(ans.complex_value);
		if (ans.complex_value == BIG_COMPLEX(0))
		{
			ans = fourierSemiring(0, 1);
			ans.isComplexValueSet = false;
		}
		else
			ans.isComplexValueSet = true;
	}
	else
	{
		double v = (double)(2 * GetVal()) / (double) GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX lhs_complex(cos_v, sin_v);
		if (GetVal() == 0 && GetRingSize() == 1)
			lhs_complex = BIG_COMPLEX(0);
		ans.complex_value = lhs_complex * p.complex_value;
		ans.complex_value = truncateToZeroIfNecessary(ans.complex_value);
		if (ans.complex_value == BIG_COMPLEX(0))
		{
			ans = fourierSemiring(0, 1);
			ans.isComplexValueSet = false;
		}
		else
			ans.isComplexValueSet = true;
	}
	return ans;
}

// curr/p
fourierSemiring fourierSemiring::operator/ (const fourierSemiring& p){
	fourierSemiring ans;
	if (isComplexValueSet == p.isComplexValueSet && !isComplexValueSet)
	{
		if (p.GetRingSize() == 1){
			if (p.GetVal() == 0) { 
				std::cout << *(this) << std::endl;
				std::cout << p << std::endl;
				abort(); 
			}
			else if (p.GetVal() == 1) { 
				ans.SetValAndRingSize(val, ringSize);
				return ans;
			}
			else {
				std::cout << *(this) << std::endl;
				std::cout << p << std::endl;
				abort();
			}
		}
		else if (ringSize == 1){
			if (val == 0) return ans;
			else if (val == 1){ 
				ans.SetValAndRingSize(-1 * p.GetVal(), p.GetRingSize());
				return ans;
			}
			else {
				std::cout << *(this) << std::endl;
				std::cout << p << std::endl;
				abort();
			}
		}
		else if (p.GetRingSize() > ringSize){
			BIG_INT factor = p.GetRingSize() / ringSize;
			ans.SetValAndRingSize(( factor * val - p.GetVal()), p.GetRingSize());
		}
		else if (p.GetRingSize() < ringSize){
			BIG_INT factor = ringSize / p.GetRingSize();
			ans.SetValAndRingSize((factor * p.GetVal() - val), ringSize);
		}
		else if (p.GetRingSize() == ringSize){
			ans.SetValAndRingSize((val - p.GetVal()), ringSize);
		}
		else{
			std::cout << *(this) << std::endl;
			std::cout << p << std::endl;
			abort();
		}
	}
	else if (isComplexValueSet == p.isComplexValueSet && isComplexValueSet)
	{
		ans.complex_value = complex_value / p.complex_value;
		ans.complex_value = truncateToZeroIfNecessary(ans.complex_value);
		if (ans.complex_value == BIG_COMPLEX(0))
		{
			ans = fourierSemiring(0, 1);
			ans.isComplexValueSet = false;
		}
		else
			ans.isComplexValueSet = true;
	}
	else if (isComplexValueSet)
	{
		double v = (double)(2 * p.GetVal()) / (double) p.GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX rhs_complex(cos_v, sin_v);
		if (p.GetVal() == 0 && p.GetRingSize() == 1)
			rhs_complex = BIG_COMPLEX(0);
		ans.complex_value = complex_value / rhs_complex;
		ans.complex_value = truncateToZeroIfNecessary(ans.complex_value);
		if (ans.complex_value == BIG_COMPLEX(0))
		{
			ans = fourierSemiring(0, 1);
			ans.isComplexValueSet = false;
		}
		else
			ans.isComplexValueSet = true;
	}
	else 
	{
		double v = (double)(2 * GetVal()) / (double) GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX lhs_complex(cos_v, sin_v);
		if (GetVal() == 0 && GetRingSize() == 1)
			lhs_complex = BIG_COMPLEX(0);
		ans.complex_value = lhs_complex / p.complex_value;
		ans.complex_value = truncateToZeroIfNecessary(ans.complex_value);
		if (ans.complex_value == BIG_COMPLEX(0))
		{
			ans = fourierSemiring(0, 1);
			ans.isComplexValueSet = false;
		}
		else
			ans.isComplexValueSet = true;
	}
	return ans;
}

fourierSemiring fourierSemiring::operator+ (const fourierSemiring& p){
	fourierSemiring ans;
	if (isComplexValueSet == p.isComplexValueSet && !isComplexValueSet)
	{
		if (p.GetRingSize() == 1 && p.GetVal() == 0){ 
			ans.SetVal(val);
			ans.SetRingSize(ringSize);
			return ans;
		}
		if (ringSize == 1 && val == 0){ return p; }
		
		double v1 = (double)(2 * p.GetVal()) / (double) p.GetRingSize();
		auto cos_v1 = boost::math::cos_pi(v1);
		auto sin_v1 = boost::math::sin_pi(v1);
		BIG_COMPLEX rhs_complex(cos_v1, sin_v1);
		if (p.GetVal() == 0 && p.GetRingSize() == 1)
			rhs_complex = BIG_COMPLEX(0);

		double v2 = (double)(2 * GetVal()) / (double) GetRingSize();
		auto cos_v2 = boost::math::cos_pi(v2);
		auto sin_v2 = boost::math::sin_pi(v2);
		BIG_COMPLEX lhs_complex(cos_v2, sin_v2);
		if (GetVal() == 0 && GetRingSize() == 1)
			lhs_complex = BIG_COMPLEX(0);

		ans.complex_value = lhs_complex + rhs_complex;
		ans.complex_value = truncateToZeroIfNecessary(ans.complex_value);
		if (ans.complex_value == BIG_COMPLEX(0))
		{
			ans = fourierSemiring(0, 1);
			ans.isComplexValueSet = false;
		}
		else
			ans.isComplexValueSet = true;	
	}
	else if (isComplexValueSet == p.isComplexValueSet && isComplexValueSet)
	{
		ans.complex_value = complex_value + p.complex_value;
		ans.complex_value = truncateToZeroIfNecessary(ans.complex_value);
		if (ans.complex_value == BIG_COMPLEX(0))
		{
			ans = fourierSemiring(0, 1);
			ans.isComplexValueSet = false;
		}
		else
			ans.isComplexValueSet = true;
	}
	else if (isComplexValueSet)
	{
		double v = (double)(2 * p.GetVal()) / (double) p.GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX rhs_complex(cos_v, sin_v);
		if (p.GetVal() == 0 && p.GetRingSize() == 1)
			rhs_complex = BIG_COMPLEX(0);
		ans.complex_value = complex_value + rhs_complex;
		ans.complex_value = truncateToZeroIfNecessary(ans.complex_value);
		if (ans.complex_value == BIG_COMPLEX(0))
		{
			ans = fourierSemiring(0, 1);
			ans.isComplexValueSet = false;
		}
		else
			ans.isComplexValueSet = true;
	}
	else
	{
		double v = (double)(2 * GetVal()) / (double) GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX lhs_complex(cos_v, sin_v);
		if (GetVal() == 0 && GetRingSize() == 1)
			lhs_complex = BIG_COMPLEX(0);
		ans.complex_value = lhs_complex + p.complex_value;
		ans.complex_value = truncateToZeroIfNecessary(ans.complex_value);
		if (ans.complex_value == BIG_COMPLEX(0))
		{
			ans = fourierSemiring(0, 1);
			ans.isComplexValueSet = false;
		}
		else
			ans.isComplexValueSet = true;
	}
	return ans;
}


// Overloaded !=
bool fourierSemiring::operator!=(const fourierSemiring& p) const
{
	if (isComplexValueSet == p.isComplexValueSet && !isComplexValueSet)
		return (val != p.GetVal()) || (ringSize != p.GetRingSize());
	else if (isComplexValueSet == p.isComplexValueSet && isComplexValueSet)
		return complex_value != p.complex_value;
	else if (isComplexValueSet)
	{
		double v = (double)(2 * p.GetVal()) / (double) p.GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX rhs_complex(cos_v, sin_v);
		if (p.GetVal() == 0 && p.GetRingSize() == 1)
			rhs_complex = BIG_COMPLEX(0);
		rhs_complex = truncateToZeroIfNecessary(rhs_complex);
		return complex_value != rhs_complex;
	}
	else 
	{
		double v = (double)(2 * GetVal()) / (double) GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX lhs_complex(cos_v, sin_v);
		if (GetVal() == 0 && GetRingSize() == 1)
			lhs_complex = BIG_COMPLEX(0);
		lhs_complex = truncateToZeroIfNecessary(lhs_complex);
		return lhs_complex != p.complex_value;
	}
}

// Overloaded ==
bool fourierSemiring::operator==(const fourierSemiring& p)
{
	if (isComplexValueSet == p.isComplexValueSet && !isComplexValueSet)
		return (val == p.GetVal()) && (ringSize == p.GetRingSize());
	else if (isComplexValueSet == p.isComplexValueSet && isComplexValueSet)
		return complex_value == p.complex_value;
	else if (isComplexValueSet)
	{
		double v = (double)(2 * p.GetVal()) / (double) p.GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX rhs_complex(cos_v, sin_v);
		if (p.GetVal() == 0 && p.GetRingSize() == 1)
			rhs_complex = BIG_COMPLEX(0);
		rhs_complex = truncateToZeroIfNecessary(rhs_complex);
		return complex_value == rhs_complex;
	}
	else 
	{
		double v = (double)(2 * GetVal()) / (double) GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX lhs_complex(cos_v, sin_v);
		if (GetVal() == 0 && GetRingSize() == 1)
			lhs_complex = BIG_COMPLEX(0);
		lhs_complex = truncateToZeroIfNecessary(lhs_complex);
		return lhs_complex == p.complex_value;
	}
}

fourierSemiring operator*(const BIG_INT lhs, const fourierSemiring& rhs)
{
	fourierSemiring ans;
	if (lhs == 0)
		return ans;
	if (rhs.GetVal() == 0 && rhs.GetRingSize() == 1)
		return ans;
	if (lhs == 1){
		ans.SetVal(rhs.GetVal());
		ans.SetRingSize(rhs.GetRingSize());
		return ans;
	}
	std::cout << lhs << std::endl;
	std::cout << rhs << std::endl;
	abort();
}

void fourierSemiring::ResetComplexValue()
{
	if (isComplexValueSet)
	{
		if (complex_value == BIG_COMPLEX(0))
		{
			val = 0;
			ringSize = 1;
		}
		else if (complex_value == BIG_COMPLEX(1))
		{
			val = 1;
			ringSize = 1;
		}
		else{
			std::cout << "complex_value: " << complex_value << std::endl;
			abort();
		}
		isComplexValueSet = false;
	}
}

void fourierSemiring::SetComplexValue()
{
	if (!isComplexValueSet)
	{
		double v = (double)(2 * GetVal()) / (double) GetRingSize();
		auto cos_v = boost::math::cos_pi(v);
		auto sin_v = boost::math::sin_pi(v);
		BIG_COMPLEX lhs_complex(cos_v, sin_v);
		if (GetVal() == 0 && GetRingSize() == 1)
			lhs_complex = BIG_COMPLEX(0);
		complex_value = lhs_complex;
		complex_value = truncateToZeroIfNecessary(complex_value);
		isComplexValueSet = true;	
	}
}

fourierSemiring operator*(const fourierSemiring& lhs, const BIG_INT rhs)
{
	fourierSemiring ans;
	if (rhs == 0)
		return ans;
	if (lhs.GetVal() == 0 && lhs.GetRingSize() == 1)
		return ans;
	if (rhs == 1){
		ans.SetVal(lhs.GetVal());
		ans.SetRingSize(lhs.GetRingSize());
		return ans;
	}
	std::cout << lhs << std::endl;
	std::cout << rhs << std::endl;
	abort();
}

std::size_t hash_value(const fourierSemiring& p)
{
	if (!p.isComplexValueSet)
	{
		boost::hash<BIG_INT> boost_hash;
		return 117 * (boost_hash(p.GetVal()) + 1) + boost_hash(p.GetRingSize());
	}
	else
	{
		boost::hash<BIG_COMPLEX> boost_hash;
		return boost_hash(p.complex_value);
	}
}