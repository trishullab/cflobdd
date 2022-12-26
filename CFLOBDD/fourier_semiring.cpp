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

fourierSemiring::fourierSemiring()
	: val(0), ringSize(1)
{
}

fourierSemiring::fourierSemiring(const BIG_INT i1, const BIG_INT i2)
{
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
}

std::ostream& operator<< (std::ostream & out, const fourierSemiring &p)
{
	out << "(" << p.GetVal() << ", " << p.GetRingSize() << ")";
	return(out);
}

bool operator==(const fourierSemiring& lhs, const fourierSemiring& rhs)
{
	return (lhs.GetVal() == rhs.GetVal()) && (lhs.GetRingSize() == rhs.GetRingSize());
}

fourierSemiring& fourierSemiring::operator= (const fourierSemiring& i)
{
	if (this != &i)      // don't assign to self!
	{
		val = i.GetVal();
		ringSize = i.GetRingSize();
	}
	return *this;
}

fourierSemiring fourierSemiring::operator* (const fourierSemiring& p){
	fourierSemiring ans;
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
	else if (p.GetRingSize() > p.GetVal()){
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
	return ans;
}

// curr/p
fourierSemiring fourierSemiring::operator/ (const fourierSemiring& p){
	fourierSemiring ans;
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
	return ans;
}

fourierSemiring fourierSemiring::operator+ (const fourierSemiring& p){
	fourierSemiring ans;
	if (p.GetRingSize() == 1 && p.GetVal() == 0){ 
		ans.SetVal(val);
		ans.SetRingSize(ringSize);
		return ans;
	}
	if (ringSize == 1 && val == 0){ return p; }
	
	std::cout << *(this) << std::endl;
	std::cout << p << std::endl;
	abort();
}


// Overloaded !=
bool fourierSemiring::operator!=(const fourierSemiring& p) const
{
	return (val != p.GetVal()) || (ringSize != p.GetRingSize());
}

// Overloaded ==
bool fourierSemiring::operator==(const fourierSemiring& p)
{
	return (val == p.GetVal()) && (ringSize == p.GetRingSize());
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
	boost::hash<BIG_INT> boost_hash;
	return 117 * (boost_hash(p.GetVal()) + 1) + boost_hash(p.GetRingSize());
}