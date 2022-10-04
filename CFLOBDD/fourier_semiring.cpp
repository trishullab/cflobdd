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

fourierSemiring::fourierSemiring(const int i1, const int i2)
	: val(i1), ringSize(i2)
{
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
			ans.SetRingSize(ringSize);
			ans.SetVal(val);
			return ans;
		}
		else abort();
	}
	else if (ringSize == 1){
		if (val == 0) return ans;
		else if (val == 1){ 
			ans.SetVal(p.GetVal());
			ans.SetRingSize(p.GetRingSize());
			return ans;
		}
		else abort();
	}
	else if (log2(p.GetRingSize()) - log2(ringSize) > 0){
		ans.SetVal((ringSize * val + p.GetVal()) % p.GetRingSize());
		ans.SetRingSize(p.GetRingSize());
	}
	else if (log2(p.GetRingSize()) - log2(ringSize) < 0){
		ans.SetVal((p.GetRingSize() * p.GetVal() + val) % ringSize);
		ans.SetRingSize(ringSize);
	}
	else if (log2(p.GetRingSize()) - log2(ringSize) == 0){
		ans.SetRingSize(ringSize);
		ans.SetVal((val + p.GetVal()) % ringSize);
	}
	else{
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

fourierSemiring operator*(const unsigned long long int lhs, const fourierSemiring& rhs)
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
	abort();
}

fourierSemiring operator*(const fourierSemiring& lhs, const unsigned long long int rhs)
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
	abort();
}

std::size_t hash_value(const fourierSemiring& p)
{
	return 117 * (p.GetVal() + 1) + p.GetRingSize();
}