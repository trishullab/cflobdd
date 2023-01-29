#ifndef FOURIER_SEMIRING_GUARD
#define FOURIER_SEMIRING_GUARD

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

#include <iostream>
#include <fstream>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/functional/hash.hpp>

typedef boost::multiprecision::cpp_int BIG_INT;
typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX;

class fourierSemiring {
public:
	fourierSemiring();                              // Default constructor
	fourierSemiring(const BIG_INT val, const BIG_INT ringSize);    // Constructor
	fourierSemiring& operator= (const fourierSemiring& p);  // Overloaded assignment
	bool operator!= (const fourierSemiring& p) const;     // Overloaded !=
	bool operator== (const fourierSemiring& p);     // Overloaded ==
	BIG_INT GetVal() const { return val; }     // Access function
	BIG_INT GetRingSize() const { return ringSize; }   // Access function
	void SetVal(BIG_INT v) { val = v; }
	void SetValAndRingSize(BIG_INT v, BIG_INT r);
	void SetRingSize(BIG_INT r){ ringSize = r; }
	fourierSemiring operator* (const fourierSemiring& p);
	fourierSemiring operator/ (const fourierSemiring& p);
	fourierSemiring operator+ (const fourierSemiring& p);
	bool isComplexValueSet;
	BIG_COMPLEX complex_value;
	void ResetComplexValue();
	void SetComplexValue();
	struct fourierSemiring_hash {
		size_t operator()(const fourierSemiring& p) const {
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
	};
private:
	BIG_INT val;
	BIG_INT ringSize;
};

bool operator==(const fourierSemiring& lhs, const fourierSemiring& rhs);
fourierSemiring operator*(const BIG_INT lhs, const fourierSemiring& rhs);
fourierSemiring operator*(const fourierSemiring& lhs, const BIG_INT rhs);
std::ostream& operator<< (std::ostream & out, const fourierSemiring &p);
std::size_t hash_value(const fourierSemiring& val);

#endif
