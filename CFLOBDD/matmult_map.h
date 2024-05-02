#ifndef MATMULT_MAP_GUARD
#define MATMULT_MAP_GUARD

//
//    Copyright (c) 2017 Thomas W. Reps
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
#include <unordered_map>
#include <boost/multiprecision/cpp_int.hpp>
#include <map>

#include "hashset.h"

class MatMultMapHandle;
class MatMultMapBody;

typedef std::pair<long int, long int> INT_PAIR;

typedef boost::multiprecision::cpp_int VAL_TYPE;
//typedef unsigned long long int VAL_TYPE;


//***************************************************************
// MatMultMapHandle
//***************************************************************

class MatMultMapHandle {
public:
	MatMultMapHandle();                               // Default constructor
	~MatMultMapHandle();                              // Destructor
	MatMultMapHandle(const MatMultMapHandle &r);             // Copy constructor
	MatMultMapHandle(const int size);             // constructor with size for vector
	MatMultMapHandle& operator= (const MatMultMapHandle &r); // Overloaded assignment
	bool operator!= (const MatMultMapHandle &r) const;      // Overloaded !=
	bool operator== (const MatMultMapHandle &r) const;      // Overloaded ==
	VAL_TYPE& operator[](INT_PAIR& p);                        // Overloaded []
	unsigned int Hash(unsigned int modsize) const;
	void Add(const INT_PAIR& p, VAL_TYPE& v);
	void ForceAdd(const INT_PAIR& p, VAL_TYPE& v);
	bool Member(INT_PAIR& p);
	VAL_TYPE Lookup(INT_PAIR& x);
	std::string ToString();
	unsigned int Size();
	void Canonicalize();
	MatMultMapBody *mapContents;
	static Hashset<MatMultMapBody> *canonicalMatMultMapBodySet;
	std::ostream& print(std::ostream & out = std::cout) const;
	MatMultMapHandle operator+ (const MatMultMapHandle&) const; // map concatenation with summation as merge operation
};

std::ostream& operator<< (std::ostream & out, const MatMultMapHandle &r);

extern MatMultMapHandle operator* (const VAL_TYPE&, const MatMultMapHandle&);
extern MatMultMapHandle operator* (const MatMultMapHandle&, const VAL_TYPE&);
extern std::size_t hash_value(const MatMultMapHandle& val);

//***************************************************************
// MatMultMapBody
//***************************************************************

class MatMultMapBody {

	friend void MatMultMapHandle::Canonicalize();
	friend unsigned int MatMultMapHandle::Hash(unsigned int modsize) const;

public:
	MatMultMapBody();    // Constructor
	void IncrRef();
	void DecrRef();
	unsigned int Hash(unsigned int modsize) const;
	void setHashCheck();
	unsigned int refCount;         // reference-count value
	std::map<INT_PAIR, VAL_TYPE> map;
	bool operator==(const MatMultMapBody &o) const;
	VAL_TYPE& operator[](INT_PAIR& i);                        // Overloaded []
	long int hashCheck;
	bool contains_zero_val;
	bool isCanonical;              // Is this MatMultMapBody in *canonicalMatMultMapBodySet?

};

#endif
