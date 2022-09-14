#ifndef ZERO_INDICES_MAP_GUARD
#define ZERO_INDICES_MAP_GUARD

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
#include <vector>

#include "hashset.h"

class ZeroIndicesMapHandle;
class ZeroIndicesMapBody;

//***************************************************************
// ZeroIndicesMapHandle
//***************************************************************

class ZeroIndicesMapHandle {
public:
	ZeroIndicesMapHandle();                               // Default constructor
	~ZeroIndicesMapHandle();                              // Destructor
	ZeroIndicesMapHandle(const ZeroIndicesMapHandle &r);             // Copy constructor
	ZeroIndicesMapHandle(const int size);             // constructor with size for vector
	ZeroIndicesMapHandle& operator= (const ZeroIndicesMapHandle &r); // Overloaded assignment
	bool operator!= (const ZeroIndicesMapHandle &r) const;      // Overloaded !=
	bool operator== (const ZeroIndicesMapHandle &r) const;      // Overloaded ==
	int operator[](int i);                        // Overloaded []
	unsigned int Hash(unsigned int modsize);
	void Add_BIndex(const int i);
	bool Member(int i);
	int Lookup(int x);
	int Get_AIndex();
	void Set_AIndex(int i);
	unsigned int Size();
	void Canonicalize();
	ZeroIndicesMapBody *mapContents;
	static Hashset<ZeroIndicesMapBody> *canonicalZeroIndicesMapBodySet;
	std::ostream& print(std::ostream & out = std::cout) const;
};

std::ostream& operator<< (std::ostream & out, const ZeroIndicesMapHandle &r);
extern std::size_t hash_value(const ZeroIndicesMapHandle& val);

//***************************************************************
// ZeroIndicesMapBody
//***************************************************************

class ZeroIndicesMapBody {

	friend void ZeroIndicesMapHandle::Canonicalize();
	friend unsigned int ZeroIndicesMapHandle::Hash(unsigned int modsize);

public:
	ZeroIndicesMapBody();    // Constructor
	void IncrRef();
	void DecrRef();
	unsigned int Hash(unsigned int modsize);
	void setHashCheck();
	unsigned int refCount;         // reference-count value
	std::vector<int> b_indices;
	int a_index;
	bool operator==(const ZeroIndicesMapBody &o) const;
	int operator[](int i);                        // Overloaded []
	unsigned int hashCheck;
	bool isCanonical;              // Is this ZeroIndicesMapBody in *canonicalZeroIndicesMapBodySet?

};

#endif
