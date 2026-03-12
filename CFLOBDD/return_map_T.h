#ifndef RETURN_MAP_GUARD
#define RETURN_MAP_GUARD

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
#include <unordered_set>
#include <vector>
#include "list_T.h"
#include "list_TPtr.h"
#include "hashset.h"
#include "reduction_map.h"
#include "intpair.h"
#include <complex>
#include <cstdint>
#include <typeinfo>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif
//#include <boost/multiprecision/cpp_int.hpp>
//#include "hash_functions.h"

inline void printProcessMemoryUsage() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        std::cout << "WorkingSetSize: " << (pmc.WorkingSetSize / (1024*1024)) << " MB"
                  << ", PeakWorkingSetSize: " << (pmc.PeakWorkingSetSize / (1024*1024)) << " MB"
                  << std::endl;
    }
#endif
}

template <typename T> class ReturnMapHandle;
template <typename T> class ReturnMapBody;

// Content-based hash and equality functors for ReturnMapBody* pointers,
// used by the std::unordered_set canonical store.
template <typename T>
struct ReturnMapBodyPtrHash {
    size_t operator()(ReturnMapBody<T>* p) const { return p->Hash(); }
};

template <typename T>
struct ReturnMapBodyPtrEq {
    bool operator()(ReturnMapBody<T>* a, ReturnMapBody<T>* b) const { return a == b || *a == *b; }
};

//using namespace boost::multiprecision;

//***************************************************************
// ReturnMapHandle
//***************************************************************

template <typename T>
class ReturnMapHandle {
 public:
  ReturnMapHandle();                               // Default constructor
  ReturnMapHandle(unsigned int capacity);
  ~ReturnMapHandle();                              // Destructor
  ReturnMapHandle(const ReturnMapHandle<T> &r);    // Copy constructor
  ReturnMapHandle<T>& operator= (const ReturnMapHandle<T> &r); // Overloaded assignment
  bool operator!= (const ReturnMapHandle<T> &r);      // Overloaded !=
  bool operator== (const ReturnMapHandle<T> &r);      // Overloaded ==
  T& operator[](unsigned int i);                       // Overloaded []
  size_t Hash();
  unsigned int Size();
  void AddToEnd(T y);
  bool Member(T y);
  T Lookup(int x);
  intpair convertValue(int y);
  int LookupInv(T y);
  void Canonicalize();
  ReturnMapHandle<T> Complement();
  ReturnMapHandle<T> Compose(ReductionMapHandle redMapHandle);
  void InducedReductionAndReturnMap(ReductionMapHandle &inducedReductionMapHandle,
	                                ReturnMapHandle<T> &inducedReturnMapHandle);
  ReturnMapBody<T> *mapContents;
  using CanonicalReturnMapBodySet = std::unordered_set<ReturnMapBody<T>*,
                                                        ReturnMapBodyPtrHash<T>,
                                                        ReturnMapBodyPtrEq<T>>;
  static CanonicalReturnMapBodySet *canonicalReturnMapBodySet;
  std::ostream& print(std::ostream & out = std::cout) const;
 private:
  static CanonicalReturnMapBodySet *initCanonicalSet();
};

template <typename T>
std::ostream& operator<< (std::ostream & out, const ReturnMapHandle<T> &r);

//***************************************************************
// ReturnMapBody
//***************************************************************

template <typename T>
class ReturnMapBody {

  friend void ReturnMapHandle<T>::Canonicalize();
  friend size_t ReturnMapHandle<T>::Hash();

 public:
  ReturnMapBody();    // Constructor
  ReturnMapBody(unsigned int capacity);    // Constructor
  //~ReturnMapBody();
  void IncrRef();
  void DecrRef();
  size_t Hash();
  void setHashCheck();
  unsigned int refCount;         // reference-count value
  std::vector<T> mapArray;
  bool operator==(const ReturnMapBody &o) const;
  T& operator[](unsigned int i);                       // Overloaded []
  unsigned int hashCheck;

  static ReturnMapBody<T>* Create();
  static ReturnMapBody<T>* Create(unsigned int capacity);

 protected:
  bool isCanonical;              // Is this ReturnMapBody in *canonicalReturnMapBodySet?

 private:
  // Heap-allocated so it is never destroyed at program exit, avoiding static
  // destruction-order issues when DecrRef is called from late static destructors.
  static std::vector<ReturnMapBody<T>*>& getFreeList() {
    static std::vector<ReturnMapBody<T>*>* const freeList = new std::vector<ReturnMapBody<T>*>();
    return *freeList;
  }
};


//***************************************************************
// ReturnMapBody
//***************************************************************

// Constructor
template <typename T>
ReturnMapBody<T>::ReturnMapBody()
  : refCount(0), isCanonical(false)
{
	  hashCheck = (unsigned int) NULL;
	  //IncrRef();
}

template <typename T>
ReturnMapBody<T>::ReturnMapBody(unsigned int capacity)
	: refCount(0), isCanonical(false)
{
	hashCheck = NULL;
	mapArray.reserve(capacity);
}

template <typename T>
void ReturnMapBody<T>::IncrRef()
{
  refCount++;    // Warning: Saturation not checked
}

template <typename T>
void ReturnMapBody<T>::DecrRef()
{
  if (--refCount == 0) {    // Warning: Saturation not checked
    if (isCanonical) {
      ReturnMapHandle<T>::canonicalReturnMapBodySet->erase(this);
    }
    mapArray.clear();
    hashCheck = 0;
    isCanonical = false;
    getFreeList().push_back(this);
  }
}

template <typename T>
ReturnMapBody<T>* ReturnMapBody<T>::Create()
{
  auto& fl = getFreeList();
  if (!fl.empty()) {
    ReturnMapBody<T>* p = fl.back();
    fl.pop_back();
    return p;
  }
  return new ReturnMapBody<T>();
}

template <typename T>
ReturnMapBody<T>* ReturnMapBody<T>::Create(unsigned int capacity)
{
  auto& fl = getFreeList();
  if (!fl.empty()) {
    ReturnMapBody<T>* p = fl.back();
    fl.pop_back();
    if (p->mapArray.capacity() < capacity)
      p->mapArray.reserve(capacity);
    return p;
  }
  return new ReturnMapBody<T>(capacity);
}

template <typename T>
bool ReturnMapBody<T>::operator==(const ReturnMapBody<T> &o) const
{
	if (hashCheck != o.hashCheck)
	{
		return false;
	}
	else if (mapArray.size() != o.mapArray.size())
	{
		return false;
	} else {
	  for (unsigned i = 0; i < mapArray.size(); i++)
	  {
		  if (mapArray[i] != o.mapArray[i])
		  {
			  return false;
		  }
	  }
	}
	return true;
}

// Overloaded []
template <typename T>
T& ReturnMapBody<T>::operator[](unsigned int i)
{
	return mapArray[i];
}


template <typename T>
std::ostream& operator<< (std::ostream & out, const ReturnMapBody<T> &r)
{
	out << "{RMB<T>: ";
	size_t last = r.mapArray.size() - 1;
	for(size_t i = 0; i <= last; ++i) {
		out << r.mapArray[i];
        if (i != last) 
            out << ", ";
    }
	out << " RMB<T>}";
    return out;
}

//***************************************************************
// ReturnMapHandle
//***************************************************************

// Factory method for static canonical-set member
template <typename T>
typename ReturnMapHandle<T>::CanonicalReturnMapBodySet *ReturnMapHandle<T>::initCanonicalSet()
{
    auto *s = new CanonicalReturnMapBodySet(RETURN_MAP_NUM_BUCKETS);
    s->max_load_factor(0.8f);
    return s;
}

// Initializations of static members ---------------------------------
template <typename T> typename ReturnMapHandle<T>::CanonicalReturnMapBodySet
    *ReturnMapHandle<T>::canonicalReturnMapBodySet = ReturnMapHandle<T>::initCanonicalSet();

// Default constructor
template <typename T>
ReturnMapHandle<T>::ReturnMapHandle()
  :  mapContents(ReturnMapBody<T>::Create())
{
  mapContents->IncrRef();
}

template <typename T>
ReturnMapHandle<T>::ReturnMapHandle(unsigned int capacity)
	: mapContents(ReturnMapBody<T>::Create(capacity))
{
	mapContents->IncrRef();
}

// Destructor
template <typename T>
ReturnMapHandle<T>::~ReturnMapHandle()
{
  mapContents->DecrRef();
}

// Copy constructor
template <typename T>
ReturnMapHandle<T>::ReturnMapHandle(const ReturnMapHandle<T> &r)
  :  mapContents(r.mapContents)
{
  mapContents->IncrRef();
}

// Overloaded assignment
template <typename T>
ReturnMapHandle<T>& ReturnMapHandle<T>::operator= (const ReturnMapHandle<T> &r)
{
  if (this != &r)      // don't assign to self!
  {
    ReturnMapBody<T> *temp = mapContents;
    mapContents = r.mapContents;
    mapContents->IncrRef();
    temp->DecrRef();
  }
  return *this;        
}

template <typename T>
intpair ReturnMapHandle<T>::convertValue(int y)
{
		return *(new intpair(y,y));
}

// Overloaded !=
template <typename T>
bool ReturnMapHandle<T>::operator!=(const ReturnMapHandle<T> &r)
{
  return (mapContents != r.mapContents);
}

// Overloaded ==
template <typename T>
bool ReturnMapHandle<T>::operator==(const ReturnMapHandle<T> &r)
{

  return (mapContents == r.mapContents);
}

// Overloaded []
template <typename T>
T& ReturnMapHandle<T>::operator[](unsigned int i)
{
	return (*(this->mapContents))[i];
}

// print
template <typename T>
std::ostream& ReturnMapHandle<T>::print(std::ostream & out) const
{
  out << *mapContents << std::endl;
  return out;
}

template <typename T>
std::ostream& operator<< (std::ostream & out, const ReturnMapHandle<T> &r)
{
  r.print(out);
  return(out);
}

template <typename T>
size_t ReturnMapHandle<T>::Hash()
{
	if (!(mapContents->isCanonical)) {
		std::cout << "Hash of a non-canonical ReturnMapHandle occurred" << std::endl;
		abort();
		this->Canonicalize();
	}
	assert(mapContents->isCanonical);
	return reinterpret_cast<uintptr_t>(mapContents) >> PTR_ALIGN_SHIFT;
}

template <typename T>
unsigned int ReturnMapHandle<T>::Size()
{
	try{
		return mapContents->mapArray.size();
	}
	catch (const std::exception& e){
		std::cout << "Exception type: " << typeid(e).name() << std::endl;
		std::cout << e.what() << std::endl;
		std::cout << "Size error in ReturnMapHandle<" << typeid(T).name() << ">" << std::endl;
		throw;
	}
}

template <typename T>
void ReturnMapHandle<T>::AddToEnd(T y)
{
	try{
		assert(mapContents->refCount <= 1);
		mapContents->mapArray.push_back(y);
	}
	catch (const std::exception& e){
		std::cout << "Exception type: " << typeid(e).name() << std::endl;
		std::cout << e.what() << std::endl;
		std::cout << "AddToEnd in ReturnMapHandle<" << typeid(T).name() << ">" << std::endl;
		std::cout << "refCount: " << mapContents->refCount << ", Size: " << Size() << std::endl;
		printProcessMemoryUsage();
		throw;
	}
}

template <typename T>
bool ReturnMapHandle<T>::Member(T y)
{
	for (unsigned i = 0; i < mapContents->mapArray.size(); i++)
	{
		if (mapContents->mapArray[i] == y) {
			return true;
		}
	}
	return false;
}

template <typename T>
T ReturnMapHandle<T>::Lookup(int x)
{
	return mapContents->mapArray[x];
}

template <typename T>
int ReturnMapHandle<T>::LookupInv(T y)
{
	for (unsigned i = 0; i < mapContents->mapArray.size(); i++)
	{
		if (mapContents->mapArray[i] == y)
		{
			return i;
		}
	}
	return -1;
}

template <typename T>
void ReturnMapHandle<T>::Canonicalize()
{ 
	try{
		ReturnMapBody<T> *answerContents;

		if (!mapContents->isCanonical) {
			mapContents->setHashCheck();
			auto it = canonicalReturnMapBodySet->find(mapContents);
			if (it == canonicalReturnMapBodySet->end()) {
				canonicalReturnMapBodySet->insert(mapContents);
				mapContents->isCanonical = true;
			}
			else {
				answerContents = *it;
				answerContents->IncrRef();
				mapContents->DecrRef();
				mapContents = answerContents;
			}
		}
	}
	catch (const std::exception& e){
		std::cout << "Exception type: " << typeid(e).name() << std::endl;
		std::cout << e.what() << std::endl;
		std::cout << "Canonicalize in ReturnMapHandle<" << typeid(T).name() << ">" << std::endl;
		throw;
	}
}

template <typename T>
ReturnMapHandle<T> ReturnMapHandle<T>::Compose(ReductionMapHandle redMapHandle)
{
  T c2, c3;
  ReturnMapHandle<T> answer;
  int size = mapContents->mapArray.size();
  for (int i = 0; i < size; i++)
  {
	  c2 = mapContents->mapArray[i];
	  c3 = redMapHandle.Lookup(c2);
	  answer.mapContents->mapArray.push_back(c3); 	  // Why not answer.AddToEnd(c3);
  }
  answer.Canonicalize();
  return answer;
}


template <typename T>
void ReturnMapHandle<T>::InducedReductionAndReturnMap(ReductionMapHandle &inducedReductionMapHandle,
	                                                  ReturnMapHandle<T> &inducedReturnMapHandle)
{
	try{
		unsigned int numEntries = 0;
		T c2;
		int d, e;

		int c1 = 0;
		int size = mapContents->mapArray.size();
		for (int i = 0; i < size; i++)
		{
			c2 = mapContents->mapArray[i];
			d = LookupInv(c2);
			if (d < c1) { // c1 and d are in same range-value equivalence class of this ReturnMap (i.e., [c2])
				e = inducedReductionMapHandle.Lookup(d);
				inducedReductionMapHandle.AddToEnd(e);
			}
			else { // c1 is the representative of a new range-value equivalence class of this ReturnMap (i.e., [c2])
				inducedReductionMapHandle.AddToEnd(numEntries);
				inducedReturnMapHandle.AddToEnd(c2);
				numEntries++;
			}
			c1++;
		}
		inducedReturnMapHandle.Canonicalize();
		inducedReductionMapHandle.Canonicalize();
	}
	catch (const std::exception& e){
		std::cout << "Exception type: " << typeid(e).name() << std::endl;
		std::cout << e.what() << std::endl;
		std::cout << "InducedReductionAndReturnMap in ReturnMapHandle<" << typeid(T).name() << ">" << std::endl;
		throw;
	}
}

// Left scalar-multiplication
template<typename T, typename T1>
ReturnMapHandle<T> operator*(T1 c, ReturnMapHandle<T> rmh)
{
	T v;
	ReturnMapHandle<T> answer;
	if (c == 1) return rmh;
	int size = rmh.mapContents->mapArray.size();
	for (int i = 0; i < size; i++)
	{
		v = rmh.mapContents->mapArray[i];
		T val = c * v;
		answer.AddToEnd(val);
		// Formerly, answer.mapContents->mapArray.push_back();  // When written as c * v, the compiler gave the message "'operator *' is ambiguous"
	}
	answer.Canonicalize();
	return answer;	
}

// Right scalar-multiplication
template<typename T>
ReturnMapHandle<T> operator*(ReturnMapHandle<T> rmh, int c)
{
	T v;
	ReturnMapHandle<T> answer;
	if (c == 1) return rmh;
	int size = rmh.mapContents->mapArray.size();
	for (int i = 0; i < size; i++)
	{
		v = rmh.mapContents->mapArray[i];
		answer.AddToEnd(v * c);
		// Formerly, answer.mapContents->mapArray.push_back(v * c);
	}
	answer.Canonicalize();
	return answer;
}

#endif
