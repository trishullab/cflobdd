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
#include <vector>
#include "list_T.h"
#include "list_TPtr.h"
#include "hashset.h"
#include "reduction_map.h"
#include "intpair.h"
#include <complex>
//#include <boost/multiprecision/cpp_int.hpp>
//#include "hash_functions.h"

template <typename T> class ReturnMapHandle;
template <typename T> class ReturnMapBody;

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
  unsigned int Hash(unsigned int modsize);
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
  static Hashset<ReturnMapBody<T>> *canonicalReturnMapBodySet;
  std::ostream& print(std::ostream & out = std::cout) const;
};

template <typename T>
std::ostream& operator<< (std::ostream & out, const ReturnMapHandle<T> &r);

//***************************************************************
// ReturnMapBody
//***************************************************************

template <typename T>
class ReturnMapBody {

  friend void ReturnMapHandle<T>::Canonicalize();
  friend unsigned int ReturnMapHandle<T>::Hash(unsigned int modsize);

 public:
  ReturnMapBody();    // Constructor
  ReturnMapBody(unsigned int capacity);    // Constructor
  //~ReturnMapBody();
  void IncrRef();
  void DecrRef();
  unsigned int Hash(unsigned int modsize);
  void setHashCheck();
  unsigned int refCount;         // reference-count value
  std::vector<T> mapArray;
  bool operator==(const ReturnMapBody &o) const;
  T& operator[](unsigned int i);                       // Overloaded []
  unsigned int hashCheck;

 protected:
  bool isCanonical;              // Is this ReturnMapBody in *canonicalReturnMapBodySet?

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
      ReturnMapHandle<T>::canonicalReturnMapBodySet->DeleteEq(this);
    }
    delete this;
  }
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

// Initializations of static members ---------------------------------
template <typename T> Hashset<ReturnMapBody<T>> *ReturnMapHandle<T>::canonicalReturnMapBodySet = new Hashset<ReturnMapBody<T>>(HASHSET_NUM_BUCKETS);

// Default constructor
template <typename T>
ReturnMapHandle<T>::ReturnMapHandle()
  :  mapContents(new ReturnMapBody<T>)
{
  mapContents->IncrRef();
}

template <typename T>
ReturnMapHandle<T>::ReturnMapHandle(unsigned int capacity)
	: mapContents(new ReturnMapBody<T>(capacity))
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
unsigned int ReturnMapHandle<T>::Hash(unsigned int modsize)
{
	if (!(mapContents->isCanonical)) {
		std::cout << "Hash of a non-canonical ReturnMapHandle occurred" << std::endl;
		abort();
		this->Canonicalize();
	}
	assert(mapContents->isCanonical);
	return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
}

template <typename T>
unsigned int ReturnMapHandle<T>::Size()
{
	try{
		return mapContents->mapArray.size();
	}
	catch (std::exception e){
		std::cout << e.what() << std::endl;
		std::cout << "Size error" << std::endl;
		throw e;
	}
}

template <typename T>
void ReturnMapHandle<T>::AddToEnd(T y)
{
	try{
		assert(mapContents->refCount <= 1);
		mapContents->mapArray.push_back(y);
	}
	catch (std::exception e){
		std::cout << e.what() << std::endl;
		//std::cout << mapContents->refCount << " " << y << " " << Size() << std::endl;
		std::cout << mapContents->refCount << " " << Size() << std::endl;
		std::cout << "ReturnMapHandle" << std::endl;
		//std::cout << y << std::endl;
		throw e;
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
		mapContents->setHashCheck();

		if (!mapContents->isCanonical) {
			unsigned int hash = canonicalReturnMapBodySet->GetHash(mapContents);
			answerContents = canonicalReturnMapBodySet->Lookup(mapContents, hash);
			if (answerContents == NULL) {
				canonicalReturnMapBodySet->Insert(mapContents, hash);
				mapContents->isCanonical = true;
			}
			else {
				answerContents->IncrRef();
				mapContents->DecrRef();
				mapContents = answerContents;
			}
		}
	}
	catch (std::exception e){
		std::cout << e.what() << std::endl;
		std::cout << "Canonicalize" << std::endl;
		throw e;
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
	catch (std::exception e){
		std::cout << e.what() << std::endl;
		std::cout << "InducedReductionAndReturnMap" << std::endl;
		throw e;
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
