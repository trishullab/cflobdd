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
#include "return_map_T.h"
#include "reduction_map.h"
#include "list_T.h"
#include "list_TPtr.h"
#include "intpair.h"
#include "hashset.h"

//***************************************************************
// ReductionMapBodyIterator
//***************************************************************

typedef ListIterator<int> ReductionMapBodyIterator;

//***************************************************************
// ReductionMapBody
//***************************************************************

// Constructor
ReductionMapBody::ReductionMapBody()
  : refCount(0), isIdentityMap(true), isCanonical(false)
{
}

ReductionMapBody::ReductionMapBody(unsigned int capacity)
	: refCount(0), isIdentityMap(true), isCanonical(false)
{
	mapArray.reserve(capacity);
}

void ReductionMapBody::IncrRef()
{
  refCount++;    // Warning: Saturation not checked
}

void ReductionMapBody::DecrRef()
{
  if (--refCount == 0) {    // Warning: Saturation not checked
    if (isCanonical) {
      ReductionMapHandle::canonicalReductionMapBodySet->DeleteEq(this);
		//canonicalReductionMapBodySet.erase(this);
    }
    delete this;
  }
}

unsigned int ReductionMapBody::Hash(unsigned int modsize)
{
  unsigned int hvalue = 0;
  /*ReductionMapBodyIterator mi(*this);

  mi.Reset();
  while (!mi.AtEnd()) {
    hvalue = (hvalue + (unsigned int)mi.Current()) % modsize;
    mi.Next();
  }*/
  for (unsigned int i = 0; i < mapArray.size(); i++){
	  hvalue = (117 * (hvalue + 1) + (unsigned int)mapArray[i]) % modsize;
  }

  return hvalue;
}

void ReductionMapBody::setHashCheck()
{
	unsigned int hvalue = 0;

	for (auto &i : mapArray) {
		hvalue = (117 * (hvalue + 1) + (int)(i));
	}
	hashCheck = hvalue;
}

void ReductionMapBody::AddToEnd(int y)
{
  /*isIdentityMap = isIdentityMap && (y == Length());
  List<int>::AddToEnd(y);*/
	isIdentityMap = isIdentityMap && (y == mapArray.size());
	mapArray.push_back(y);
}

bool ReductionMapBody::operator==(const ReductionMapBody &o) const
{
	if (hashCheck != o.hashCheck)
		return false;

	if (mapArray.size() != o.mapArray.size())
		return false;

	for (unsigned int i = 0; i < mapArray.size(); i++){
		if (mapArray[i] != o.mapArray[i])
			return false;
	}
	return true;
}
int& ReductionMapBody::operator[](unsigned int i){                       // Overloaded []
	return mapArray[i];
}

unsigned int ReductionMapBody::Size(){
	return (unsigned int)mapArray.size();
}

std::ostream& operator<< (std::ostream & out, const ReductionMapBody &r)
{
  //out << (List<int>&)r;
	for (unsigned int i = 0; i < r.mapArray.size(); i++)
	{
		out << r.mapArray[i] << " ";
	}
  return(out);
}

//***************************************************************
// ReductionMapHandle
//***************************************************************

// Initializations of static members ---------------------------------
Hashset<ReductionMapBody> *ReductionMapHandle::canonicalReductionMapBodySet = new Hashset<ReductionMapBody>(HASHSET_NUM_BUCKETS);

// Default constructor
ReductionMapHandle::ReductionMapHandle()
  :  mapContents(new ReductionMapBody)
{
  mapContents->IncrRef();
}

// Destructor
ReductionMapHandle::~ReductionMapHandle()
{
  mapContents->DecrRef();
}

// Copy constructor
ReductionMapHandle::ReductionMapHandle(const ReductionMapHandle &r)
  :  mapContents(r.mapContents)
{
  mapContents->IncrRef();
}

ReductionMapHandle::ReductionMapHandle(unsigned int capacity)
	: mapContents(new ReductionMapBody(capacity))
{
	mapContents->IncrRef();
}

// Overloaded assignment
ReductionMapHandle& ReductionMapHandle::operator= (const ReductionMapHandle &r)
{
  if (this != &r)      // don't assign to self!
  {
    ReductionMapBody *temp = mapContents;
    mapContents = r.mapContents;
    mapContents->IncrRef();
    temp->DecrRef();
  }
  return *this;        
}

// Overloaded !=
bool ReductionMapHandle::operator!=(const ReductionMapHandle &r)
{
  return (mapContents != r.mapContents);
}

// Overloaded ==
bool ReductionMapHandle::operator==(const ReductionMapHandle &r)
{
  return (mapContents == r.mapContents);
}

// print
std::ostream& ReductionMapHandle::print(std::ostream & out) const
{
  out << *mapContents << std::endl;
  return out;
}

std::ostream& operator<< (std::ostream & out, const ReductionMapHandle &r)
{
  r.print(out);
  return(out);
}

unsigned int ReductionMapHandle::Hash(unsigned int modsize)
{
  return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
}

unsigned int ReductionMapHandle::Size()
{
  //return mapContents->Length();
	return mapContents->Size();
}

void ReductionMapHandle::AddToEnd(int y)
{
  assert(mapContents->refCount <= 1);
  mapContents->AddToEnd(y);
}

/*
intpair ReductionMapHandle::Lookup(intpair x)
{
	ReductionMapBodyIterator mi(*mapContents);
	int xx = 0;
	int x1 = x.First();
	int x2 = x.Second();
	int temp;
	mi.Reset();
	bool done = false;
	while (!mi.AtEnd()) {
		if (xx == x1) {
			if (done) {
				return intpair(mi.Current(),temp);
			}
			else {
				done = true;
				temp = mi.Current();
			}
		} 
		if (xx == x2) {
			if (done) {
				return intpair(temp,mi.Current());
			} else {
				done = true;
				temp = mi.Current();
			}
		}
		xx++;
		mi.Next();
	}
	std::cerr << "Failure in ReductionMapHandle::Lookup: " << x << " not found" << std::endl;
	return intpair(-1,-1);
}

int ReductionMapHandle::Lookup(int x)
{
  ReductionMapBodyIterator mi(*mapContents);

  int xx = 0;
  mi.Reset();
  while (!mi.AtEnd()) {
    if (xx == x) {
      return mi.Current();
    }
    xx++;
    mi.Next();
  }
  std::cerr << "Failure in ReductionMapHandle::Lookup: " << x << " not found" << std::endl;
  return -1;
}

int ReductionMapHandle::LookupInv(int y)
{
  ReductionMapBodyIterator mi(*mapContents);

  int x = 0;
  mi.Reset();
  while (!mi.AtEnd()) {
    if (mi.Current() == y) {
      return x;
    }
    x++;
    mi.Next();
  }
  return -1;
}
*/

intpair ReductionMapHandle::Lookup(intpair& x)
{
	if ((unsigned int)x.First() < Size() && (unsigned int)x.Second() < Size())
		return intpair(mapContents->mapArray[x.First()], mapContents->mapArray[x.Second()]);
	return intpair(-1, -1);
}

int ReductionMapHandle::Lookup(int x)
{
	if ((unsigned int)x < Size())
		return mapContents->mapArray[x];
	return -1;
}

int ReductionMapHandle::LookupInv(int y)
{
	for (unsigned int i = 0; i < Size(); i++){
		if (mapContents->mapArray[i] == y)
			return i;
	}
	return -1;
}


void ReductionMapHandle::Canonicalize()
{
  ReductionMapBody *answerContents;
  mapContents->setHashCheck();

  if (!mapContents->isCanonical) {
	unsigned int hash = canonicalReductionMapBodySet->GetHash(mapContents);
    answerContents = canonicalReductionMapBodySet->Lookup(mapContents, hash);
    if (answerContents == NULL) {
      canonicalReductionMapBodySet->Insert(mapContents, hash);
      mapContents->isCanonical = true;
    }
    else {
      answerContents->IncrRef();
      mapContents->DecrRef();
      mapContents = answerContents;
    }
	/*
	  auto it = canonicalReductionMapBodySet.find(mapContents);
	  if (it == canonicalReductionMapBodySet.end()) {
		  mapContents->isCanonical = true;
		  canonicalReductionMapBodySet.insert(mapContents);
	  }
	  else {
		  answerContents = *it;
		answerContents->IncrRef();
		mapContents->DecrRef();
		mapContents = answerContents;
	  }
	  */
  }
}

std::size_t hash_value(const ReductionMapHandle& val)
{
	return val.mapContents->hashCheck;
}



