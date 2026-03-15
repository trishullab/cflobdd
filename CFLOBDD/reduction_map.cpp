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
#include <cstdint>

//***************************************************************
// ReductionMapBodyIterator
//***************************************************************

typedef ListIterator<int> ReductionMapBodyIterator;

//***************************************************************
// ReductionMapBody
//***************************************************************

// Constructor
ReductionMapBody::ReductionMapBody()
  : refCount(0), isIdentityMap(true), isCanonical(false), hashCheck(0)
{
}

ReductionMapBody::ReductionMapBody(unsigned int capacity)
	: refCount(0), isIdentityMap(true), isCanonical(false), hashCheck(0)
{
	mapArray.reserve(capacity);
}

ReductionMapBody* ReductionMapBody::Create()
{
    auto& fl = getFreeList();
    if (!fl.empty()) {
        ReductionMapBody* p = fl.back();
        fl.pop_back();
        return p;
    }
    return new ReductionMapBody();
}

ReductionMapBody* ReductionMapBody::Create(unsigned int capacity)
{
    auto& fl = getFreeList();
    if (!fl.empty()) {
        ReductionMapBody* p = fl.back();
        fl.pop_back();
        if (p->mapArray.capacity() < capacity)
            p->mapArray.reserve(capacity);
        return p;
    }
    return new ReductionMapBody(capacity);
}

void ReductionMapBody::IncrRef()
{
  refCount++;    // Warning: Saturation not checked
}

void ReductionMapBody::DecrRef()
{
  if (--refCount == 0) {    // Warning: Saturation not checked
    if (isCanonical) {
      ReductionMapHandle::canonicalReductionMapBodySet->erase(this);
    }
    mapArray.clear();
    refCount = 0;
    isCanonical = false;
    hashCheck = 0;
    isIdentityMap = true;
    getFreeList().push_back(this);
  }
}

// Murmur3 finalizer — ensures full avalanche (each output bit depends on all input bits)
static inline size_t fmix64(size_t h) {
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
}

size_t ReductionMapBody::Hash()
{
  return fmix64(hashCheck);
}

void ReductionMapBody::setHashCheck()
{
  unsigned int hvalue = 0;
  const auto* data = mapArray.data();
  unsigned int sz = mapArray.size();
  for (unsigned int i = 0; i < sz; i++) {
      hvalue = (131 * (hvalue + 1) + data[i]);
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
	unsigned int sz = mapArray.size();
	if (sz != o.mapArray.size())
		return false;
	const auto* d1 = mapArray.data();
	const auto* d2 = o.mapArray.data();
	for (unsigned int i = 0; i < sz; i++){
		if (d1[i] != d2[i])
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
ReductionMapHandle::CanonicalReductionMapBodySet *ReductionMapHandle::initCanonicalSet()
{
    auto *s = new CanonicalReductionMapBodySet(REDUCTION_MAP_NUM_BUCKETS);
    s->max_load_factor(0.8f);
    return s;
}

ReductionMapHandle::CanonicalReductionMapBodySet
    *ReductionMapHandle::canonicalReductionMapBodySet = ReductionMapHandle::initCanonicalSet();

// Default constructor
ReductionMapHandle::ReductionMapHandle()
  :  mapContents(ReductionMapBody::Create())
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
	: mapContents(ReductionMapBody::Create(capacity))
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

size_t ReductionMapHandle::Hash()
{
  return reinterpret_cast<uintptr_t>(mapContents) >> PTR_ALIGN_SHIFT;
}

void ReductionMapHandle::AddToEnd(int y)
{
  assert(mapContents->refCount <= 1);
  mapContents->AddToEnd(y);
}

int ReductionMapHandle::LookupInv(int y)
{
	const auto* data = mapContents->mapArray.data();
	unsigned int sz = Size();
	for (unsigned int i = 0; i < sz; i++){
		if (data[i] == y)
			return i;
	}
	return -1;
}


void ReductionMapHandle::Canonicalize()
{
  ReductionMapBody *answerContents;

  if (!mapContents->isCanonical) {
    mapContents->setHashCheck();
    auto it = canonicalReductionMapBodySet->find(mapContents);
    if (it == canonicalReductionMapBodySet->end()) {
      canonicalReductionMapBodySet->insert(mapContents);
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

