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
#include <unordered_map>
#include "cflobdd_node.h"
#include "list_T.h"
#include "list_TPtr.h"
#include "intpair.h"
#include "inttriple.h"
#include "cross_product.h"

using namespace CFL_OBDD;

// ********************************************************************
// 2-Way Cross Product
// ********************************************************************

//***************************************************************
// PairProductMapBody
//***************************************************************

// Initializations of static members ---------------------------------
Hashset<PairProductMapBody> *PairProductMapBody::canonicalPairProductMapBodySet = new Hashset<PairProductMapBody>(HASHSET_NUM_BUCKETS);

// Constructor
PairProductMapBody::PairProductMapBody()
  : refCount(0), isCanonical(false)
{
}

void PairProductMapBody::IncrRef()
{
  refCount++;    // Warning: Saturation not checked
}

void PairProductMapBody::DecrRef()
{
  if (--refCount == 0) {    // Warning: Saturation not checked
    if (isCanonical) {
      PairProductMapBody::canonicalPairProductMapBodySet->DeleteEq(this);
    }
    delete this;
  }
}

unsigned int PairProductMapBody::Hash(unsigned int modsize)
{
  unsigned int hvalue = 0;

  /*
  PairProductMapBodyIterator mi(*this);

  mi.Reset();
  while (!mi.AtEnd()) {
    hvalue = (997 * (997 * hvalue + (unsigned int)mi.Current().First())
                                  + (unsigned int)mi.Current().Second() ) % modsize;
    mi.Next();
  }
  */

  for (unsigned int i = 0; i < mapArray.size(); i++){
	  hvalue = (997*hvalue + (unsigned int)97*mapArray[i].First() + (unsigned int)mapArray[i].Second()) % modsize;
  }

  return hvalue;
}

void PairProductMapBody::setHashCheck()
{
	unsigned int hvalue = 0;

	for (auto &i : mapArray) {
		hvalue = (117 * (hvalue + 1) + (int)(97 * i.First()) + i.Second());
	}
	hashCheck = hvalue;
}

void PairProductMapBody::AddToEnd(const intpair& y)
{
	mapArray.push_back(y);
}

bool PairProductMapBody::operator==(const PairProductMapBody &o) const
{
	if (mapArray.size() != o.mapArray.size())
		return false;

	for (unsigned int i = 0; i < mapArray.size(); i++){
		if (mapArray[i] != o.mapArray[i])
			return false;
	}
	return true;
}
intpair& PairProductMapBody::operator[](unsigned int i){                       // Overloaded []
	return mapArray[i];
}

unsigned int PairProductMapBody::Size(){
	return (unsigned int)mapArray.size();
}


namespace CFL_OBDD {
std::ostream& operator<< (std::ostream & out, const PairProductMapBody &r)
{
  //out << (List<int>&)r;
	for (unsigned int i = 0; i < r.mapArray.size(); i++)
	{
		out << r.mapArray[i] << " ";
	}
  return(out);
}
}
//***************************************************************
// PairProductMapHandle
//***************************************************************


// Default constructor
PairProductMapHandle::PairProductMapHandle()
  :  mapContents(new PairProductMapBody)
{
  mapContents->IncrRef();
}

// Destructor
PairProductMapHandle::~PairProductMapHandle()
{
  mapContents->DecrRef();
}

// Copy constructor
PairProductMapHandle::PairProductMapHandle(const PairProductMapHandle &r)
  :  mapContents(r.mapContents)
{
  mapContents->IncrRef();
}

// Overloaded assignment
PairProductMapHandle& PairProductMapHandle::operator= (const PairProductMapHandle &r)
{
  if (this != &r)      // don't assign to self!
  {
    PairProductMapBody *temp = mapContents;
    mapContents = r.mapContents;
    mapContents->IncrRef();
    temp->DecrRef();
  }
  return *this;        
}

// Overloaded !=
bool PairProductMapHandle::operator!=(const PairProductMapHandle &r)
{
  return (mapContents != r.mapContents);
}

// Overloaded ==
bool PairProductMapHandle::operator==(const PairProductMapHandle &r)
{
  return (mapContents == r.mapContents);
}

std::ostream& operator<< (std::ostream & out, const PairProductMapHandle &r)
{
  out << "[" << *r.mapContents << "]";
  return(out);
}

unsigned int PairProductMapHandle::Hash(unsigned int modsize)
{
  return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
}

unsigned int PairProductMapHandle::Size()
{
  return mapContents->Size();
}

intpair& PairProductMapHandle::operator[](unsigned int i)
{
	return mapContents->mapArray[i];
}

void PairProductMapHandle::AddToEnd(const intpair& p)
{
  assert(mapContents->refCount <= 1);
  mapContents->AddToEnd(p);
}

bool PairProductMapHandle::Member(intpair& p)
{
	/*
  PairProductMapBodyIterator mi(*mapContents);

  mi.Reset();
  while (!mi.AtEnd()) {
    if (mi.Current() == p) {
      return true;
    }
    mi.Next();
  }
  return false;
  */
	for (auto& i : mapContents->mapArray){
		if (i == p)
			return true;
	}
	return false;
}

int PairProductMapHandle::Lookup(intpair& p)
{
	/*  
	int index = 0;
  PairProductMapBodyIterator mi(*mapContents);

  mi.Reset();
  while (!mi.AtEnd()) {
    if (mi.Current() == p) {
      return index;
    }
    mi.Next();
    index++;
  }
  //std::cerr << "Failure in PairProductMapHandle::Lookup: " << p << " not found" << std::endl;
  std::cout << "Failure in PairProductMapHandle::Lookup: " << p << " not found" << std::endl;
  */
	for (unsigned int i = 0; i < mapContents->mapArray.size(); i++){
		if (mapContents->mapArray[i] == p)
			return i;
	}
  return -1;
}

void PairProductMapHandle::Canonicalize()
{
  PairProductMapBody *answerContents;
  unsigned int hash = PairProductMapBody::canonicalPairProductMapBodySet->GetHash(mapContents);
  answerContents = PairProductMapBody::canonicalPairProductMapBodySet->Lookup(mapContents, hash);
  if (answerContents == NULL) {
    PairProductMapBody::canonicalPairProductMapBodySet->Insert(mapContents, hash);
    mapContents->isCanonical = true;
  }
  else {
    answerContents->IncrRef();
    mapContents->DecrRef();
    mapContents = answerContents;
  }
}

// Create map with reversed entries
PairProductMapHandle PairProductMapHandle::Flip()
{
  //PairProductMapBodyIterator mi(*mapContents);
  PairProductMapHandle answer;
  /*intpair p;

  mi.Reset();
  while (!mi.AtEnd()) {
    answer.AddToEnd(intpair(mi.Current().Second(), mi.Current().First()));
    mi.Next();
  }
  */

  for (auto& i : mapContents->mapArray){
	  intpair p(i.Second(), i.First());
	  answer.AddToEnd(p);
  }
  return answer;
}

//***************************************************************
// PairProductKey
//***************************************************************

// Constructor
PairProductKey::PairProductKey(CFLOBDDNodeHandle nodeHandle1, CFLOBDDNodeHandle nodeHandle2)
  :  nodeHandle1(nodeHandle1), nodeHandle2(nodeHandle2)
{
}

// Hash
unsigned int PairProductKey::Hash(unsigned int modsize)
{
  unsigned int hvalue = 0;
  hvalue = (997 * nodeHandle1.Hash(modsize) + nodeHandle2.Hash(modsize)) % modsize;
  return hvalue;
}

// print
std::ostream& PairProductKey::print(std::ostream & out) const
{
  out << "(" << nodeHandle1 << ", " << nodeHandle2 << ")";
  return out;
}

template <typename T>
std::ostream& operator<< (std::ostream & out, const PairProductKey &p)
{
  p.print(out);
  return(out);
}

PairProductKey& PairProductKey::operator= (const PairProductKey& i)
{
  if (this != &i)      // don't assign to self!
  {
    nodeHandle1 = i.nodeHandle1;
    nodeHandle2 = i.nodeHandle2;
  }
  return *this;        
}

// Overloaded !=
bool PairProductKey::operator!=(const PairProductKey& p)
{
  return (nodeHandle1 != p.nodeHandle1) || (nodeHandle2 != p.nodeHandle2);
}

// Overloaded ==
bool PairProductKey::operator==(const PairProductKey& p)
{
  return (nodeHandle1 == p.nodeHandle1) && (nodeHandle2 == p.nodeHandle2);
}

//***************************************************************
// PairProductMemo
//***************************************************************

// Default constructor
PairProductMemo::PairProductMemo()
  :  nodeHandle(CFLOBDDNodeHandle()), pairProductMapHandle(PairProductMapHandle())
{
}

// Constructor
PairProductMemo::PairProductMemo(CFLOBDDNodeHandle nodeHandle, PairProductMapHandle pairProductMapHandle)
  :  nodeHandle(nodeHandle), pairProductMapHandle(pairProductMapHandle)
{
}

std::ostream& operator<< (std::ostream & out, const PairProductMemo &p)
{
  out << "(" << p.nodeHandle << ", " << p.pairProductMapHandle << ")";
  return(out);
}

PairProductMemo& PairProductMemo::operator= (const PairProductMemo& i)
{
  if (this != &i)      // don't assign to self!
  {
    nodeHandle = i.nodeHandle;
    pairProductMapHandle = i.pairProductMapHandle;
  }
  return *this;        
}

// Overloaded !=
bool PairProductMemo::operator!=(const PairProductMemo& p)
{
  return (nodeHandle != p.nodeHandle) || (pairProductMapHandle != p.pairProductMapHandle);
}

// Overloaded ==
bool PairProductMemo::operator==(const PairProductMemo& p)
{
  return (nodeHandle == p.nodeHandle) && (pairProductMapHandle == p.pairProductMapHandle);
}

// --------------------------------------------------------------------
// PairProduct
//
// Returns a new CFLOBDDNodeHandle, and (in pairProductMapHandle) a descriptor of the
// node's exits
// --------------------------------------------------------------------

static Hashtable<PairProductKey, PairProductMemo> *pairProductCache = NULL;

namespace CFL_OBDD {

CFLOBDDNodeHandle PairProduct(CFLOBDDInternalNode *n1,
                              CFLOBDDInternalNode *n2,
                              PairProductMapHandle &pairProductMapHandle
                             )
{
  if (n1 == CFLOBDDNodeHandle::NoDistinctionNode[n1->level].handleContents) {
    if (n2 == CFLOBDDNodeHandle::NoDistinctionNode[n2->level].handleContents) {   // ND, ND
      pairProductMapHandle.AddToEnd(intpair(0,0));
      pairProductMapHandle.Canonicalize();
      return CFLOBDDNodeHandle(n1);
    }
    else {                                                                        // ND, XX
      for (unsigned int kk = 0; kk < n2->numExits; kk++) {
        pairProductMapHandle.AddToEnd(intpair(0,kk));
      }
      pairProductMapHandle.Canonicalize();
      return CFLOBDDNodeHandle(n2);
    }
  }
  else {
    if (n2 == CFLOBDDNodeHandle::NoDistinctionNode[n2->level].handleContents) {   // XX, ND
      for (unsigned int kk = 0; kk < n1->numExits; kk++) {
        pairProductMapHandle.AddToEnd(intpair(kk,0));
      }
      pairProductMapHandle.Canonicalize();
      return CFLOBDDNodeHandle(n1);
    }
    else {                                                                        // XX, XX
      PairProductMapHandle AMap;
      CFLOBDDInternalNode *n;
      unsigned int j;
      unsigned int curExit;
      int b1, b2;
    
      n = new CFLOBDDInternalNode(n1->level);
    
      // Perform the cross product of the AConnections
      CFLOBDDNodeHandle aHandle =
                 PairProduct(*(n1->AConnection.entryPointHandle),
                             *(n2->AConnection.entryPointHandle),
                             AMap
                            );
      // Fill in n->AConnection.returnMapHandle
      // Correctness relies on AMap having no duplicates
      CFLOBDDReturnMapHandle aReturnHandle;
        for (unsigned int k = 0; k < AMap.Size(); k++) {
            aReturnHandle.AddToEnd(k);
        }
        aReturnHandle.Canonicalize();
        n->AConnection = Connection(aHandle, aReturnHandle);
      // Perform the appropriate cross products of the BConnections
         j = 0;
         curExit = 0;
         n->numBConnections = AMap.Size();
         n->BConnection = new Connection[n->numBConnections];
         //PairProductMapBodyIterator AMapIterator(*AMap.mapContents);
         //AMapIterator.Reset();
		 unsigned int Aiterator = 0;
		 std::unordered_map<intpair, unsigned int, intpair::intpair_hash> pair_to_index;
         //while (!AMapIterator.AtEnd()) {
		 while (Aiterator < AMap.Size()) {
           PairProductMapHandle BMap;
		   b1 = AMap[Aiterator].First();//AMapIterator.Current().First();
		   b2 = AMap[Aiterator].Second();//AMapIterator.Current().Second();
       CFLOBDDNodeHandle bHandle = 
                 PairProduct(*(n1->BConnection[b1].entryPointHandle),
                             *(n2->BConnection[b2].entryPointHandle),
                             BMap
                            );
        CFLOBDDReturnMapHandle bReturnHandle;
           // Fill in n->BConnection[j].returnMapHandle and add new pairs (as appropriate)
           // to pairProductMapHandle
              //PairProductMapBodyIterator BMapIterator(*BMap.mapContents);
              //BMapIterator.Reset();
		   //while (!BMapIterator.AtEnd()) {
		   unsigned int Biterator = 0;
		   while (Biterator < BMap.Size()){
                int c1, c2;
					//c1 = n1->BConnection[b1].returnMapHandle.Lookup(BMapIterator.Current().First());
					//c2 = n2->BConnection[b2].returnMapHandle.Lookup(BMapIterator.Current().Second());
				c1 = n1->BConnection[b1].returnMapHandle.Lookup(BMap[Biterator].First());
				c2 = n2->BConnection[b2].returnMapHandle.Lookup(BMap[Biterator].Second());
                // Test whether the pair (c1,c2) occurs in pairProductMapHandle
			  intpair p = intpair(c1, c2);
				auto it = pair_to_index.find(p);
                   //if (pairProductMapHandle.Member(intpair(c1,c2))) {
				if (it != pair_to_index.end()){
					bReturnHandle.AddToEnd(it->second);
                     //int index = pairProductMapHandle.Lookup(intpair(c1,c2));
                     //n->BConnection[j].returnMapHandle.AddToEnd(index);
                     // std::cout << "[PairProduct] Duplicate found: j = " << j << "; index = " << index << std::endl;
                   }
                   else {   // New pair found (i.e., new exit node found)
                     pairProductMapHandle.AddToEnd(p);
                     bReturnHandle.AddToEnd(curExit);
					 pair_to_index.emplace(p, curExit);
                     curExit++;
                   }
                //BMapIterator.Next();
				   Biterator++;
              }
           bReturnHandle.Canonicalize();
           //AMapIterator.Next();
           n->BConnection[j] = Connection(bHandle, bReturnHandle);
		   Aiterator++;
           j++;
         }
         n->numExits = curExit;
#ifdef PATH_COUNTING_ENABLED
         n->InstallPathCounts();
#endif
      pairProductMapHandle.Canonicalize();
      return CFLOBDDNodeHandle(n);
    }
  }
}

CFLOBDDNodeHandle PairProduct(CFLOBDDNodeHandle n1,
                              CFLOBDDNodeHandle n2,
                              PairProductMapHandle &pairProductMapHandle
                             )
{
  PairProductMemo cachedPairProductMemo;
  bool isCached = pairProductCache->Fetch(PairProductKey(n1,n2), cachedPairProductMemo);
  if (isCached) {
    pairProductMapHandle = cachedPairProductMemo.pairProductMapHandle;
    return cachedPairProductMemo.nodeHandle;
  }
  else if (pairProductCache->Fetch(PairProductKey(n2,n1), cachedPairProductMemo)) {
    pairProductMapHandle = cachedPairProductMemo.pairProductMapHandle.Flip();
    return cachedPairProductMemo.nodeHandle;
    }
  else {
    CFLOBDDNodeHandle answer;
  
    if (n1.handleContents->NodeKind() == CFLOBDD_INTERNAL) {
      answer = PairProduct((CFLOBDDInternalNode *)n1.handleContents,
                           (CFLOBDDInternalNode *)n2.handleContents,
                           pairProductMapHandle
                          );
    }
    else if (n1.handleContents->NodeKind() == CFLOBDD_FORK) {
      if (n2.handleContents->NodeKind() == CFLOBDD_FORK) {                 // CFLOBDD_FORK, CFLOBDD_FORK
        pairProductMapHandle.AddToEnd(intpair(0,0));
        pairProductMapHandle.AddToEnd(intpair(1,1));
        pairProductMapHandle.Canonicalize();
        answer = n1;
      }
      else { /* n2.handleContents->NodeKind() == CFLOBDD_DONTCARE */       // CFLOBDD_FORK, CFLOBDD_DONTCARE
        pairProductMapHandle.AddToEnd(intpair(0,0));
        pairProductMapHandle.AddToEnd(intpair(1,0));
        pairProductMapHandle.Canonicalize();
        answer = n1;
      }
    }
    else { /* n1.handleContents->NodeKind() == CFLOBDD_DONTCARE */
      if (n2.handleContents->NodeKind() == CFLOBDD_FORK) {                 // CFLOBDD_DONTCARE, CFLOBDD_FORK
        pairProductMapHandle.AddToEnd(intpair(0,0));
        pairProductMapHandle.AddToEnd(intpair(0,1));
        pairProductMapHandle.Canonicalize();
        answer = n2;
      }
      else { /* n2.handleContents->NodeKind() == CFLOBDD_DONTCARE */       // CFLOBDD_DONTCARE, CFLOBDD_DONTCARE
        pairProductMapHandle.AddToEnd(intpair(0,0));
        pairProductMapHandle.Canonicalize();
        answer = n1;
      }
    }
		pairProductCache->Insert(PairProductKey(n1, n2),
		PairProductMemo(answer, pairProductMapHandle));
    return answer;
  }
}


void InitPairProductCache()
{
  pairProductCache = new Hashtable<PairProductKey, PairProductMemo>(HASH_NUM_BUCKETS);
}

void DisposeOfPairProductCache()
{
	//std::cout << "PairProductCache Size: " << pairProductCache->Size() << std::endl;
	delete pairProductCache;
	pairProductCache = NULL;
}
}
// ********************************************************************
// 3-Way Cross Product
// ********************************************************************

//***************************************************************
// TripleProductMapBody
//***************************************************************

// Initializations of static members ---------------------------------
Hashset<TripleProductMapBody> *TripleProductMapBody::canonicalTripleProductMapBodySet = new Hashset<TripleProductMapBody>;

// Constructor
TripleProductMapBody::TripleProductMapBody()
  : refCount(0), isCanonical(false)
{
}

void TripleProductMapBody::IncrRef()
{
  refCount++;    // Warning: Saturation not checked
}

void TripleProductMapBody::DecrRef()
{
  if (--refCount == 0) {    // Warning: Saturation not checked
    if (isCanonical) {
      TripleProductMapBody::canonicalTripleProductMapBodySet->DeleteEq(this);
    }
    delete this;
  }
}

unsigned int TripleProductMapBody::Hash(unsigned int modsize)
{
  unsigned int hvalue = 0;
  TripleProductMapBodyIterator mi(*this);

  mi.Reset();
  while (!mi.AtEnd()) {
    hvalue = (997 * (997 * (997 * hvalue + (unsigned int)mi.Current().First())
                                         + (unsigned int)mi.Current().Second())
                                         + (unsigned int)mi.Current().Third()) % modsize;
    mi.Next();
  }
  return hvalue;
}

std::ostream& operator<< (std::ostream & out, const TripleProductMapBody &r)
{
  out << (List<int>&)r;
  return(out);
}

//***************************************************************
// TripleProductMapHandle
//***************************************************************

// Default constructor
TripleProductMapHandle::TripleProductMapHandle()
  :  mapContents(new TripleProductMapBody)
{
  mapContents->IncrRef();
}

// Destructor
TripleProductMapHandle::~TripleProductMapHandle()
{
  mapContents->DecrRef();
}

// Copy constructor
TripleProductMapHandle::TripleProductMapHandle(const TripleProductMapHandle &r)
  :  mapContents(r.mapContents)
{
  mapContents->IncrRef();
}

// Overloaded assignment
TripleProductMapHandle& TripleProductMapHandle::operator= (const TripleProductMapHandle &r)
{
  if (this != &r)      // don't assign to self!
  {
    TripleProductMapBody *temp = mapContents;
    mapContents = r.mapContents;
    mapContents->IncrRef();
    temp->DecrRef();
  }
  return *this;        
}

// Overloaded !=
bool TripleProductMapHandle::operator!=(const TripleProductMapHandle &r)
{
  return (mapContents != r.mapContents);
}

// Overloaded ==
bool TripleProductMapHandle::operator==(const TripleProductMapHandle &r)
{
  return (mapContents == r.mapContents);
}

std::ostream& operator<< (std::ostream & out, const TripleProductMapHandle &r)
{
  out << "[" << *r.mapContents << "]";
  return(out);
}

unsigned int TripleProductMapHandle::Hash(unsigned int modsize)
{
  return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
}

unsigned int TripleProductMapHandle::Size()
{
  return mapContents->Length();
}

void TripleProductMapHandle::AddToEnd(inttriple t)
{
  assert(mapContents->refCount <= 1);
  mapContents->AddToEnd(t);
}

bool TripleProductMapHandle::Member(inttriple t)
{
  TripleProductMapBodyIterator mi(*mapContents);

  mi.Reset();
  while (!mi.AtEnd()) {
    if (mi.Current() == t) {
      return true;
    }
    mi.Next();
  }
  return false;
}

int TripleProductMapHandle::Lookup(inttriple t)
{
  int index = 0;
  TripleProductMapBodyIterator mi(*mapContents);

  mi.Reset();
  while (!mi.AtEnd()) {
    if (mi.Current() == t) {
      return index;
    }
    mi.Next();
    index++;
  }
  std::cerr << "Failure in TripleProductMapHandle::Lookup: " << t << " not found" << std::endl;
  return -1;
}

void TripleProductMapHandle::Canonicalize()
{
  TripleProductMapBody *answerContents;
  unsigned int hash = TripleProductMapBody::canonicalTripleProductMapBodySet->GetHash(mapContents);
  answerContents = TripleProductMapBody::canonicalTripleProductMapBodySet->Lookup(mapContents, hash);
  if (answerContents == NULL) {
    TripleProductMapBody::canonicalTripleProductMapBodySet->Insert(mapContents, hash);
    mapContents->isCanonical = true;
  }
  else {
    answerContents->IncrRef();
    mapContents->DecrRef();
    mapContents = answerContents;
  }
}

//***************************************************************
// TripleProductKey
//***************************************************************

// Constructor
TripleProductKey::TripleProductKey(CFLOBDDNodeHandle nodeHandle1, CFLOBDDNodeHandle nodeHandle2, CFLOBDDNodeHandle nodeHandle3)
  :  nodeHandle1(nodeHandle1), nodeHandle2(nodeHandle2), nodeHandle3(nodeHandle3)
{
}

// Hash
unsigned int TripleProductKey::Hash(unsigned int modsize)
{
  unsigned int hvalue = 0;
  hvalue = (997 * (997 * nodeHandle1.Hash(modsize) + nodeHandle2.Hash(modsize)) + nodeHandle3.Hash(modsize)) % modsize;
  return hvalue;
}

// print
std::ostream& TripleProductKey::print(std::ostream & out) const
{
  out << "(" << nodeHandle1 << ", " << nodeHandle2 << ", " << nodeHandle3 << ")";
  return out;
}

std::ostream& operator<< (std::ostream & out, const TripleProductKey &p)
{
  p.print(out);
  return(out);
}

TripleProductKey& TripleProductKey::operator= (const TripleProductKey& i)
{
  if (this != &i)      // don't assign to self!
  {
    nodeHandle1 = i.nodeHandle1;
    nodeHandle2 = i.nodeHandle2;
    nodeHandle2 = i.nodeHandle3;
  }
  return *this;        
}

// Overloaded !=
bool TripleProductKey::operator!=(const TripleProductKey& p)
{
  return (nodeHandle1 != p.nodeHandle1) || (nodeHandle2 != p.nodeHandle2) || (nodeHandle3 != p.nodeHandle3);
}

// Overloaded ==
bool TripleProductKey::operator==(const TripleProductKey& p)
{
  return (nodeHandle1 == p.nodeHandle1) && (nodeHandle2 == p.nodeHandle2) && (nodeHandle3 == p.nodeHandle3);
}

//***************************************************************
// TripleProductMemo
//***************************************************************

// Default constructor
TripleProductMemo::TripleProductMemo()
  :  nodeHandle(CFLOBDDNodeHandle()), tripleProductMapHandle(TripleProductMapHandle())
{
}

// Constructor
TripleProductMemo::TripleProductMemo(CFLOBDDNodeHandle nodeHandle, TripleProductMapHandle tripleProductMapHandle)
  :  nodeHandle(nodeHandle), tripleProductMapHandle(tripleProductMapHandle)
{
}

std::ostream& operator<< (std::ostream & out, const TripleProductMemo &p)
{
  out << "(" << p.nodeHandle << ", " << p.tripleProductMapHandle << ")";
  return(out);
}

TripleProductMemo& TripleProductMemo::operator= (const TripleProductMemo& i)
{
  if (this != &i)      // don't assign to self!
  {
    nodeHandle = i.nodeHandle;
    tripleProductMapHandle = i.tripleProductMapHandle;
  }
  return *this;        
}

// Overloaded !=
bool TripleProductMemo::operator!=(const TripleProductMemo& p)
{
  return (nodeHandle != p.nodeHandle) || (tripleProductMapHandle != p.tripleProductMapHandle);
}

// Overloaded ==
bool TripleProductMemo::operator==(const TripleProductMemo& p)
{
  return (nodeHandle == p.nodeHandle) && (tripleProductMapHandle == p.tripleProductMapHandle);
}

// --------------------------------------------------------------------
// TripleProduct
//
// Returns a new CFLOBDDNodeHandle, and (in tripleProductMap) a descriptor of the
// node's exits
// --------------------------------------------------------------------

static Hashtable<TripleProductKey, TripleProductMemo> *tripleProductCache = NULL;

namespace CFL_OBDD {

CFLOBDDNodeHandle TripleProduct(CFLOBDDInternalNode *n1,
                                CFLOBDDInternalNode *n2,
                                CFLOBDDInternalNode *n3,
                                TripleProductMapHandle &tripleProductMapHandle
                               )
{
  if (n1 == CFLOBDDNodeHandle::NoDistinctionNode[n1->level].handleContents) {
    if (n2 == CFLOBDDNodeHandle::NoDistinctionNode[n2->level].handleContents) {
      if (n3 == CFLOBDDNodeHandle::NoDistinctionNode[n3->level].handleContents) {   // ND, ND, ND
        tripleProductMapHandle.AddToEnd(inttriple(0,0,0));
        tripleProductMapHandle.Canonicalize();
        return CFLOBDDNodeHandle(n1);
      }
      else {                                                   // ND, ND, XX
        for (unsigned int kk = 0; kk < n3->numExits; kk++) {
          tripleProductMapHandle.AddToEnd(inttriple(0,0,kk));
        }
        tripleProductMapHandle.Canonicalize();
        return CFLOBDDNodeHandle(n3);
      }
    }
    else {
      if (n3 == CFLOBDDNodeHandle::NoDistinctionNode[n3->level].handleContents) {   // ND, XX, ND
        for (unsigned int kk = 0; kk < n2->numExits; kk++) {
          tripleProductMapHandle.AddToEnd(inttriple(0,kk,0));
        }
        tripleProductMapHandle.Canonicalize();
        return CFLOBDDNodeHandle(n2);
      }
      else {                                                   // ND, XX, XX
        PairProductMapHandle PPMapHandle;
        CFLOBDDNodeHandle n = PairProduct(n2, n3, PPMapHandle);
        //PairProductMapBodyIterator PPMapIterator(*PPMapHandle.mapContents);
        //PPMapIterator.Reset();
        //while (!PPMapIterator.AtEnd()) {
          //int c1 = PPMapIterator.Current().First();
          //int c2 = PPMapIterator.Current().Second();
		unsigned int iterator = 0;
		while (iterator < PPMapHandle.Size()){
			int c1 = PPMapHandle[iterator].First();
			int c2 = PPMapHandle[iterator].Second();
          tripleProductMapHandle.AddToEnd(inttriple(0,c1,c2));
          //PPMapIterator.Next();
		  iterator++;
        }
        tripleProductMapHandle.Canonicalize();
        return n;
      }
    }
  }
  else {
    if (n2 == CFLOBDDNodeHandle::NoDistinctionNode[n2->level].handleContents) {
      if (n3 == CFLOBDDNodeHandle::NoDistinctionNode[n3->level].handleContents) {   // XX, ND, ND
        for (unsigned int kk = 0; kk < n1->numExits; kk++) {
          tripleProductMapHandle.AddToEnd(inttriple(kk,0,0));
        }
        tripleProductMapHandle.Canonicalize();
        return CFLOBDDNodeHandle(n1);
      }
      else {                                                   // XX, ND, XX
        PairProductMapHandle PPMapHandle;
        CFLOBDDNodeHandle n = PairProduct(n1, n3, PPMapHandle);
        //PairProductMapBodyIterator PPMapIterator(*PPMapHandle.mapContents);
        //PPMapIterator.Reset();
        //while (!PPMapIterator.AtEnd()) {
          //int c1 = PPMapIterator.Current().First();
          //int c2 = PPMapIterator.Current().Second();
		unsigned int iterator = 0;
		while (iterator < PPMapHandle.Size()){
			int c1 = PPMapHandle[iterator].First();
			int c2 = PPMapHandle[iterator].Second();
          tripleProductMapHandle.AddToEnd(inttriple(c1,0,c2));
          //PPMapIterator.Next();
		  iterator++;
        }
        tripleProductMapHandle.Canonicalize();
        return n;
      }
    }
    else {
      if (n3 == CFLOBDDNodeHandle::NoDistinctionNode[n3->level].handleContents) {   // XX, XX, ND
        PairProductMapHandle PPMapHandle;
        CFLOBDDNodeHandle n = PairProduct(n1, n2, PPMapHandle);
        //PairProductMapBodyIterator PPMapIterator(*PPMapHandle.mapContents);
        //PPMapIterator.Reset();
        //while (!PPMapIterator.AtEnd()) {
        //  int c1 = PPMapIterator.Current().First();
        //  int c2 = PPMapIterator.Current().Second();
		unsigned int iterator = 0;
		while (iterator < PPMapHandle.Size()){
			int c1 = PPMapHandle[iterator].First();
			int c2 = PPMapHandle[iterator].Second();
          tripleProductMapHandle.AddToEnd(inttriple(c1,c2,0));
        //  PPMapIterator.Next();
		  iterator++;
        }
        tripleProductMapHandle.Canonicalize();
        return n;
      }
      else {                                                   // XX, XX, XX
        TripleProductMapHandle AMap;
        unsigned int j;
        unsigned int curExit;
        int b1, b2, b3;

        CFLOBDDInternalNode *n = new CFLOBDDInternalNode(n1->level);
      
        // Perform the cross product of the AConnections
           *(n->AConnection.entryPointHandle) =
                   TripleProduct(*(n1->AConnection.entryPointHandle),
                                 *(n2->AConnection.entryPointHandle),
                                 *(n3->AConnection.entryPointHandle),
                                 AMap
                                );
           // Fill in n->AConnection.returnMapHandle
           // Correctness relies on AMap having no duplicates
              for (unsigned int k = 0; k < AMap.Size(); k++) {
                 n->AConnection.returnMapHandle.AddToEnd(k);
              }
              n->AConnection.returnMapHandle.Canonicalize();
      
        // Perform the appropriate cross products of the BConnections
           j = 0;
           curExit = 0;
           n->numBConnections = AMap.Size();
           n->BConnection = new Connection[n->numBConnections];
           TripleProductMapBodyIterator AMapIterator(*AMap.mapContents);
           AMapIterator.Reset();
           while (!AMapIterator.AtEnd()) {
             TripleProductMapHandle BMap;
             b1 = AMapIterator.Current().First();
             b2 = AMapIterator.Current().Second();
             b3 = AMapIterator.Current().Third();
             *(n->BConnection[j].entryPointHandle) =
                   TripleProduct(*(n1->BConnection[b1].entryPointHandle),
                                 *(n2->BConnection[b2].entryPointHandle),
                                 *(n3->BConnection[b3].entryPointHandle),
                                 BMap
                                );
      
             // Fill in n->BConnection[j].returnMapHandle and add new triples (as appropriate)
             // to tripleProductMap
                TripleProductMapBodyIterator BMapIterator(*BMap.mapContents);
                BMapIterator.Reset();
                while (!BMapIterator.AtEnd()) {
                  int c1, c2, c3;
                  c1 = n1->BConnection[b1].returnMapHandle.Lookup(BMapIterator.Current().First());
                  c2 = n2->BConnection[b2].returnMapHandle.Lookup(BMapIterator.Current().Second());
                  c3 = n3->BConnection[b3].returnMapHandle.Lookup(BMapIterator.Current().Third());
                  // Test whether the triple (c1,c2,c3) occurs in tripleProductMapHandle
                     if (tripleProductMapHandle.Member(inttriple(c1,c2,c3))) {
                       int index = tripleProductMapHandle.Lookup(inttriple(c1,c2,c3));
                       n->BConnection[j].returnMapHandle.AddToEnd(index);
                         // std::cout << "[TripleProduct] Duplicate found: j = " << j << "; index = " << index << std::endl;
                     }
                     else {   // New triple found (i.e., new exit node found)
                       tripleProductMapHandle.AddToEnd(inttriple(c1,c2,c3));
                       n->BConnection[j].returnMapHandle.AddToEnd(curExit);
                       curExit++;
                     }
                  BMapIterator.Next();
                }
             n->BConnection[j].returnMapHandle.Canonicalize();
             AMapIterator.Next();
             j++;
           }
           n->numExits = curExit;
#ifdef PATH_COUNTING_ENABLED
           n->InstallPathCounts();
#endif
        tripleProductMapHandle.Canonicalize();
        return CFLOBDDNodeHandle(n);
      }
    }
  }
}

CFLOBDDNodeHandle TripleProduct(CFLOBDDNodeHandle n1,
                                CFLOBDDNodeHandle n2,
                                CFLOBDDNodeHandle n3,
                                TripleProductMapHandle &tripleProductMapHandle
                               )
{
  TripleProductMemo cachedTripleProductMemo;
  bool isCached = tripleProductCache->Fetch(TripleProductKey(n1,n2,n3), cachedTripleProductMemo);
  if (isCached) {
    tripleProductMapHandle = cachedTripleProductMemo.tripleProductMapHandle;
    return cachedTripleProductMemo.nodeHandle;
  }
  else {
    CFLOBDDNodeHandle answer;
  
    if (n1.handleContents->NodeKind() == CFLOBDD_INTERNAL) {
      answer = TripleProduct((CFLOBDDInternalNode *)n1.handleContents,
                             (CFLOBDDInternalNode *)n2.handleContents,
                             (CFLOBDDInternalNode *)n3.handleContents,
                             tripleProductMapHandle
                            );
    }
    else if (n1.handleContents->NodeKind() == CFLOBDD_FORK) {
      if (n2.handleContents->NodeKind() == CFLOBDD_FORK) {
        if (n3.handleContents->NodeKind() == CFLOBDD_FORK) {                      // CFLOBDD_FORK, CFLOBDD_FORK, CFLOBDD_FORK
          tripleProductMapHandle.AddToEnd(inttriple(0,0,0));
          tripleProductMapHandle.AddToEnd(inttriple(1,1,1));
          tripleProductMapHandle.Canonicalize();
          answer = n1;
        }
        else { /* n3.handleContents->NodeKind() == CFLOBDD_DONTCARE */            // CFLOBDD_FORK, CFLOBDD_FORK, CFLOBDD_DONTCARE
          tripleProductMapHandle.AddToEnd(inttriple(0,0,0));
          tripleProductMapHandle.AddToEnd(inttriple(1,1,0));
          tripleProductMapHandle.Canonicalize();
          answer = n1;
        }
      }
      else { /* n2.handleContents->NodeKind() == CFLOBDD_DONTCARE */
        if (n3.handleContents->NodeKind() == CFLOBDD_FORK) {                      // CFLOBDD_FORK, CFLOBDD_DONTCARE, CFLOBDD_FORK
          tripleProductMapHandle.AddToEnd(inttriple(0,0,0));
          tripleProductMapHandle.AddToEnd(inttriple(1,0,1));
          tripleProductMapHandle.Canonicalize();
          answer = n1;
        }
        else { /* n3.handleContents->NodeKind() == CFLOBDD_DONTCARE */            // CFLOBDD_FORK, CFLOBDD_DONTCARE, CFLOBDD_DONTCARE
          tripleProductMapHandle.AddToEnd(inttriple(0,0,0));
          tripleProductMapHandle.AddToEnd(inttriple(1,0,0));
          tripleProductMapHandle.Canonicalize();
          answer = n1;
        }
      }
    }
    else { /* n1.handleContents->NodeKind() == CFLOBDD_DONTCARE */
      if (n2.handleContents->NodeKind() == CFLOBDD_FORK) {
        if (n3.handleContents->NodeKind() == CFLOBDD_FORK) {                      // CFLOBDD_DONTCARE, CFLOBDD_FORK, CFLOBDD_FORK
          tripleProductMapHandle.AddToEnd(inttriple(0,0,0));
          tripleProductMapHandle.AddToEnd(inttriple(0,1,1));
          tripleProductMapHandle.Canonicalize();
          answer = n2;
        }
        else { /* n3.handleContents->NodeKind() == CFLOBDD_DONTCARE */            // CFLOBDD_DONTCARE, CFLOBDD_FORK, CFLOBDD_DONTCARE
          tripleProductMapHandle.AddToEnd(inttriple(0,0,0));
          tripleProductMapHandle.AddToEnd(inttriple(0,1,0));
          tripleProductMapHandle.Canonicalize();
          answer = n2;
        }
      }
      else { /* n2.handleContents->NodeKind() == CFLOBDD_DONTCARE */
        if (n3.handleContents->NodeKind() == CFLOBDD_FORK) {                      // CFLOBDD_DONTCARE, CFLOBDD_DONTCARE, CFLOBDD_FORK
          tripleProductMapHandle.AddToEnd(inttriple(0,0,0));
          tripleProductMapHandle.AddToEnd(inttriple(0,0,1));
          tripleProductMapHandle.Canonicalize();
          answer = n3;
        }
        else { /* n3.handleContents->NodeKind() == CFLOBDD_DONTCARE */            // CFLOBDD_DONTCARE, CFLOBDD_DONTCARE, CFLOBDD_DONTCARE
          tripleProductMapHandle.AddToEnd(inttriple(0,0,0));
          answer = n1;
        }
      }
    }
    tripleProductCache->Insert(TripleProductKey(n1,n2,n3),
                               TripleProductMemo(answer,tripleProductMapHandle));
    return answer;
  }
}

void InitTripleProductCache()
{
  tripleProductCache = new Hashtable<TripleProductKey, TripleProductMemo>(40000);
}

void DisposeOfTripleProductCache()
{
	delete tripleProductCache;
	tripleProductCache = NULL;
}
}
