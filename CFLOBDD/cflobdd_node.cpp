//
//    Copyright (c) 1999, 2017 Thomas W. Reps
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
#include <algorithm>
#include <cstdarg>
#include <unordered_set>
#include <map>
//#include <mpirxx.h>

#include "ntz_T.h"
#include "cflobdd_node.h"
#include "connectionT.h"
#include "list_T.h"
#include "list_TPtr.h"
#include "intpair.h"
#include "conscell.h"
#include "assignment.h"
#include "bool_op.h"
#include "return_map_T.h"
#include "reduction_map.h"
#include "cross_product.h"
#include "traverse_state_cfl.h"
#include "hash.h"
#include "hashset.h"

using namespace CFL_OBDD;

//********************************************************************
// CFLOBDDNodeHandle
//
// Contains a canonical CFLOBDDNode*
//********************************************************************

// Initializations of static members ---------------------------------

Hashset<CFLOBDDNode> *CFLOBDDNodeHandle::canonicalNodeTable = new Hashset<CFLOBDDNode>(HASHSET_NUM_BUCKETS);
CFLOBDDNodeHandle *CFLOBDDNodeHandle::NoDistinctionNode = NULL;
CFLOBDDNodeHandle CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
CFLOBDDNodeHandle CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle;
CFLOBDDNodeHandle *CFLOBDDNodeHandle::AdditionInterleavedNoCarryNode = NULL;
CFLOBDDNodeHandle *CFLOBDDNodeHandle::AdditionInterleavedCarryNode = NULL;
std::vector<ReturnMapHandle<int>> commonly_used_return_maps;// m0, m1, m01, m10

void InitReturnMapHandles(){
	ReturnMapHandle<int> m0, m1, m01, m10;
	m0.AddToEnd(0);
	m0.Canonicalize();
	m1.AddToEnd(1);
	m1.Canonicalize();
	m01.AddToEnd(0);
	m01.AddToEnd(1);
	m01.Canonicalize();
	m10.AddToEnd(1);
	m10.AddToEnd(0);
	m10.Canonicalize();
	commonly_used_return_maps.push_back(m0);
	commonly_used_return_maps.push_back(m1);
	commonly_used_return_maps.push_back(m01);
	commonly_used_return_maps.push_back(m10);
}

// InitNoDistinctionTable
//
// Create and record the "NoDistinctionNodes" (one per level).
// As a side effect, these nodes are entered in the canonicalNodeTable.
//
void CFLOBDDNodeHandle::InitNoDistinctionTable()
{
	// TODO: Place this in a better place
  InitReturnMapHandles();
  NoDistinctionNode = new CFLOBDDNodeHandle [CFLOBDDMaxLevel+1];

  for (unsigned int i = 0; i <= CFLOBDDMaxLevel; i++) {
    if (i == 0) {
      NoDistinctionNode[0] = CFLOBDDNodeHandle(new CFLOBDDDontCareNode());
    }
    else {
      CFLOBDDInternalNode *n;
      //CFLOBDDReturnMapHandle m1, m2;

      n = new CFLOBDDInternalNode(i);
      n->AConnection.entryPointHandle = &(NoDistinctionNode[i-1]);
      //m1.AddToEnd(0);
      //m1.Canonicalize();
	  n->AConnection.returnMapHandle = commonly_used_return_maps[0];//m1
  
      n->numBConnections = 1;
      n->BConnection = new Connection[1];
      n->BConnection[0].entryPointHandle = &(NoDistinctionNode[i-1]);
      //m2.AddToEnd(0);
      //m2.Canonicalize();
	  n->BConnection[0].returnMapHandle = commonly_used_return_maps[0];//m2
      n->numExits = 1;
//#ifdef PATH_COUNTING_ENABLED
      n->InstallPathCounts();
//#endif
      NoDistinctionNode[i] = CFLOBDDNodeHandle(n);
    }
  }

  CFLOBDDForkNodeHandle = CFLOBDDNodeHandle(new CFLOBDDForkNode());
  CFLOBDDDontCareNodeHandle = NoDistinctionNode[0];
  return;
} // InitNoDistinctionTable


// InitAdditionInterleavedTable
//
// Create and record the entries for "AdditionInterleavedNoCarryNode" and "AdditionInterleavedCarryNode"(each have one entry per level).
// As a side effect, these nodes are entered in the canonicalNodeTable.
//
// Note that the tables have CFLOBDDMaxLevel-1 entries, corresponding to levels 2, ..., CFLOBDDMaxLevel
//
void CFLOBDDNodeHandle::InitAdditionInterleavedTable()
{
	AdditionInterleavedNoCarryNode = new CFLOBDDNodeHandle[CFLOBDDMaxLevel - 1];
	AdditionInterleavedCarryNode   = new CFLOBDDNodeHandle[CFLOBDDMaxLevel - 1];

	CFLOBDDInternalNode *n_no_carry;
	CFLOBDDInternalNode *n_carry;

	//CFLOBDDReturnMapHandle m0, m1, m01, m10;
	CFLOBDDReturnMapHandle m02, m12, m21, m20, m012, m102;
	/*m0.AddToEnd(0);
	m0.Canonicalize();

	m1.AddToEnd(1);
	m1.Canonicalize();

	m01.AddToEnd(0);
	m01.AddToEnd(1);
	m01.Canonicalize();
	*/
	m02.AddToEnd(0);
	m02.AddToEnd(2);
	m02.Canonicalize();
	/*
	m10.AddToEnd(1);
	m10.AddToEnd(0);
	m10.Canonicalize();
	*/
	m12.AddToEnd(1);
	m12.AddToEnd(2);
	m12.Canonicalize();

	m21.AddToEnd(2);
	m21.AddToEnd(1);
	m21.Canonicalize();

	m20.AddToEnd(2);
	m20.AddToEnd(0);
	m20.Canonicalize();

	m012.AddToEnd(0);
	m012.AddToEnd(1);
	m012.AddToEnd(2);
	m012.Canonicalize();

	m102.AddToEnd(1);
	m102.AddToEnd(0);
	m102.AddToEnd(2);
	m102.Canonicalize();

	for (unsigned int i = 0; i < CFLOBDDMaxLevel-1; i++) {  // Note: i+2 is the node level
		unsigned int level = i + 2;
		if (i == 0) {
			CFLOBDDInternalNode *n1 = new CFLOBDDInternalNode(level-1);
			n1->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
			n1->numBConnections = 2;
			n1->BConnection = new Connection[n1->numBConnections];
			n1->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
			n1->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);
			n1->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			n1->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nh1(n1);

			CFLOBDDInternalNode *n2 = new CFLOBDDInternalNode(level-1);
			n2->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
			n2->numBConnections = 2;
			n2->BConnection = new Connection[n2->numBConnections];
			n2->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
			n2->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[1]);//m1
			n2->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
			n2->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nh2(n2);

			// Create the no-carry node
			n_no_carry = new CFLOBDDInternalNode(level);
			n_no_carry->AConnection = Connection(nh1, m012);
			n_no_carry->numBConnections = 3;
			n_no_carry->BConnection = new Connection[n_no_carry->numBConnections];
			n_no_carry->BConnection[0] = Connection(nh2, commonly_used_return_maps[2]);//01
			n_no_carry->BConnection[1] = Connection(nh2, commonly_used_return_maps[3]);//10
			n_no_carry->BConnection[2] = Connection(nh2, m21);
			n_no_carry->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			n_no_carry->InstallPathCounts();
#endif
			AdditionInterleavedNoCarryNode[0] = CFLOBDDNodeHandle(n_no_carry);

			// Create the carry node
			n_carry = new CFLOBDDInternalNode(level);
			n_carry->AConnection = Connection(nh1, m012);
			n_carry->numBConnections = 3;
			n_carry->BConnection = new Connection[n_carry->numBConnections];
			n_carry->BConnection[0] = Connection(nh2, commonly_used_return_maps[2]);//m01
			n_carry->BConnection[1] = Connection(nh2, m20);
			n_carry->BConnection[2] = Connection(nh2, m02);
			n_carry->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			n_carry->InstallPathCounts();
#endif
			AdditionInterleavedCarryNode[0] = CFLOBDDNodeHandle(n_carry);
		} // if (i == 0)
		else {  // i > 0
			CFLOBDDNodeHandle c0 = AdditionInterleavedNoCarryNode[i-1];
			CFLOBDDNodeHandle c1 = AdditionInterleavedCarryNode[i-1];

			// Create the no-carry node
			n_no_carry = new CFLOBDDInternalNode(level);
			n_no_carry->AConnection = Connection(c0, m012);
			n_no_carry->numBConnections = 3;
			n_no_carry->BConnection = new Connection[n_no_carry->numBConnections];
			n_no_carry->BConnection[0] = Connection(c0, m012);
			n_no_carry->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], commonly_used_return_maps[1]);//m1
			n_no_carry->BConnection[2] = Connection(c1, m102);
			n_no_carry->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			n_no_carry->InstallPathCounts();
#endif
			AdditionInterleavedNoCarryNode[i] = CFLOBDDNodeHandle(n_no_carry);

			// Create the carry node
			n_carry = new CFLOBDDInternalNode(level);
			n_carry->AConnection = Connection(c1, m012);
			n_carry->numBConnections = 3;
			n_carry->BConnection = new Connection[n_carry->numBConnections];
			n_carry->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], commonly_used_return_maps[0]);//m0
			n_carry->BConnection[1] = Connection(c0, m102);
			n_carry->BConnection[2] = Connection(c1, m012);
			n_carry->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			n_carry->InstallPathCounts();
#endif
			AdditionInterleavedCarryNode[i] = CFLOBDDNodeHandle(n_carry);
		}
	}
	return;
} // InitAdditionInterleavedTable

// Constructors/Destructor -------------------------------------------

// Default constructor
CFLOBDDNodeHandle::CFLOBDDNodeHandle()
  : handleContents(NULL)
{
}

// Constructor
//
// Construct and canonicalize
//
CFLOBDDNodeHandle::CFLOBDDNodeHandle(CFLOBDDNode *n)
  : handleContents(n)
{
  assert(n != NULL);
  handleContents->IncrRef();
  Canonicalize();
}

// Copy constructor
CFLOBDDNodeHandle::CFLOBDDNodeHandle(const CFLOBDDNodeHandle &c)
{
  handleContents = c.handleContents;
  if (handleContents != NULL) {
	// std::cout << handleContents << std::endl;
    handleContents->IncrRef();
	// std::cout << "handleContents" << std::endl;
  }
}

CFLOBDDNodeHandle::~CFLOBDDNodeHandle()
{
  if (handleContents != NULL) {
    handleContents->DecrRef();
  }
}

// Hash
unsigned int CFLOBDDNodeHandle::Hash(unsigned int modsize)
{
  return ((unsigned int) reinterpret_cast<uintptr_t>(handleContents) >> 2) % modsize;
}

// Overloaded !=
bool CFLOBDDNodeHandle::operator!= (const CFLOBDDNodeHandle & C)
{
  return handleContents != C.handleContents;
}

// Overloaded ==
bool CFLOBDDNodeHandle::operator== (const CFLOBDDNodeHandle & C) const
{
  return handleContents == C.handleContents;
}

// Overloaded assignment
CFLOBDDNodeHandle & CFLOBDDNodeHandle::operator= (const CFLOBDDNodeHandle &c)
{
  if (this != &c)      // don't assign to self!
  {
    CFLOBDDNode *temp = handleContents;
    handleContents = c.handleContents;
    if (handleContents != NULL) {
      handleContents->IncrRef();
    }
    if (temp != NULL) {
      temp->DecrRef();
    }
  }
  return *this;        
}

// Reduce and its associated cache ----------------------------------

static Hashtable<CFLReduceKey, CFLOBDDNodeHandle> *reduceCache = NULL;

CFLOBDDNodeHandle CFLOBDDNodeHandle::Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, bool forceReduce)
{
	if (replacementNumExits == 1 && !forceReduce) {
		return CFLOBDDNodeHandle::NoDistinctionNode[handleContents->Level()];
	}

	if (redMapHandle.mapContents->isIdentityMap && !forceReduce) {
		return *this;
	}

	CFLOBDDNodeHandle cachedNodeHandle;
	bool isCached = reduceCache->Fetch(CFLReduceKey(*this, redMapHandle), cachedNodeHandle);
	if (isCached) {
		// std::cout << "Hit : " << handleContents->Level() << std::endl;
		return cachedNodeHandle;
	}
	else {
		// std::cout << "Miss: " << handleContents->Level() << std::endl;
		CFLOBDDNodeHandle temp;
		temp = handleContents->Reduce(redMapHandle, replacementNumExits, forceReduce);
		reduceCache->Insert(CFLReduceKey(*this, redMapHandle), temp);
		return temp;
	}
}

void CFLOBDDNodeHandle::InitReduceCache()
{
  reduceCache = new Hashtable<CFLReduceKey, CFLOBDDNodeHandle>(HASH_NUM_BUCKETS);
}

void CFLOBDDNodeHandle::DisposeOfReduceCache()
{
	//std::cout << "Reduce cache size: " << reduceCache->Size() << std::endl;
	delete reduceCache;
	reduceCache = NULL;
}

// Canonicalization --------------------------------------------
void CFLOBDDNodeHandle::Canonicalize()
{
  CFLOBDDNode *answerContents;

  if (!handleContents->IsCanonical()) {
	  unsigned int hash = canonicalNodeTable->GetHash(handleContents);
    answerContents = canonicalNodeTable->Lookup(handleContents, hash);
    if (answerContents == NULL) {
      canonicalNodeTable->Insert(handleContents, hash);
      handleContents->SetCanonical();
    }
    else {
      answerContents->IncrRef();
      handleContents->DecrRef();
      handleContents = answerContents;
    }
  }
}

// print
std::ostream& CFLOBDDNodeHandle::print(std::ostream & out) const
{
  out << *handleContents << std::endl;
  return out;
}

bool CFLOBDDNodeHandle::IsValid()
{
	return handleContents->IsValid();
}

namespace CFL_OBDD {

std::ostream& operator<< (std::ostream & out, const CFLOBDDNodeHandle &d)
{
  d.print(out);
  return(out);
}

CFLOBDDNodeHandle MkDistinction(unsigned int level, unsigned int i)
{
  if (level == 0) {
    assert(i == 0);
    return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
  }
  else {  // Create an appropriate CFLOBDDInternalNode
    CFLOBDDInternalNode *n = new CFLOBDDInternalNode(level);
    CFLOBDDReturnMapHandle m1, m2, m3;
    if (i < (unsigned int)(1 << (level-1))) { // i falls in AConnection range
      m1.AddToEnd(0);
      m1.AddToEnd(1);
      m1.Canonicalize();
      CFLOBDDNodeHandle temp;
	  temp = MkDistinction(level-1, i);
      n->AConnection = Connection(temp, m1);

      n->numBConnections = 2;
      n->BConnection = new Connection[n->numBConnections];
      m2.AddToEnd(0);
      m2.Canonicalize();
      n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m2);
      m3.AddToEnd(1);
      m3.Canonicalize();
      n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m3);
      n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
      n->InstallPathCounts();
#endif
    }
    else {         // i falls in BConnection range
      m1.AddToEnd(0);
      m1.Canonicalize();
      n->AConnection = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);

      i = i ^ (1 << (level-1));  // Mask off high-order bit for recursive call
      n->numBConnections = 1;
      n->BConnection = new Connection[n->numBConnections];
      m2.AddToEnd(0);
      m2.AddToEnd(1);
      m2.Canonicalize();
      CFLOBDDNodeHandle temp = MkDistinction(level-1, i);
      n->BConnection[0] = Connection(temp, m2);
      n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
      n->InstallPathCounts();
#endif
    }
    return CFLOBDDNodeHandle(n);
  }
} // MkDistinction


CFLOBDDNodeHandle MkAdditionInterleavedRecursive(unsigned int level, bool carry)
{
  CFLOBDDInternalNode *n;

  //CFLOBDDReturnMapHandle m0, m1, m01, m10;
  CFLOBDDReturnMapHandle m02, m12, m21, m20, m012, m102;
  /*m0.AddToEnd(0);
  m0.Canonicalize();

  m1.AddToEnd(1);
  m1.Canonicalize();

  m01.AddToEnd(0);
  m01.AddToEnd(1);
  m01.Canonicalize();*/

  m02.AddToEnd(0);
  m02.AddToEnd(2);
  m02.Canonicalize();

  /*m10.AddToEnd(1);
  m10.AddToEnd(0);
  m10.Canonicalize();*/

  m12.AddToEnd(1);
  m12.AddToEnd(2);
  m12.Canonicalize();

  m21.AddToEnd(2);
  m21.AddToEnd(1);
  m21.Canonicalize();

  m20.AddToEnd(2);
  m20.AddToEnd(0);
  m20.Canonicalize();

  m012.AddToEnd(0);
  m012.AddToEnd(1);
  m012.AddToEnd(2);
  m012.Canonicalize();

  m102.AddToEnd(1);
  m102.AddToEnd(0);
  m102.AddToEnd(2);
  m102.Canonicalize();

  assert(level >= 2);
  if (level == 2) {
    CFLOBDDInternalNode *n1 = new CFLOBDDInternalNode(1);
    n1->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
    n1->numBConnections = 2;
    n1->BConnection = new Connection[n1->numBConnections];
    n1->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
    n1->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);
    n1->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
    n1->InstallPathCounts();
#endif
    CFLOBDDNodeHandle nh1(n1);

    CFLOBDDInternalNode *n2 = new CFLOBDDInternalNode(1);
    n2->AConnection =  Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
    n2->numBConnections = 2;
    n2->BConnection = new Connection[n2->numBConnections];
    n2->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
    n2->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[1]);//m1
    n2->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
    n2->InstallPathCounts();
#endif
    CFLOBDDNodeHandle nh2(n2);

    n = new CFLOBDDInternalNode(2);
    n->AConnection = Connection(nh1, m012);
    n->numBConnections = 3;
    n->BConnection = new Connection[n->numBConnections];
    if (!carry) {
      n->BConnection[0] = Connection(nh2, commonly_used_return_maps[2]);//m01
      n->BConnection[1] = Connection(nh2, commonly_used_return_maps[3]);//m10
      n->BConnection[2] = Connection(nh2, m21);
    }
    else {
      assert(carry == 1);
      n->BConnection[0] = Connection(nh2, commonly_used_return_maps[2]);//m01
      n->BConnection[1] = Connection(nh2, m20);
      n->BConnection[2] = Connection(nh2, m02);
    }
    n->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
    n->InstallPathCounts();
#endif
  }
  else {  // level > 2; Create an appropriate CFLOBDDInternalNode
    CFLOBDDNodeHandle c0 = MkAdditionInterleavedRecursive(level-1, false);
    CFLOBDDNodeHandle c1 = MkAdditionInterleavedRecursive(level-1, true);

    n = new CFLOBDDInternalNode(level);
    n->AConnection = Connection((carry ? c1 : c0), m012);
    n->numBConnections = 3;
    n->BConnection = new Connection[n->numBConnections];
    if (!carry) {
      n->BConnection[0] = Connection(c0, m012);
      n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], commonly_used_return_maps[1]);//m1
      n->BConnection[2] = Connection(c1, m102);
    }
    else {
      assert(carry == 1);
      n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], commonly_used_return_maps[0]);//m0
      n->BConnection[1] = Connection(c0, m102);
      n->BConnection[2] = Connection(c1, m012);
    }
    n->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
    n->InstallPathCounts();
#endif
  }
  return CFLOBDDNodeHandle(n);
} // MkAdditionInterleavedRecursive

CFLOBDDNodeHandle MkParity(unsigned int level)
{
  if (level == 0) {
    return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
  }
  else {  // Create an appropriate CFLOBDDInternalNode
    CFLOBDDInternalNode *n = new CFLOBDDInternalNode(level);
    /*CFLOBDDReturnMapHandle m1, m2;
    m1.AddToEnd(0);
    m1.AddToEnd(1);
    m1.Canonicalize();*/
    CFLOBDDNodeHandle temp = MkParity(level-1);
    n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01

    n->numBConnections = 2;
    n->BConnection = new Connection[n->numBConnections];
    n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
    //m2.AddToEnd(1);
    //m2.AddToEnd(0);
    //m2.Canonicalize();
    n->BConnection[1] = Connection(temp, commonly_used_return_maps[2]);//m10

    n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
    n->InstallPathCounts();
#endif
    return CFLOBDDNodeHandle(n);
  }
} // MkParity

CFLOBDDNodeHandle MkStepOneFourth(unsigned int level)
{
  assert(level >= 1);

  CFLOBDDInternalNode *n = new CFLOBDDInternalNode(level);
  
  if (level == 1) {
    CFLOBDDReturnMapHandle m1, m2;
    m1.AddToEnd(0);
    m1.AddToEnd(1);
    m1.Canonicalize();
    n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m1);

    n->numBConnections = 2;
    n->BConnection = new Connection[n->numBConnections];
    n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m1);
    m2.AddToEnd(1);
    m2.Canonicalize();
    n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m2);
  }
  else {  // Create an appropriate CFLOBDDInternalNode
    CFLOBDDReturnMapHandle m1, m2, m3;
    m1.AddToEnd(0);
    m1.AddToEnd(1);
    m1.Canonicalize();
    CFLOBDDNodeHandle temp = MkStepOneFourth(level-1);
    n->AConnection = Connection(temp, m1);

    n->numBConnections = 2;
    n->BConnection = new Connection[n->numBConnections];
    m2.AddToEnd(0);
    m2.Canonicalize();
    n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m2);
    m3.AddToEnd(1);
    m3.Canonicalize();
    n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m3);
  }
  n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
  n->InstallPathCounts();
#endif
  return CFLOBDDNodeHandle(n);
} // MkStepOneFourth

double ComputeProbabilityNode(CFLOBDDNodeHandle g, std::vector<double>& var_probs, std::vector<double>& path_probs, int start, int end){
	if (g.handleContents->level == 0){
		if (g == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle){
				return (1 - var_probs[start]) * path_probs[0] + var_probs[start] * path_probs[1];
		} else {
			return path_probs[0];
		}
	} else if (g == CFLOBDDNodeHandle::NoDistinctionNode[g.handleContents->level]){
			return path_probs[0];
	} else {
		CFLOBDDInternalNode* gh = (CFLOBDDInternalNode *) g.handleContents;
		std::vector<double> AConnection_PathProbs;
		for (int i = 0; i < gh->numBConnections; i++){
			std::vector<double> BConnection_PathProbs;
			for (int j = 0; j < gh->BConnection[i].returnMapHandle.Size(); j++){
				BConnection_PathProbs.push_back(path_probs[gh->BConnection[i].returnMapHandle[j]]);
			}
			double prob = ComputeProbabilityNode(*(gh->BConnection[i].entryPointHandle), var_probs, BConnection_PathProbs, (end - start)/2 + 1 + start, end);
			AConnection_PathProbs.push_back(prob);
		}

		double AProb = ComputeProbabilityNode(*(gh->AConnection.entryPointHandle), var_probs, AConnection_PathProbs, start, (end - start)/2 + start);
		return AProb;
	}
}

std::vector<double> ComputeProbabilityOfListNode(CFLOBDDNodeHandle g, std::vector<std::vector<double>>& var_probs, std::vector<std::vector<double>>& path_probs, int start, int end){
	if (g.handleContents->level == 0){
		if (g == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle){
			std::vector<double> ans (path_probs[0].size(), 0);
			for (int i = 0; i < path_probs[0].size(); i++)
				ans[i] = (1 - var_probs[start][i]) * path_probs[0][i] + var_probs[start][i] * path_probs[1][i];
			return ans;
		} else {
			std::vector<double> ans (path_probs[0].size(), 0);
			for (int i = 0; i < path_probs[0].size(); i++)
				ans[i] = path_probs[0][i];
			return ans;
		}
	} else if (g == CFLOBDDNodeHandle::NoDistinctionNode[g.handleContents->level]){
		std::vector<double> ans (path_probs[0].size(), 0);
		for (int i = 0; i <  path_probs[0].size(); i++)
			ans[i] = path_probs[0][i];
		return ans;
	} else {
		CFLOBDDInternalNode* gh = (CFLOBDDInternalNode *) g.handleContents;
		std::vector<std::vector<double>> AConnection_PathProbs;
		for (int i = 0; i < gh->numBConnections; i++){
			std::vector<std::vector<double>> BConnection_PathProbs;
			for (int j = 0; j < gh->BConnection[i].returnMapHandle.Size(); j++){
				BConnection_PathProbs.push_back(path_probs[gh->BConnection[i].returnMapHandle[j]]);
			}
			std::vector<double> prob = ComputeProbabilityOfListNode(*(gh->BConnection[i].entryPointHandle), var_probs, BConnection_PathProbs, (end - start)/2 + 1 + start, end);
			AConnection_PathProbs.push_back(prob);
		}

		std::vector<double> AProb = ComputeProbabilityOfListNode(*(gh->AConnection.entryPointHandle), var_probs, AConnection_PathProbs, start, (end - start)/2 + start);
		return AProb;
	}
}

// TODO
// std::vector<double> ComputeEntropyOfListNode(CFLOBDDNodeHandle g, std::vector<std::vector<double>>& var_probs, std::vector<std::vector<double>>& path_probs, std::vector<std::vector<double>>& entropy, int start, int end){
// 	if (g.handleContents->level == 0){
// 		if (g == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle){
// 			std::vector<double> ans (path_probs[0].size(), 0);
// 			for (int i = 0; i < path_probs[0].size(); i++){
// 				ans[i] = (1 - var_probs[start][i]) * path_probs[0][i] + var_probs[start][i] * path_probs[1][i];
// 			}
// 			return ans;
// 		} else {
// 			std::vector<double> ans (path_probs[0].size(), 0);
// 			for (int i = 0; i < path_probs[0].size(); i++){
// 				if (path_probs[0][i] == 0)
// 					ans[i] = 0;
// 				else{
// 					double v1 = (1 - var_probs[start][i]) == 0 ? 0 : (1 - var_probs[start][i]) * log2((1 - var_probs[start][i]));
// 					double v2 = var_probs[start][i] == 0 ? 0 : var_probs[start][i] * log2(var_probs[start][i]);
// 					ans[i] = path_probs[0][i] * (log2(path_probs[0][i]) + v1 + v2 - entropy[0][i]);
// 				}
// 			}
// 			return ans;
// 		}
// 	} else if (g == CFLOBDDNodeHandle::NoDistinctionNode[g.handleContents->level]){
// 		std::vector<double> ans (path_probs[0].size(), 0);
// 		for (int i = 0; i <  path_probs[0].size(); i++)
// 			ans[i] = path_probs[0][i];
// 		return ans;
// 	} else {
// 		CFLOBDDInternalNode* gh = (CFLOBDDInternalNode *) g.handleContents;
// 		std::vector<std::vector<double>> AConnection_PathProbs;
// 		for (int i = 0; i < gh->numBConnections; i++){
// 			std::vector<std::vector<double>> BConnection_PathProbs;
// 			for (int j = 0; j < gh->BConnection[i].returnMapHandle.Size(); j++){
// 				BConnection_PathProbs.push_back(path_probs[gh->BConnection[i].returnMapHandle[j]]);
// 			}
// 			std::vector<double> prob = ComputeProbabilityOfListNode(*(gh->BConnection[i].entryPointHandle), var_probs, BConnection_PathProbs, (end - start)/2 + 1 + start, end);
// 			AConnection_PathProbs.push_back(prob);
// 		}

// 		std::vector<double> AProb = ComputeProbabilityOfListNode(*(gh->AConnection.entryPointHandle), var_probs, AConnection_PathProbs, start, (end - start)/2 + start);
// 		return AProb;
// 	}
// }


}

//
// binQuotient
//
// Return v/(2**(2**level))
//
unsigned int binQuotient(unsigned int v, unsigned int level) {
  assert(level < 5);
  //  return v >> (1 << level);
  return v/(1 << (1 << level));
}

//
// binRemainder
//
// Return v mod (2**(2**level))
//
unsigned int binRemainder(unsigned int v, unsigned int level) {
  assert(level < 5);
  //  return v & (~((~0) << (1 << level)));
  return v % (1 << (1 << level));
}

namespace CFL_OBDD {
#ifdef ARBITRARY_STEP_FUNCTIONS
	CFLOBDDNodeHandle MkStepNode(unsigned int level, unsigned int left, unsigned int middle, unsigned int right)
	{
		assert(level <= 5);    // Need LONG ints
		assert(middle == 0 || middle == 1);

		if (level == 0) {
			if (left == 2 || right == 2) {
				return CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle;
			}
			else {
				return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
			}
		}
		else {  // Create an appropriate CFLOBDDInternalNode
			CFLOBDDInternalNode *n = new CFLOBDDInternalNode(level);
			int numberOfExits = (left != 0) + (middle != 0) + (right != 0);
			unsigned int a = binQuotient(left, level-1);
			unsigned int b = binQuotient(binRemainder(left, level-1) + (middle * (1<<(1<<(level-1)))) + binRemainder(right, level-1), level-1);
			unsigned int c = binQuotient(right, level-1);

			// Create AConnection -----------------------------------
			int numberOfAExits = (a != 0) + (b != 0) + (c != 0);
			assert(numberOfAExits != 0);

			CFLOBDDReturnMapHandle m1;
			for (int i = 0; i < numberOfAExits; i++) {
				m1.AddToEnd(i);
			}
			m1.Canonicalize();
			CFLOBDDNodeHandle temp = MkStepNode(level-1, a, b, c);
			n->AConnection = Connection(temp, m1);

			// Create the BConnections ------------------------------
			n->numBConnections = numberOfAExits;
			n->BConnection = new Connection[n->numBConnections];
			int curConnection = 0;
			if (a != 0) {
				CFLOBDDReturnMapHandle m2;
				m2.AddToEnd(0);  // Connect to first exit
				m2.Canonicalize();
				n->BConnection[curConnection++] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m2);
			}
			if (b != 0) {
				unsigned int aa = binRemainder(left, level-1);
				unsigned int cc = binRemainder(right, level-1);
				unsigned int bb = middle;
				CFLOBDDReturnMapHandle m3;
				if (aa != 0) {
					m3.AddToEnd(0);  // Connect to first exit
				}
				if (bb != 0) {
					if (left != 0) {
						m3.AddToEnd(1);
					}
					else {
						m3.AddToEnd(0);
					}
				}
				if (cc != 0) {
					m3.AddToEnd(numberOfExits-1);  // Connect to last exit
				}
				m3.Canonicalize();
				CFLOBDDNodeHandle temp = MkStepNode(level-1, aa, bb, cc);
				n->BConnection[curConnection++] = Connection(temp, m3);
			}
			if (c != 0) {
				CFLOBDDReturnMapHandle m4;
				m4.AddToEnd(numberOfExits-1);  // Connect to last exit
				m4.Canonicalize();
				n->BConnection[curConnection++] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m4);
			}
			n->numExits = numberOfExits;
#ifdef PATH_COUNTING_ENABLED
			n->InstallPathCounts();
#endif
			return CFLOBDDNodeHandle(n);
		}
	} // MkStepNode
#endif

	//
	// shiftAtoB
	//
	// Special operation used for Karatsuba multiplication with the
	// ASCENDING_INTERLEAVED variable ordering
	//
	// Using recursion, find (in the AConnection) the node n at level "levelAtWhichToShift",
	// and perform the following transfer of information from n's AConnection to
	// the n's BConnection: create a node whose AConnection is
	// NoDistinctionNode[levelAtWhichToShift-1] and the BConnection is the
	// original AConnection of n.
	//
	// During the recursion down AConnections to n, at each node m,
	// check that each BConnection of m is a NoDistinctionNode (of the appropriate level).
	//
	// Recursive version, so works for shifts at levels other than at top level.
	//
	// TODO: Checks on g == NoDistinctionNode[g->level] and function caching
	//
	CFLOBDDNodeHandle shiftAtoB(CFLOBDDNodeHandle f, const unsigned int levelAtWhichToShift)
	{
		assert(f.handleContents->level > 0);
		CFLOBDDInternalNode *g = (CFLOBDDInternalNode *)(f.handleContents);

		// Check that each BConnection of *g is a NoDistinctionNode[level-1]
		for (unsigned int j = 0; j < g->numBConnections; j++) {
			if (*(g->BConnection[j].entryPointHandle) != CFLOBDDNodeHandle::NoDistinctionNode[g->level - 1]) {
				std::cout << "g->level = " << g->level << std::endl;
				std::cout << *(g->BConnection[j].entryPointHandle) << std::endl << std::endl;
				std::cout << CFLOBDDNodeHandle::NoDistinctionNode[g->level - 1] << std::endl << std::endl;
			}
			assert(*(g->BConnection[j].entryPointHandle) == CFLOBDDNodeHandle::NoDistinctionNode[g->level - 1]);
		}

		// Create an appropriate CFLOBDDInternalNode
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(g->level);
		if (g->level == levelAtWhichToShift) {  // Base case: the current level is consistent with shift
			// Put a NoDistinctionNode in the AConnection
			CFLOBDDReturnMapHandle m1;
			m1.AddToEnd(0);
			m1.Canonicalize();
			n->AConnection = Connection(CFLOBDDNodeHandle::NoDistinctionNode[g->level - 1], m1);

			// Transfer g's AConnection to the BConnection
			n->numBConnections = 1;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = g->AConnection;
		}
		else {    // Haven't reached the correct level yet, so . . .
			// Apply shiftAtoB recursively to g's AConnection
			CFLOBDDNodeHandle temp = shiftAtoB(*(g->AConnection.entryPointHandle), levelAtWhichToShift);
			n->AConnection = Connection(temp, g->AConnection.returnMapHandle);

			// Copy over all BConnections from g
			n->numBConnections = g->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				n->BConnection[j] = g->BConnection[j];
			}
		}
		n->numExits = g->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	//
	// shiftBtoA
	//
	// Special operation used for Karatsuba multiplication with the
	// DESCENDING_INTERLEAVED variable ordering
	//
	// Using recursion, find (in the BConnection) the node n at level "levelAtWhichToShift",
	// and perform the following transfer of information from n's BConnection to
	// the n's AConnection: create a node whose BConnection is
	// NoDistinctionNode[levelAtWhichToShift-1] and the AConnection is the
	// original BConnection of n.
	//
	// During the recursion down to n, at each node m,
	// check that the AConnection of m is a NoDistinctionNode (of the appropriate level).
	//
	// Recursive version, so works for shifts at levels other than at top level.
	//
	// TODO: Checks on g == NoDistinctionNode[g->level] and function caching
	//
	CFLOBDDNodeHandle shiftBtoA(CFLOBDDNodeHandle f, const unsigned int levelAtWhichToShift)
	{
		assert(f.handleContents->level > 0);
		CFLOBDDInternalNode *g = (CFLOBDDInternalNode *)(f.handleContents);

		if (f == CFLOBDDNodeHandle::NoDistinctionNode[g->level]) {
			return f;
		}

		// Check that the AConnection of *g is a NoDistinctionNode[level-1]
		if (*(g->AConnection.entryPointHandle) != CFLOBDDNodeHandle::NoDistinctionNode[g->level - 1]) {
			std::cout << "g->level = " << g->level << std::endl;
			std::cout << *(g->AConnection.entryPointHandle) << std::endl << std::endl;
			std::cout << CFLOBDDNodeHandle::NoDistinctionNode[g->level - 1] << std::endl << std::endl;
		}
		assert(*(g->AConnection.entryPointHandle) == CFLOBDDNodeHandle::NoDistinctionNode[g->level - 1]);

		// Create an appropriate CFLOBDDInternalNode
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(g->level);
		if (g->level == levelAtWhichToShift) {  // Base case: the current level is consistent with shift
			// Transfer g's BConnection[0] to the AConnection
			n->AConnection = g->BConnection[0];

			// Put a NoDistinctionNode in each BConnection
			n->numBConnections = g->numExits;
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int j = 0; j < n->numBConnections; j++) {
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(j);
				m1.Canonicalize();
				n->BConnection[j] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[g->level - 1], m1);
			}
		}
		else {    // Haven't reached the correct level yet, so . . .
			// Copy over the AConnection from g
			n->AConnection = g->AConnection;

			// Apply shiftBtoA recursively to g's BConnection[0]
			CFLOBDDNodeHandle temp = shiftBtoA(*(g->BConnection[0].entryPointHandle), levelAtWhichToShift);
			n->numBConnections = g->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, g->BConnection[0].returnMapHandle);
		}
		n->numExits = g->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	//
	// shiftAtoBAtLevelOne
	//
	// Given a CFLOBDDNodeHandle for which all BConnections at level 1 are
	// NoDistinctionNode[0], create a corresponding CFLOBDDNodeHandle
	// in which
	//  1. The AConnection is moved to BConnection[0], and
	//  2. All AConnections at level 1 are NoDistinctionNode[0]
	//
	// Precondition: *this is a CFLOBDDNodeHandle for which all BConnections at level 1 are
	//               NoDistinctionNode[0]
	// Postcondition: *ans is a CFLOBDDNodeHandle that satisfies conditions 1 and 2
	//
	// TODO: Checks on g == NoDistinctionNode[g->level] and function caching
	//
	// Hashtable<CFLOBDDNodeHandle, CFLOBDDNodeHandle> *memoTable
	//
	CFLOBDDNodeHandle shiftAtoBAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, CFLOBDDNodeHandle f)
	{
		assert(f.handleContents->level > 0);
		CFLOBDDInternalNode *g = (CFLOBDDInternalNode *)(f.handleContents);

		if (f == CFLOBDDNodeHandle::NoDistinctionNode[g->level]) {
			return f;
		}
		totalVisitCount++;
		if (visitedNodes->Lookup(g) == NULL) {
			visitedNodes->Insert(g);
		}
		else {
			redundantVisitCount++;
		}

		// Create an appropriate CFLOBDDInternalNode
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(g->level);
		if (g->level == 1) {  // Base case: the current level is consistent with shift
			// Check that each BConnection of *g is a NoDistinctionNode[0]
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				if (*(g->BConnection[j].entryPointHandle) != CFLOBDDNodeHandle::NoDistinctionNode[0]) {
					std::cout << "g->level = " << g->level << std::endl;
					std::cout << *(g->BConnection[j].entryPointHandle) << std::endl << std::endl;
					std::cout << CFLOBDDNodeHandle::NoDistinctionNode[0] << std::endl << std::endl;
				}
				assert(*(g->BConnection[j].entryPointHandle) == CFLOBDDNodeHandle::NoDistinctionNode[0]);
			}

			// Put a NoDistinctionNode[0] in the AConnection
			CFLOBDDReturnMapHandle m1;
			m1.AddToEnd(0);
			m1.Canonicalize();
			n->AConnection = Connection(CFLOBDDNodeHandle::NoDistinctionNode[0], m1);

			// Transfer g's AConnection to BConnection[0]
			n->numBConnections = 1;
			n->BConnection = new Connection[1];
			n->BConnection[0] = g->AConnection;
		}
		else {    // Haven't reached the correct level yet, so apply shiftAtoBAtLevelOne recursively
			// Create the AConnection
			CFLOBDDNodeHandle temp1 = shiftAtoBAtLevelOne(visitedNodes, totalVisitCount, redundantVisitCount, *(g->AConnection.entryPointHandle));
			n->AConnection = Connection(temp1, g->AConnection.returnMapHandle);

			// Create the BConnections
			n->numBConnections = g->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				CFLOBDDNodeHandle temp2 = shiftAtoBAtLevelOne(visitedNodes, totalVisitCount, redundantVisitCount, *(g->BConnection[j].entryPointHandle));
				n->BConnection[j] = Connection(temp2, g->BConnection[j].returnMapHandle);
			}
		}
		n->numExits = g->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	//
	// shiftBtoAAtLevelOne
	//
	// Given a CFLOBDDNodeHandle for which all AConnections at level 1 are
	// NoDistinctionNode[0], create a corresponding CFLOBDDNodeHandle
	// in which
	//  1. The BConnection[0] is moved to AConnection, and
	//  2. All BConnections at level 1 are NoDistinctionNode[0]
	//
	// Precondition: *this is a CFLOBDDNodeHandle for which all AConnections at level 1 are
	//               NoDistinctionNode[0]
	// Postcondition: *ans is a CFLOBDDNodeHandle that satisfies conditions 1 and 2
	//
	// TODO: Checks on g == NoDistinctionNode[g->level] and function caching
	//
	CFLOBDDNodeHandle shiftBtoAAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, CFLOBDDNodeHandle f)
	{
		assert(f.handleContents->level > 0);
		CFLOBDDInternalNode *g = (CFLOBDDInternalNode *)(f.handleContents);

		if (f == CFLOBDDNodeHandle::NoDistinctionNode[g->level]) {
			return f;
		}
		totalVisitCount++;
		if (visitedNodes->Lookup(g) == NULL) {
			visitedNodes->Insert(g);
		}
		else {
			redundantVisitCount++;
		}

		// Create an appropriate CFLOBDDInternalNode
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(g->level);
		if (g->level == 1) {  // Base case: the current level is consistent with shift
			// Check that the AConnection of *g is a NoDistinctionNode[0]
			if (*(g->AConnection.entryPointHandle) != CFLOBDDNodeHandle::NoDistinctionNode[0]) {
				std::cout << "g->level = " << g->level << std::endl;
				std::cout << *(g->AConnection.entryPointHandle) << std::endl << std::endl;
				std::cout << CFLOBDDNodeHandle::NoDistinctionNode[0] << std::endl << std::endl;
			}
			assert(*(g->AConnection.entryPointHandle) == CFLOBDDNodeHandle::NoDistinctionNode[0]);

			// Transfer g's BConnection[0] to the AConnection
			n->AConnection = g->BConnection[0];

			// Put a NoDistinctionNode in each BConnection
			n->numBConnections = g->numExits;
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int j = 0; j < n->numBConnections; j++) {
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(j);
				m1.Canonicalize();
				n->BConnection[j] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[0], m1);
			}
		}
		else {    // Haven't reached the correct level yet, so apply shiftBtoAAtLevelOne recursively
			// Create the AConnection
			CFLOBDDNodeHandle temp1 = shiftBtoAAtLevelOne(visitedNodes, totalVisitCount, redundantVisitCount, *(g->AConnection.entryPointHandle));
			n->AConnection = Connection(temp1, g->AConnection.returnMapHandle);

			// Create the BConnections
			n->numBConnections = g->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				CFLOBDDNodeHandle temp2 = shiftBtoAAtLevelOne(visitedNodes, totalVisitCount, redundantVisitCount, *(g->BConnection[j].entryPointHandle));
				n->BConnection[j] = Connection(temp2, g->BConnection[j].returnMapHandle);
			}
		}
		n->numExits = g->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	//
	// duplicateAinBAtLevelOne
	//
	// Given a CFLOBDDNodeHandle for which all BConnections at level 1 are
	// NoDistinctionNode[0], create a corresponding CFLOBDDNodeHandle
	// in which
	//  1. The AConnection is duplicated in BConnection[0], and
	//  2. The AConnection and BConnection[0] "match"
	// "Match" means that if the AConnection's handle is CFLOBDDNodeHandle::CFLOBDDForkNodeHandle
	// then the three assignments of the level-one Boolean variables (0,0), (0,1), and (1,0)
	// are routed to exit[0] and (1,1) is routed to exit[1].
	//
	// Precondition: *this is a CFLOBDDNodeHandle for which all BConnections at level 1 are
	//               NoDistinctionNode[0]
	// Postcondition: *ans is a CFLOBDDNodeHandle that satisfies conditions 1 and 2
	//
	// TODO: Checks on g == NoDistinctionNode[g->level] and function caching
	//
	CFLOBDDNodeHandle duplicateAinBAtLevelOne(CFLOBDDNodeHandle f)
	{
		assert(f.handleContents->level > 0);
		CFLOBDDInternalNode *g = (CFLOBDDInternalNode *)(f.handleContents);

		// Create an appropriate CFLOBDDInternalNode
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(g->level);
		if (g->level == 1) {  // Base case: the level at which to create a node with duplicated AConnection
			// Check that each BConnection of *g is a NoDistinctionNode[0]
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				if (*(g->BConnection[j].entryPointHandle) != CFLOBDDNodeHandle::NoDistinctionNode[0]) {
					std::cout << "g->level = " << g->level << std::endl;
					std::cout << g->BConnection[j].entryPointHandle << std::endl << std::endl;
					std::cout << CFLOBDDNodeHandle::NoDistinctionNode[0] << std::endl << std::endl;
				}
				assert(*(g->BConnection[j].entryPointHandle) == CFLOBDDNodeHandle::NoDistinctionNode[0]);
			}

			// Transfer g's AConnection to n
			n->AConnection = g->AConnection;

			// Create an appropriate set of BConnections
			n->numBConnections = g->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			if (g->numBConnections == 1) {
				assert(*(g->AConnection.entryPointHandle) == CFLOBDDNodeHandle::NoDistinctionNode[0]);
				// Transfer g's BConnection[0] to n; n should be equal to CFLOBDDNodeHandle::NoDistinctionNode[1]
				n->BConnection[0] = g->AConnection;
			}
			else {
				assert(g->numBConnections == 2 && *(g->AConnection.entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle);
				// Put a NoDistinctionNode[0] in BConnection[0]
				CFLOBDDReturnMapHandle m0;
				m0.AddToEnd(0);
				m0.Canonicalize();
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[0], m0);

				// Put a ForkNode in BConnection[1]
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m01.Canonicalize();
				n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				assert(n->BConnection[1] == g->AConnection);
			}
		}
		else {    // Haven't reached the correct level yet, so apply duplicateAinBAtLevelOne recursively
			// Create the AConnection
			CFLOBDDNodeHandle temp1 = duplicateAinBAtLevelOne(*(g->AConnection.entryPointHandle));
			n->AConnection = Connection(temp1, g->AConnection.returnMapHandle);

			// Create the BConnections
			n->numBConnections = g->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				CFLOBDDNodeHandle temp2 = duplicateAinBAtLevelOne(*(g->BConnection[j].entryPointHandle));
				n->BConnection[j] = Connection(temp2, g->BConnection[j].returnMapHandle);
			}
		}
		n->numExits = g->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

}  // namespace CFL_OBDD

//********************************************************************
// CFLReduceKey
//********************************************************************

// Constructor
CFLReduceKey::CFLReduceKey(CFLOBDDNodeHandle nodeHandle, ReductionMapHandle redMapHandle)
  :  nodeHandle(nodeHandle), redMapHandle(redMapHandle)
{
}

// Hash
unsigned int CFLReduceKey::Hash(unsigned int modsize)
{
  unsigned int hvalue = 0;
  hvalue = (997 * nodeHandle.Hash(modsize) + redMapHandle.Hash(modsize)) % modsize;
  return hvalue;
}

// print
std::ostream& CFLReduceKey::print(std::ostream & out) const
{
  out << "(" << nodeHandle << ", " << redMapHandle << ")";
  return out;
}

std::ostream& operator<< (std::ostream & out, const CFLReduceKey &p)
{
  p.print(out);
  return(out);
}

CFLReduceKey& CFLReduceKey::operator= (const CFLReduceKey& i)
{
  if (this != &i)      // don't assign to self!
  {
    nodeHandle = i.nodeHandle;
    redMapHandle = i.redMapHandle;
  }
  return *this;        
}

// Overloaded !=
bool CFLReduceKey::operator!=(const CFLReduceKey& p)
{
  return (nodeHandle != p.nodeHandle) || (redMapHandle != p.redMapHandle);
}

// Overloaded ==
bool CFLReduceKey::operator==(const CFLReduceKey& p)
{
  return (nodeHandle == p.nodeHandle) && (redMapHandle == p.redMapHandle);
}

//********************************************************************
// CFLOBDDNode
//********************************************************************

// Initializations of static members ---------------------------------

unsigned int const CFLOBDDNode::maxLevel = CFLOBDDMaxLevel;

// Constructors/Destructor -------------------------------------------

// Default constructor
CFLOBDDNode::CFLOBDDNode()
	: level(maxLevel), refCount(0), isCanonical(false), isValid(false), isNumPathsMemAllocated(false)
{
}

// Constructor
CFLOBDDNode::CFLOBDDNode(const unsigned int l)
	: level(l), refCount(0), isCanonical(false), isValid(false), isNumPathsMemAllocated(false)
{
}

CFLOBDDNode::~CFLOBDDNode()
{
}

// print
namespace CFL_OBDD {
std::ostream& operator<< (std::ostream & out, const CFLOBDDNode &n)
{
  n.print(out);
  return(out);
}
}

//********************************************************************
// CFLOBDDInternalNode
//********************************************************************

// Constructors/Destructor -------------------------------------------

CFLOBDDInternalNode::CFLOBDDInternalNode(const unsigned int l)
  :  CFLOBDDNode(l)
{
}

CFLOBDDInternalNode::~CFLOBDDInternalNode()
{
  delete [] BConnection;
//#ifdef PATH_COUNTING_ENABLED
  if (isNumPathsMemAllocated)
	delete [] numPathsToExit;
//#endif
}

// print
std::ostream& CFLOBDDInternalNode::print(std::ostream & out) const
{
    unsigned int i, j;
    for (i = level; i < maxLevel; i++) {  // Indentation
        out << "  ";
    }
    out << "A: " << std::endl;
    if (CFLOBDDNodeHandle::NoDistinctionNode[level - 1] == *AConnection.entryPointHandle) {
	    for (i = level-1; i < maxLevel; i++) {  // Indentation
		    out << "  ";
	    }
	    out << "NoDistinctionNode[" << level - 1 << "]" << std::endl;
    }
    else {
	    out << *AConnection.entryPointHandle;
    }
    for (i = level; i < maxLevel; i++) {  // Indentation
        out << "  ";
    }
    out << AConnection.returnMapHandle;

    for (j = 0; j < numBConnections; j++) {
	    for (i = level; i < maxLevel; i++) {  // Indentation
	        out << "  ";
        }
        out << "B[" << j << "]:" << std::endl;
	    if (CFLOBDDNodeHandle::NoDistinctionNode[level - 1] == *BConnection[j].entryPointHandle) {
		    for (i = level - 1; i < maxLevel; i++) {  // Indentation
			    out << "  ";
		    }
			out << "NoDistinctionNode[" << level - 1 << "]" << std::endl;
        }
        else {
		    out << *BConnection[j].entryPointHandle;
	    }
		for (i = level; i < maxLevel; i++) {  // Indentation
            out << "  ";
        }
        out << BConnection[j].returnMapHandle;
    }
    return out;
}

bool CFLOBDDInternalNode::IsValid()
{
	if (isValid)
		return true;
	
	if (!(level == (AConnection.entryPointHandle->handleContents->level + 1)))
		return false;
	for (unsigned int i = 0; i < numBConnections; i++){
		if (!(level == (BConnection[i].entryPointHandle->handleContents->level + 1)))
			return false;
	}

	for (unsigned int i = 0; i < AConnection.returnMapHandle.Size(); i++){
		if (AConnection.returnMapHandle[i] != i)
			return false;
	}

	unsigned int tempNumExits = 0;
	for (unsigned int i = 0; i < numBConnections; i++){
		std::unordered_set<int> s;
		for (unsigned int j = 0; j < BConnection[i].returnMapHandle.Size(); j++){
			if (s.find(BConnection[i].returnMapHandle[j]) != s.end())
				return false;
			s.insert(BConnection[i].returnMapHandle[j]);
			if (BConnection[i].returnMapHandle[j] > tempNumExits)
				return false;
			else if (BConnection[i].returnMapHandle[j] == tempNumExits)
				tempNumExits++;
		}
	}

	if (numExits != tempNumExits)
		return false;

	if (!AConnection.entryPointHandle->handleContents->IsValid())
		return false;
	for (unsigned int i = 0; i < numBConnections; i++){
		if (!BConnection[i].entryPointHandle->handleContents->IsValid())
			return false;
	}
	
	isValid = true;
	return true;
}

void CFLOBDDInternalNode::FillSatisfyingAssignment(unsigned int exitNumber, SH_OBDD::Assignment &assignment, unsigned int &index)
{
  for (unsigned int i = 0; i < AConnection.entryPointHandle->handleContents->numExits; i++) {
    for (unsigned int j = 0; j < BConnection[i].entryPointHandle->handleContents->numExits; j++) {
      unsigned int k = BConnection[i].returnMapHandle.Lookup(j);
      if (k == exitNumber) { // Found it
        BConnection[i].entryPointHandle->handleContents->FillSatisfyingAssignment(j, assignment, index);
        AConnection.entryPointHandle->handleContents->FillSatisfyingAssignment(i, assignment, index);
        return;
      }
    }
  }
  std::cerr << "Failure in CFLOBDDInternalNode::FillSatisfyingAssignment:" << std::endl;
  std::cerr << "  exitNumber = " << exitNumber << std::endl;
  //std::cerr << "  assignment = " << assignment << std::endl;  ETTODO - Fix
  std::cerr << "  index = " << index << std::endl;
  abort();
}

int CFLOBDDInternalNode::Traverse(SH_OBDD::AssignmentIterator &ai)
{
  int ans;

  if (this == CFLOBDDNodeHandle::NoDistinctionNode[level].handleContents) {
    ai.Advance((unsigned int)(((unsigned int)1) << level));
    ans = 0;
  }
  else {
    int i, j, k;
    i = AConnection.entryPointHandle->handleContents->Traverse(ai);
    j = AConnection.returnMapHandle.Lookup(i);
    k = BConnection[j].entryPointHandle->handleContents->Traverse(ai);
    ans = BConnection[j].returnMapHandle.Lookup(k);
  }
  return ans;
}

CFLOBDDReturnMapHandle ComposeAndReduce(CFLOBDDReturnMapHandle& mapHandle, ReductionMapHandle& redMapHandle, ReductionMapHandle& inducedRedMapHandle)
{
	int c2, c3;
	int size = mapHandle.mapContents->mapArray.size();
	CFLOBDDReturnMapHandle answer;// (size);
	if (redMapHandle.mapContents->isIdentityMap){
		inducedRedMapHandle = redMapHandle;
		return mapHandle;
	}
	std::unordered_map<int, unsigned int> reductionMap (size);
	for (int i = 0; i < size; i++)
	{
		c2 = mapHandle.mapContents->mapArray[i];
		c3 = redMapHandle.Lookup(c2);
		if (reductionMap.find(c3) == reductionMap.end()){
			answer.AddToEnd(c3); 	  // Why not answer.AddToEnd(c3);
			reductionMap.emplace(c3, answer.Size() - 1);
			inducedRedMapHandle.AddToEnd(answer.Size() - 1);
		}
		else{
			inducedRedMapHandle.AddToEnd(reductionMap[c3]);
		}
	}
	inducedRedMapHandle.Canonicalize();
	answer.Canonicalize();
	return answer;
}

CFLOBDDNodeHandle CFLOBDDInternalNode::Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, bool forceReduce)
{
  CFLOBDDInternalNode *n = new CFLOBDDInternalNode(level);
  ReductionMapHandle AReductionMapHandle;        // To record duplicate BConnections

  // Reduce the B connections
     n->BConnection = new Connection[numBConnections];   // May create shorter version later
     n->numBConnections = 0;
     for (unsigned int i = 0; i < numBConnections; i++) {
        ReductionMapHandle inducedReductionMapHandle(redMapHandle.Size());
        CFLOBDDReturnMapHandle inducedReturnMap;
		inducedReturnMap = ComposeAndReduce(BConnection[i].returnMapHandle, redMapHandle, inducedReductionMapHandle);
        //CFLOBDDReturnMapHandle reducedReturnMap = BConnection[i].returnMapHandle.Compose(redMapHandle);
        //reducedReturnMap.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
        CFLOBDDNodeHandle temp = BConnection[i].entryPointHandle->Reduce(inducedReductionMapHandle, inducedReturnMap.Size(), forceReduce);
        Connection c(temp, inducedReturnMap);
        unsigned int position = n->InsertBConnection(n->numBConnections, c);
        AReductionMapHandle.AddToEnd(position);
     }
     AReductionMapHandle.Canonicalize();
     if (n->numBConnections < numBConnections) {  // Shorten
       Connection *temp = n->BConnection;
       n->BConnection = new Connection[n->numBConnections];
       for (unsigned int j = 0; j < n->numBConnections; j++) {
         n->BConnection[j] = temp[j];
       }
       delete [] temp;
     }

  // Reduce the A connection
     ReductionMapHandle inducedAReductionMapHandle;
	 CFLOBDDReturnMapHandle inducedAReturnMap;
	 inducedAReturnMap = ComposeAndReduce(AConnection.returnMapHandle, AReductionMapHandle, inducedAReductionMapHandle);
     //CFLOBDDReturnMapHandle reducedAReturnMap = AConnection.returnMapHandle.Compose(AReductionMapHandle);
     //reducedAReturnMap.InducedReductionAndReturnMap(inducedAReductionMapHandle, inducedAReturnMap);
     CFLOBDDNodeHandle tempHandle = AConnection.entryPointHandle->Reduce(inducedAReductionMapHandle, inducedAReturnMap.Size(), forceReduce);
     n->AConnection = Connection(tempHandle, inducedAReturnMap);

  // Other material that has to be filled in
     n->numExits = replacementNumExits;
#ifdef PATH_COUNTING_ENABLED
     n->InstallPathCounts();
#endif
  return CFLOBDDNodeHandle(n);
} // CFLOBDDInternalNode::Reduce

unsigned int CFLOBDDInternalNode::Hash(unsigned int modsize)
{
  unsigned int hvalue = AConnection.Hash(modsize);
  for (unsigned int j = 0; j < numBConnections; j++) {
    hvalue = (997 * hvalue + BConnection[j].Hash(modsize)) % modsize;
  }
  return hvalue;
}

void CFLOBDDInternalNode::DumpConnections(Hashset<CFLOBDDNodeHandle> *visited, std::ostream & out /* = std::cout */)
{
	if (visited->Lookup(new CFLOBDDNodeHandle(this)) == NULL) {
    unsigned int i;
	visited->Insert(new CFLOBDDNodeHandle(this));
    AConnection.entryPointHandle->handleContents->DumpConnections(visited, out);
    for (i = 0; i < numBConnections; i++) {
      BConnection[i].entryPointHandle->handleContents->DumpConnections(visited, out);
    }
    out << AConnection << std::endl;
    for (i = 0; i < numBConnections; i++) {
      out << BConnection[i] << std::endl;
    }
  }
}

void CFLOBDDInternalNode::CountNodesAndEdges(Hashset<CFLOBDDNodeHandle> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
	unsigned int &nodeCount, unsigned int &edgeCount, unsigned int &returnEdgesCount)
{
  if (visitedNodes->Lookup(new CFLOBDDNodeHandle(this)) == NULL) {
    visitedNodes->Insert(new CFLOBDDNodeHandle(this));
    nodeCount++;
    edgeCount += 2* (1 + numBConnections);
	/*AConnection.entryPointHandle->handleContents->CountNodesAndEdges(visitedNodes, visitedEdges, nodeCount, edgeCount);
	for (unsigned int i = 0; i < numBConnections; i++){
		BConnection[i].entryPointHandle->handleContents->CountNodesAndEdges(visitedNodes, visitedEdges, nodeCount, edgeCount);
		edgeCount += BConnection[i].returnMapHandle.Size();
	}*/
    if (visitedEdges->Lookup(AConnection.returnMapHandle.mapContents) == NULL) {
      visitedEdges->Insert(AConnection.returnMapHandle.mapContents);
      edgeCount += AConnection.returnMapHandle.Size();
	  returnEdgesCount += AConnection.returnMapHandle.Size();
    }
    AConnection.entryPointHandle->handleContents->CountNodesAndEdges(visitedNodes, visitedEdges, nodeCount, edgeCount, returnEdgesCount);
    for (unsigned int i = 0; i < numBConnections; i++) {
      BConnection[i].entryPointHandle->handleContents->CountNodesAndEdges(visitedNodes, visitedEdges, nodeCount, edgeCount, returnEdgesCount);
      if (visitedEdges->Lookup(BConnection[i].returnMapHandle.mapContents) == NULL) {
        visitedEdges->Insert(BConnection[i].returnMapHandle.mapContents);
        edgeCount += BConnection[i].returnMapHandle.Size();
		returnEdgesCount += BConnection[i].returnMapHandle.Size();
      }
    }
	
  }
}

void CFLOBDDInternalNode::CountNodes(Hashset<CFLOBDDNodeHandle> *visitedNodes, unsigned int &nodeCount)
{
	if (visitedNodes->Lookup(new CFLOBDDNodeHandle(this)) == NULL) {
		visitedNodes->Insert(new CFLOBDDNodeHandle(this));
		nodeCount++;
		AConnection.entryPointHandle->handleContents->CountNodes(visitedNodes, nodeCount);
		for (unsigned int i = 0; i < numBConnections; i++) {
			BConnection[i].entryPointHandle->handleContents->CountNodes(visitedNodes, nodeCount);
		}
	}
}

void CFLOBDDInternalNode::CountPaths(Hashset<CFLOBDDNodeHandle> *visitedNodes)
{
	CFLOBDDNodeHandle* handle = new CFLOBDDNodeHandle(this);
	if (visitedNodes->Lookup(handle) == NULL) {
		visitedNodes->Insert(handle);
		AConnection.entryPointHandle->handleContents->CountPaths(visitedNodes);
		for (unsigned int i = 0; i < numBConnections; i++) {
			BConnection[i].entryPointHandle->handleContents->CountPaths(visitedNodes);
		}
		InstallPathCounts();
	}
}

// Overloaded !=
bool CFLOBDDInternalNode::operator!= (const CFLOBDDNode & n)
{
  return !(*this == n);
}

// Overloaded ==
bool CFLOBDDInternalNode::operator== (const CFLOBDDNode & n)
{
  if (n.NodeKind() != CFLOBDD_INTERNAL)
    return false;
  CFLOBDDInternalNode &m = (CFLOBDDInternalNode &)n;
  if (level != m.level)
    return false;
  if (numExits != m.numExits)
    return false;
  if (numBConnections != m.numBConnections)
    return false;
  if (AConnection != m.AConnection)
    return false;
  //for (unsigned int j = numBConnections-1; j >= 0; j--) {
  for (unsigned int j = 0; j < numBConnections; j++) {
    if (BConnection[j] != m.BConnection[j])
      return false;
  }
  return true;
}

void CFLOBDDInternalNode::IncrRef()
{
  refCount++;    // Warning: Saturation not checked
}

void CFLOBDDInternalNode::DecrRef()
{
  if (--refCount == 0) {    // Warning: Saturation not checked
    if (isCanonical) {
      CFLOBDDNodeHandle::canonicalNodeTable->DeleteEq(this);
    }
    delete this;
  }
}

// Insert Connection c at BConnection position j, but only if c does not
// duplicate an existing BConnection.
// Return the position at which c was found or inserted.
unsigned int CFLOBDDInternalNode::InsertBConnection(unsigned int &j, Connection &c)
{
  for (unsigned int i = 0; i < j; i++) {
    if (BConnection[i] == c) {
      return i;
    }
  }
  BConnection[j] = c;
  j++;
  return j-1;
}

//#ifdef PATH_COUNTING_ENABLED

long double addNumPathsToExit(long double path1, long double path2){
	long double maxPath = std::max(path1, path2);
	long double minPath = std::min(path1, path2);
	if (minPath == -1.0)
		return maxPath;
	long double totalPaths = minPath;
	if ((maxPath - minPath) > 15)
		totalPaths = maxPath;
	else{
		unsigned int powerVal = pow(2, maxPath - minPath) + 1;
		totalPaths += log2l(powerVal);
	}
	return totalPaths;
}

long double addNumPathsToExit(std::vector<long double>& logOfPaths){
	if (logOfPaths.size() == 1)
		return logOfPaths.back();
	long double sum = 0.0;
	for (int i = 0; i < logOfPaths.size() - 1; i++){
		if (logOfPaths[i] != -1.0 * std::numeric_limits<long double>::infinity())
			sum += pow(2, logOfPaths[i] - logOfPaths.back());
	}
	long double logOfSum = logOfPaths.back() + log1p(sum)/((long double)log(2));
	return logOfSum;
}

// InstallPathCounts
void CFLOBDDInternalNode::InstallPathCounts()
{
  //numPathsToExit = new unsigned long long int[numExits];
  //numPathsToExit = new cpp_int[numExits];
  numPathsToExit = new long double[numExits];
  isNumPathsMemAllocated = true;
  for (unsigned int i = 0; i < numExits; i++) {
    numPathsToExit[i] = -1.0 * std::numeric_limits<long double>::infinity();
  }

  std::map<unsigned int, std::vector<long double>> storingNumPathsToExit;

  for (unsigned int i = 0; i < AConnection.entryPointHandle->handleContents->numExits; i++) {
    for (unsigned int j = 0; j < BConnection[i].entryPointHandle->handleContents->numExits; j++) {
      unsigned int k = BConnection[i].returnMapHandle.Lookup(j);
	  //std::cout << "Install Paths --------------------------------------\n";
	  //for (unsigned int l = 0; l < BConnection[i].returnMapHandle.mapContents->mapArray.size(); l++)
		 // std::cout << l << " " << BConnection[i].returnMapHandle.mapContents->mapArray[l] << " " << BConnection[i].returnMapHandle.Lookup(l) << std::endl;
	  //std::cout << i << " " << j << " " << k << std::endl;
	  ///*std::cout << (AConnection) << std::endl;
	  //std::cout << (AConnection.entryPointHandle->handleContents == NULL) << std::endl;
	  //std::cout << AConnection.entryPointHandle->handleContents->numPathsToExit[i] << std::endl;
	  //std::cout << (BConnection[i]) << std::endl;
	  //std::cout << (BConnection[i].entryPointHandle->handleContents == NULL) << std::endl;
	  //std::cout << BConnection[i].entryPointHandle->handleContents->numPathsToExit[j] << std::endl;*/
	  //std::cout << "-------------------------------------------------------\n";
	  long double numPathsValue = AConnection.entryPointHandle->handleContents->numPathsToExit[i] + BConnection[i].entryPointHandle->handleContents->numPathsToExit[j];
	  if (storingNumPathsToExit.find(k) == storingNumPathsToExit.end()){
		  std::vector<long double> logOfPaths;
		  logOfPaths.push_back(numPathsValue);
		  storingNumPathsToExit[k] = logOfPaths;
	  }
	  else{
		  storingNumPathsToExit[k].push_back(numPathsValue);
	  }
	  //numPathsToExit[k] = addNumPathsToExit(numPathsToExit[k], numPathsValue);
    }

	for (std::map<unsigned int, std::vector<long double>>::iterator it = storingNumPathsToExit.begin(); it != storingNumPathsToExit.end(); it++){
		sort(it->second.begin(), it->second.end());
		numPathsToExit[it->first] = addNumPathsToExit(it->second);
	}
  }
}
//#endif

//********************************************************************
// CFLOBDDLeafNode
//********************************************************************

// Constructors/Destructor -------------------------------------------

// Default constructor
CFLOBDDLeafNode::CFLOBDDLeafNode()
  :  CFLOBDDNode(0)
{
  refCount = 1;
}

CFLOBDDLeafNode::~CFLOBDDLeafNode()
{
}

void CFLOBDDLeafNode::DumpConnections(Hashset<CFLOBDDNodeHandle> *, std::ostream & /* = std::cout */ )
{
}

void CFLOBDDLeafNode::CountNodesAndEdges(Hashset<CFLOBDDNodeHandle> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *, unsigned int &nodeCount, 
	unsigned int &, unsigned int& )
{
  if (visitedNodes->Lookup(new CFLOBDDNodeHandle(this)) == NULL) {
    visitedNodes->Insert(new CFLOBDDNodeHandle(this));
    nodeCount++;
  }
}

void CFLOBDDLeafNode::CountNodes(Hashset<CFLOBDDNodeHandle> *visitedNodes, unsigned int &nodeCount)
{
	if (visitedNodes->Lookup(new CFLOBDDNodeHandle(this)) == NULL) {
		visitedNodes->Insert(new CFLOBDDNodeHandle(this));
		nodeCount++;
	}
}

void CFLOBDDLeafNode::CountPaths(Hashset<CFLOBDDNodeHandle> *visitedNodes)
{
	CFLOBDDNodeHandle* handle = new CFLOBDDNodeHandle(this);
	if (visitedNodes->Lookup(handle) == NULL) {
		visitedNodes->Insert(handle);
	}
}

void CFLOBDDLeafNode::IncrRef() { }
void CFLOBDDLeafNode::DecrRef() { }

//********************************************************************
// CFLOBDDForkNode
//********************************************************************

// Constructors/Destructor -------------------------------------------

// Default constructor
CFLOBDDForkNode::CFLOBDDForkNode()
  :  CFLOBDDLeafNode()
{
  numExits = 2;
//#ifdef PATH_COUNTING_ENABLED
  //numPathsToExit = new unsigned long long int[2];
  //numPathsToExit = new cpp_int[2];
  numPathsToExit = new long double[2];
  numPathsToExit[0] = 0;
  numPathsToExit[1] = 0;
//#endif
}

CFLOBDDForkNode::~CFLOBDDForkNode()
{
}

// print
std::ostream& CFLOBDDForkNode::print(std::ostream & out) const
{
  for (unsigned int i = level; i < maxLevel; i++) {
    out << "  ";
  }
  out << "Fork";
  return out;
}

bool CFLOBDDForkNode::IsValid()
{
	if (level == 0)
		isValid = true;
	return level == 0;
}

void CFLOBDDForkNode::FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index)
{
  assert(i <= 1);
  index--;
  assignment[index] = i;
}

int CFLOBDDForkNode::Traverse(SH_OBDD::AssignmentIterator &ai)
{
  bool val = ai.Current();
  ai.Next();
  return (int)val;
}

CFLOBDDNodeHandle CFLOBDDForkNode::Reduce(ReductionMapHandle&, unsigned int replacementNumExits, bool forceReduce)
{
	if (forceReduce){
		if (replacementNumExits == 1) {
			return CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle;
		}
		else {
			//assert(replacementNumExits == 2);
			return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
		}
	}
	else{
		assert(replacementNumExits == 2);
		return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
	}
}

unsigned int CFLOBDDForkNode::Hash(unsigned int modsize)
{
  return ((unsigned int)reinterpret_cast<uintptr_t>(this) >> 2) % modsize;
}

// Overloaded !=
bool CFLOBDDForkNode::operator!= (const CFLOBDDNode & n)
{
  return n.NodeKind() != CFLOBDD_FORK;
}

// Overloaded ==
bool CFLOBDDForkNode::operator== (const CFLOBDDNode & n)
{
  return n.NodeKind() == CFLOBDD_FORK;
}

//********************************************************************
// CFLOBDDDontCareNode
//********************************************************************

// Constructors/Destructor -------------------------------------------

// Default constructor
CFLOBDDDontCareNode::CFLOBDDDontCareNode()
  :  CFLOBDDLeafNode()
{
  numExits = 1;
//#ifdef PATH_COUNTING_ENABLED
  //numPathsToExit = new unsigned long long int[1];
  //numPathsToExit = new cpp_int[1];
  numPathsToExit = new long double[1];
  numPathsToExit[0] = 1;
//#endif
}

CFLOBDDDontCareNode::~CFLOBDDDontCareNode()
{
}

// print
std::ostream& CFLOBDDDontCareNode::print(std::ostream & out) const
{
  for (unsigned int i = level; i < maxLevel; i++) {
    out << "  ";
  }
  out << "Don't care";
  return out;
}

bool CFLOBDDDontCareNode::IsValid()
{
	if (level == 0)
		isValid = true;
	return (level == 0);
}

void CFLOBDDDontCareNode::FillSatisfyingAssignment(unsigned int, SH_OBDD::Assignment &assignment, unsigned int &index)
{
  index--;
  assignment[index] = 0;
}

int CFLOBDDDontCareNode::Traverse(SH_OBDD::AssignmentIterator &ai)
{
  ai.Next();
  return 0;
}

CFLOBDDNodeHandle CFLOBDDDontCareNode::Reduce(ReductionMapHandle&, unsigned int, bool)
{
  return CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle;
}

unsigned int CFLOBDDDontCareNode::Hash(unsigned int modsize)
{
  return ((unsigned int) reinterpret_cast<uintptr_t>(this) >> 2) % modsize;
}

// Overloaded !=
bool CFLOBDDDontCareNode::operator!= (const CFLOBDDNode & n)
{
  return n.NodeKind() != CFLOBDD_DONTCARE;
}

// Overloaded ==
bool CFLOBDDDontCareNode::operator== (const CFLOBDDNode & n)
{
  return n.NodeKind() == CFLOBDD_DONTCARE;
}

// Restrict -----------------------------------------------------------
namespace CFL_OBDD {
CFLOBDDNodeHandle Restrict(CFLOBDDNodeHandle g, unsigned int i, bool val,
                           CFLOBDDReturnMapHandle &MapHandle
                          )
{
    CFLOBDDNodeHandle answer;
  
    if (g.handleContents->NodeKind() == CFLOBDD_INTERNAL) {
      answer = Restrict((CFLOBDDInternalNode *)g.handleContents, i, val,
                        MapHandle
                       );
    }
    else if (g.handleContents->NodeKind() == CFLOBDD_FORK) {
      if (val == false) {
        MapHandle.AddToEnd(0);
        MapHandle.Canonicalize();
        answer = CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle;
      }
      else { /* val == true */
        MapHandle.AddToEnd(1);
        MapHandle.Canonicalize();
        answer = CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle;
      }
    }
    else { /* g.handleContents->NodeKind() == CFLOBDD_DONTCARE */
      MapHandle.AddToEnd(0);
      MapHandle.Canonicalize();
      answer = g;
    }
    return answer;
}

CFLOBDDNodeHandle Restrict(CFLOBDDInternalNode *g, unsigned int i, bool val,
                           CFLOBDDReturnMapHandle &MapHandle
                          )
{
  if (g == CFLOBDDNodeHandle::NoDistinctionNode[g->level].handleContents) {
    MapHandle.AddToEnd(0);
    MapHandle.Canonicalize();
    return CFLOBDDNodeHandle(g);
  }

  CFLOBDDReturnMapHandle AMap;
  CFLOBDDInternalNode *n = new CFLOBDDInternalNode(g->level);
  unsigned int j;
  unsigned int curExit;
  int b;

  if (i < (unsigned int)(1 << (g->level-1))) { // i falls in AConnection range
  	CFLOBDDNodeHandle aHandle = Restrict(*(g->AConnection.entryPointHandle), i, val, AMap);
    n->AConnection.entryPointHandle = &aHandle;
    for (unsigned int k = 0; k < AMap.Size(); k++) {
      n->AConnection.returnMapHandle.AddToEnd(k);
    }
    n->AConnection.returnMapHandle.Canonicalize();
    j = 0;
    curExit = 0;
    n->numBConnections = AMap.Size();
    n->BConnection = new Connection[n->numBConnections];
	unsigned AMapSize = AMap.mapContents->mapArray.size();
    for (unsigned sAI = 0; sAI < AMapSize; sAI++)
	{
      b = AMap.mapContents->mapArray[sAI];
      n->BConnection[j].entryPointHandle = g->BConnection[b].entryPointHandle;
      // Fill in n->BConnection[j].returnMapHandle and add new items (as appropriate) to MapHandle
      CFLOBDDReturnMapHandle BMap = g->BConnection[b].returnMapHandle;
	  unsigned BMapSize = BMap.mapContents->mapArray.size();
	  for (unsigned sBI = 0; sBI < BMapSize; sBI++)
	  {
		  int c = BMap.mapContents->mapArray[sBI];
        // Test whether c occurs in MapHandle
           if (MapHandle.Member(c)) {
             int index = MapHandle.LookupInv(c);
             n->BConnection[j].returnMapHandle.AddToEnd(index);
           }
           else {   // New pair found (i.e., new exit node found)
             MapHandle.AddToEnd(c);
             n->BConnection[j].returnMapHandle.AddToEnd(curExit);
             curExit++;
           }
      }
      n->BConnection[j].returnMapHandle.Canonicalize();
      j++;
    }
  }

  else { // i falls in BConnection range
    ReductionMapHandle AReductionMapHandle;  // To record pattern of duplicate BConnections

    curExit = 0;
    n->BConnection = new Connection[g->numBConnections];   // May create shorter version later
    n->numBConnections = 0;
    for (j = 0; j < g->numBConnections; j++) { // Perform a Restrict for each middle vertex
      CFLOBDDReturnMapHandle BMap;
      CFLOBDDNodeHandle m = Restrict(*(g->BConnection[j].entryPointHandle), i-(unsigned int)(1 << (g->level-1)), val, BMap);

      // Fill in inducedReturnMapHandleB and add new items (as appropriate) to MapHandle
      CFLOBDDReturnMapHandle inducedReturnMapHandleB;
      unsigned BMapSize = BMap.mapContents->mapArray.size();
	  for (unsigned sBI = 0; sBI < BMapSize; sBI++)
	  {
		  int c = g->BConnection[j].returnMapHandle.Lookup(sBI);
        // Test whether c occurs in MapHandle
           if (MapHandle.Member(c)) {
             int index = MapHandle.LookupInv(c);
             inducedReturnMapHandleB.AddToEnd(index);
           }
           else {   // New pair found (i.e., new exit node found)
             MapHandle.AddToEnd(c);
             inducedReturnMapHandleB.AddToEnd(curExit);
             curExit++;
           }
      }
      inducedReturnMapHandleB.Canonicalize();
      Connection candidate(m, inducedReturnMapHandleB);
      unsigned int position = n->InsertBConnection(n->numBConnections, candidate);
      AReductionMapHandle.AddToEnd(position);
    }
    AReductionMapHandle.Canonicalize();

    if (n->numBConnections < g->numBConnections) {  // Shorten
      Connection *temp = n->BConnection;
      n->BConnection = new Connection[n->numBConnections];
      for (unsigned int j = 0; j < n->numBConnections; j++) {
        n->BConnection[j] = temp[j];
      }
      delete [] temp;
    }

    // Reduce the A-connection w.r.t AReductionMapHandle
       ReductionMapHandle inducedAReductionMapHandle;
       CFLOBDDReturnMapHandle inducedAReturnMap;
       CFLOBDDReturnMapHandle reducedAReturnMap = g->AConnection.returnMapHandle.Compose(AReductionMapHandle);
       reducedAReturnMap.InducedReductionAndReturnMap(inducedAReductionMapHandle, inducedAReturnMap);
       CFLOBDDNodeHandle tempHandle = g->AConnection.entryPointHandle->Reduce(inducedAReductionMapHandle, inducedAReturnMap.Size());
       n->AConnection = Connection(tempHandle, inducedAReturnMap);
  }

  n->numExits = curExit;
#ifdef PATH_COUNTING_ENABLED
  n->InstallPathCounts();
#endif
  MapHandle.Canonicalize();
  return CFLOBDDNodeHandle(n);
}

}
