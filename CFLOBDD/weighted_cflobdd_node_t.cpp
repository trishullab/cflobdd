#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <cstdarg>
#include <unordered_set>
#include <map>
#include <unordered_map>
//#include <mpirxx.h>

#include "ntz_T.h"
#include "weighted_cflobdd_node_t.h"
// #include "weighted_connectionT.h"
#include "list_T.h"
#include "list_TPtr.h"
#include "intpair.h"
#include "conscell.h"
#include "assignment.h"
#include "bool_op.h"
#include "return_map_T.h"
#include "reduction_map.h"
#include "hash.h"
#include "hashset.h"
#include "weighted_values.h"

using namespace CFL_OBDD;

//********************************************************************
// WeightedCFLOBDDNodeHandleT
//
// Contains a canonical WeightedCFLOBDDNode*
//********************************************************************

// Initializations of static members ---------------------------------

template <typename T, typename Op>
Hashset<WeightedCFLOBDDNode<T, Op>> *WeightedCFLOBDDNodeHandleT<T, Op>::canonicalNodeTable = new Hashset<WeightedCFLOBDDNode<T,Op>>(HASHSET_NUM_BUCKETS);
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op> *WeightedCFLOBDDNodeHandleT<T, Op>::NoDistinctionNode = NULL;
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op> *WeightedCFLOBDDNodeHandleT<T, Op>::NoDistinctionNode_Ann = NULL;
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op> *WeightedCFLOBDDNodeHandleT<T, Op>::IdentityNode = NULL;
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op> WeightedCFLOBDDNodeHandleT<T, Op>::CFLOBDDForkNodeHandle;
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op> WeightedCFLOBDDNodeHandleT<T, Op>::CFLOBDDForkNodeHandle01;
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op> WeightedCFLOBDDNodeHandleT<T, Op>::CFLOBDDForkNodeHandle10;
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op> WeightedCFLOBDDNodeHandleT<T, Op>::CFLOBDDDontCareNodeHandle;
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op> *WeightedCFLOBDDNodeHandleT<T, Op>::AdditionInterleavedNoCarryNode = NULL;
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op> *WeightedCFLOBDDNodeHandleT<T, Op>::AdditionInterleavedCarryNode = NULL;

// InitNoDistinctionTable
//
// Create and record the "NoDistinctionNodes" (one per level).
// As a side effect, these nodes are entered in the canonicalNodeTable.
//
template <typename T, typename Op>
void WeightedCFLOBDDNodeHandleT<T, Op>::InitNoDistinctionTable()
{
	// TODO: Place this in a better place
  NoDistinctionNode = new WeightedCFLOBDDNodeHandleT<T,Op>[CFLOBDDMaxlevel+1];

  for (unsigned int i = 0; i <= CFLOBDDMaxlevel; i++) {
    if (i == 0) {
      NoDistinctionNode[0] = WeightedCFLOBDDNodeHandleT<T, Op>(new WeightedCFLOBDDDontCareNode<T, Op>());
    }
    else {
      WeightedCFLOBDDInternalNode<T,Op> *n;
      CFLOBDDReturnMapHandle m1, m2;

      n = new WeightedCFLOBDDInternalNode<T,Op>(i);
      n->AConnection.entryPointHandle = &(NoDistinctionNode[i-1]);
      m1.AddToEnd(0);
      m1.Canonicalize();
	  n->AConnection.returnMapHandle = m1;//m1
  
      n->numBConnections = 1;
      n->BConnection = new WConnection<T,Op>[1];
      n->BConnection[0].entryPointHandle = &(NoDistinctionNode[i-1]);
      m2.AddToEnd(0);
      m2.Canonicalize();
	  n->BConnection[0].returnMapHandle = m2;
      n->numExits = 1;
//#ifdef PATH_COUNTING_ENABLED
      n->InstallPathCounts();
//#endif
      NoDistinctionNode[i] = WeightedCFLOBDDNodeHandleT(n);
    }
  }

  CFLOBDDForkNodeHandle = WeightedCFLOBDDNodeHandleT(new WeightedCFLOBDDForkNode<T,Op>());
  CFLOBDDForkNodeHandle01 = WeightedCFLOBDDNodeHandleT(new WeightedCFLOBDDForkNode<T,Op>(getAnnhilatorValue<T,Op>(), getIdentityValue<T,Op>()));
  CFLOBDDForkNodeHandle10 = WeightedCFLOBDDNodeHandleT(new WeightedCFLOBDDForkNode<T,Op>(getIdentityValue<T,Op>(), getAnnhilatorValue<T,Op>()));
  CFLOBDDDontCareNodeHandle = NoDistinctionNode[0];
  return;
} // InitNoDistinctionTable

template <typename T, typename Op>
void WeightedCFLOBDDNodeHandleT<T, Op>::InitNoDistinctionTable_Ann()
{
	// TODO: Place this in a better place
  NoDistinctionNode_Ann = new WeightedCFLOBDDNodeHandleT<T,Op>[CFLOBDDMaxlevel+1];

  for (unsigned int i = 0; i <= CFLOBDDMaxlevel; i++) {
    if (i == 0) {
      NoDistinctionNode_Ann[0] = WeightedCFLOBDDNodeHandleT<T, Op>(new WeightedCFLOBDDDontCareNode<T, Op>(false));
    }
    else {
      WeightedCFLOBDDInternalNode<T,Op> *n;
      CFLOBDDReturnMapHandle m1, m2;

      n = new WeightedCFLOBDDInternalNode<T,Op>(i);
      n->AConnection.entryPointHandle = &(NoDistinctionNode_Ann[i-1]);
      m1.AddToEnd(0);
      m1.Canonicalize();
	  n->AConnection.returnMapHandle = m1;//m1
  
      n->numBConnections = 1;
      n->BConnection = new WConnection<T,Op>[1];
      n->BConnection[0].entryPointHandle = &(NoDistinctionNode_Ann[i-1]);
      m2.AddToEnd(0);
      m2.Canonicalize();
	  n->BConnection[0].returnMapHandle = m2;
      n->numExits = 1;
//#ifdef PATH_COUNTING_ENABLED
      n->InstallPathCounts();
//#endif
      NoDistinctionNode_Ann[i] = WeightedCFLOBDDNodeHandleT(n);
    }
  }
  return;
} // InitNoDistinctionTable

template <typename T, typename Op>
void WeightedCFLOBDDNodeHandleT<T, Op>::InitIdentityNodeTable()
{
	// TODO: Place this in a better place
  IdentityNode = new WeightedCFLOBDDNodeHandleT<T,Op>[CFLOBDDMaxlevel+1];

  for (unsigned int i = 1; i <= CFLOBDDMaxlevel; i++) {
    if (i == 1) {
        WeightedCFLOBDDInternalNode<T,Op> *n;
        n = new WeightedCFLOBDDInternalNode<T,Op>(i);
        CFLOBDDReturnMapHandle m01, m10;
        m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
        m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
        n->AConnection = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle, m01);
        n->numBConnections = 2;
        n->BConnection = new WConnection<T,Op>[n->numBConnections];
        WeightedCFLOBDDNodeHandleT<T,Op> f10, f01;
        n->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle10, m01);
        n->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle01, m10);
        n->numExits = 2;
        IdentityNode[1] = WeightedCFLOBDDNodeHandleT<T, Op>(n);
    }
    else {
      WeightedCFLOBDDInternalNode<T,Op> *n;
      CFLOBDDReturnMapHandle m1, m2;

      n = new WeightedCFLOBDDInternalNode<T,Op>(i);
      n->AConnection.entryPointHandle = &(IdentityNode[i-1]);
      m1.AddToEnd(0);
      m1.AddToEnd(1);
      m1.Canonicalize();
	  n->AConnection.returnMapHandle = m1;//m1
  
      n->numBConnections = 2;
      n->BConnection = new WConnection<T,Op>[2];
      n->BConnection[0].entryPointHandle = &(IdentityNode[i-1]);
      m2.AddToEnd(1);
      m2.Canonicalize();
	  n->BConnection[0].returnMapHandle = m1;
      n->BConnection[1] = WConnection<T,Op>(NoDistinctionNode_Ann[i-1], m2);
      n->numExits = 2;
//#ifdef PATH_COUNTING_ENABLED
    //   n->InstallPathCounts();
//#endif
      IdentityNode[i] = WeightedCFLOBDDNodeHandleT(n);
    }
  }
} // InitNoDistinctionTable



// InitAdditionInterleavedTable
//
// Create and record the entries for "AdditionInterleavedNoCarryNode" and "AdditionInterleavedCarryNode"(each have one entry per level).
// As a side effect, these nodes are entered in the canonicalNodeTable.
//
// Note that the tables have CFLOBDDMaxlevel-1 entries, corresponding to levels 2, ..., CFLOBDDMaxlevel
//
// TODO: Check this
template <typename T, typename Op>
void WeightedCFLOBDDNodeHandleT<T,Op>::InitAdditionInterleavedTable()
{
	AdditionInterleavedNoCarryNode = new WeightedCFLOBDDNodeHandleT[CFLOBDDMaxlevel - 1];
	AdditionInterleavedCarryNode   = new WeightedCFLOBDDNodeHandleT[CFLOBDDMaxlevel - 1];

	WeightedCFLOBDDInternalNode<T,Op> *n_no_carry;
	WeightedCFLOBDDInternalNode<T,Op> *n_carry;

	CFLOBDDReturnMapHandle m0, m1, m01, m10;
	CFLOBDDReturnMapHandle m02, m12, m21, m20, m012, m102;
	m0.AddToEnd(0);
	m0.Canonicalize();

	m1.AddToEnd(1);
	m1.Canonicalize();

	m01.AddToEnd(0);
	m01.AddToEnd(1);
	m01.Canonicalize();
	
	m02.AddToEnd(0);
	m02.AddToEnd(2);
	m02.Canonicalize();
	
	m10.AddToEnd(1);
	m10.AddToEnd(0);
	m10.Canonicalize();
	
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

	for (unsigned int i = 0; i < CFLOBDDMaxlevel-1; i++) {  // Note: i+2 is the node level
		unsigned int level = i + 2;
		if (i == 0) {
			WeightedCFLOBDDInternalNode<T,Op> *n1 = new WeightedCFLOBDDInternalNode<T,Op>(level-1);
			n1->AConnection = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m01);//m01
			n1->numBConnections = 2;
			n1->BConnection = new WConnection<T,Op>[n1->numBConnections];
			n1->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m01);//m01
			n1->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m12);
			n1->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			n1->InstallPathCounts();
#endif
			WeightedCFLOBDDNodeHandleT nh1(n1);

			WeightedCFLOBDDInternalNode<T,Op> *n2 = new WeightedCFLOBDDInternalNode<T,Op>(level-1);
			n2->AConnection = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m01);//m01
			n2->numBConnections = 2;
			n2->BConnection = new WConnection<T,Op>[n2->numBConnections];
			n2->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDDontCareNodeHandle, m0);//m0
			n2->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDDontCareNodeHandle, m1);//m1
			n2->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
			n2->InstallPathCounts();
#endif
			WeightedCFLOBDDNodeHandleT nh2(n2);

			// Create the no-carry node
			n_no_carry = new WeightedCFLOBDDInternalNode<T,Op>(level);
			n_no_carry->AConnection = WConnection<T,Op>(nh1, m012);
			n_no_carry->numBConnections = 3;
			n_no_carry->BConnection = new WConnection<T,Op>[n_no_carry->numBConnections];
			n_no_carry->BConnection[0] = WConnection<T,Op>(nh2, m01);//01
			n_no_carry->BConnection[1] = WConnection<T,Op>(nh2, m10);//10
			n_no_carry->BConnection[2] = WConnection<T,Op>(nh2, m21);
			n_no_carry->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			n_no_carry->InstallPathCounts();
#endif
			AdditionInterleavedNoCarryNode[0] = WeightedCFLOBDDNodeHandleT(n_no_carry);

			// Create the carry node
			n_carry = new WeightedCFLOBDDInternalNode<T,Op>(level);
			n_carry->AConnection = WConnection<T,Op>(nh1, m012);
			n_carry->numBConnections = 3;
			n_carry->BConnection = new WConnection<T,Op>[n_carry->numBConnections];
			n_carry->BConnection[0] = WConnection<T,Op>(nh2, m01);//m01
			n_carry->BConnection[1] = WConnection<T,Op>(nh2, m20);
			n_carry->BConnection[2] = WConnection<T,Op>(nh2, m02);
			n_carry->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			n_carry->InstallPathCounts();
#endif
			AdditionInterleavedCarryNode[0] = WeightedCFLOBDDNodeHandleT(n_carry);
		} // if (i == 0)
		else {  // i > 0
			WeightedCFLOBDDNodeHandleT c0 = AdditionInterleavedNoCarryNode[i-1];
			WeightedCFLOBDDNodeHandleT c1 = AdditionInterleavedCarryNode[i-1];

			// Create the no-carry node
			n_no_carry = new WeightedCFLOBDDInternalNode<T,Op>(level);
			n_no_carry->AConnection = WConnection<T,Op>(c0, m012);
			n_no_carry->numBConnections = 3;
			n_no_carry->BConnection = new WConnection<T,Op>[n_no_carry->numBConnections];
			n_no_carry->BConnection[0] = WConnection<T,Op>(c0, m012);
			n_no_carry->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m1);//m1
			n_no_carry->BConnection[2] = WConnection<T,Op>(c1, m102);
			n_no_carry->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			n_no_carry->InstallPathCounts();
#endif
			AdditionInterleavedNoCarryNode[i] = WeightedCFLOBDDNodeHandleT(n_no_carry);

			// Create the carry node
			n_carry = new WeightedCFLOBDDInternalNode<T,Op>(level);
			n_carry->AConnection = WConnection<T,Op>(c1, m012);
			n_carry->numBConnections = 3;
			n_carry->BConnection = new WConnection<T,Op>[n_carry->numBConnections];
			n_carry->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m0);//m0
			n_carry->BConnection[1] = WConnection<T,Op>(c0, m102);
			n_carry->BConnection[2] = WConnection<T,Op>(c1, m012);
			n_carry->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			n_carry->InstallPathCounts();
#endif
			AdditionInterleavedCarryNode[i] = WeightedCFLOBDDNodeHandleT(n_carry);
		}
	}
	return;
} // InitAdditionInterleavedTable

// Constructors/Destructor -------------------------------------------

// Default constructor
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op>::WeightedCFLOBDDNodeHandleT()
  : handleContents(NULL)
{
}

// Constructor
//
// Construct and canonicalize
//
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T, Op>::WeightedCFLOBDDNodeHandleT(WeightedCFLOBDDNode<T,Op> *n)
  : handleContents(n)
{
  assert(n != NULL);
  handleContents->IncrRef();
  Canonicalize();
}

// Copy constructor
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T,Op>::WeightedCFLOBDDNodeHandleT(const WeightedCFLOBDDNodeHandleT<T,Op> &c)
{
  handleContents = c.handleContents;
  if (handleContents != NULL) {
	// std::cout << handleContents << std::endl;
    handleContents->IncrRef();
	// std::cout << "handleContents" << std::endl;
  }
}

template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T,Op>::~WeightedCFLOBDDNodeHandleT()
{
  if (handleContents != NULL) {
    handleContents->DecrRef();
  }
}

// Hash
template <typename T, typename Op>
unsigned int WeightedCFLOBDDNodeHandleT<T,Op>::Hash(unsigned int modsize)
{
  return ((unsigned int) reinterpret_cast<uintptr_t>(handleContents) >> 2) % modsize;
}

// Overloaded !=
template <typename T, typename Op>
bool WeightedCFLOBDDNodeHandleT<T,Op>::operator!= (const WeightedCFLOBDDNodeHandleT<T,Op> & C)
{
  return handleContents != C.handleContents;
}

// Overloaded ==
template <typename T, typename Op>
bool WeightedCFLOBDDNodeHandleT<T,Op>::operator== (const WeightedCFLOBDDNodeHandleT & C) const
{
  return handleContents == C.handleContents;
}

// Overloaded assignment
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T,Op> & WeightedCFLOBDDNodeHandleT<T,Op>::operator= (const WeightedCFLOBDDNodeHandleT<T,Op> &c)
{
  if (this != &c)      // don't assign to self!
  {
    WeightedCFLOBDDNode<T,Op> *temp = handleContents;
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

template <typename T, typename Op>
class ReducePair{
	public:
		WeightedCFLOBDDNodeHandleT<T,Op> m;
		T value;

        ReducePair() {}

		ReducePair(WeightedCFLOBDDNodeHandleT<T,Op> p, T v)
		{
			m = p;
			value = v;
		}

		struct ReducePairHash {
			size_t operator()(const ReducePair& p) const
			{
                boost::hash<T> boost_hash;
				auto hash1 = p.m.Hash(997);
				auto hash2 = boost_hash(p.value);
				return 117 * (hash1 + 1) + hash2;
			}
		};

		bool operator==(const ReducePair& p) const
		{
			return (m == p.m) && (value == p.value);
		}
};

template <typename T, typename Op> 
static Hashtable<WeightedCFLReduceKey<T, Op>, ReducePair<T,Op>> *reduceCache = NULL;

// TODO: check this
template <typename T, typename Op> 
std::pair<WeightedCFLOBDDNodeHandleT<T,Op>, T> WeightedCFLOBDDNodeHandleT<T,Op>::Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valueTuple, bool addCheck)
{
	if (redMapHandle.mapContents->isIdentityMap && valueTuple.mapContents->isAllSame) {
        if (valueTuple.mapContents->value == getAnnhilatorValue<T,Op>())
            return std::make_pair(WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[this->handleContents->level], getAnnhilatorValue<T,Op>());
		return std::make_pair(*this, valueTuple.mapContents->value);
	}

    if (addCheck && redMapHandle.mapContents->isIdentityMap && valueTuple.mapContents->isOneOrZero && this->handleContents->level != 0)
    {
        return std::make_pair(*this, getIdentityValue<T,Op>());
    }

	ReducePair<T,Op> cachedNodeHandle;
	bool isCached = reduceCache<T,Op>->Fetch(WeightedCFLReduceKey<T,Op>(*this, redMapHandle, valueTuple), cachedNodeHandle);
	if (isCached) {
		// std::cout << "Hit : " << handleContents->Level() << std::endl;
		return std::make_pair(cachedNodeHandle.m, cachedNodeHandle.value);
	}
	else {
		// std::cout << "Miss: " << handleContents->Level() << std::endl;
		std::pair<WeightedCFLOBDDNodeHandleT<T,Op>, T> temp;
		temp = handleContents->Reduce(redMapHandle, replacementNumExits, valueTuple, addCheck);
		reduceCache<T,Op>->Insert(WeightedCFLReduceKey(*this, redMapHandle, valueTuple), ReducePair<T,Op>(temp.first, temp.second));
		return temp;
	}
}

template <typename T, typename Op> 
void WeightedCFLOBDDNodeHandleT<T,Op>::InitReduceCache()
{
  reduceCache<T,Op> = new Hashtable<WeightedCFLReduceKey<T,Op>, ReducePair<T,Op>>(HASH_NUM_BUCKETS);
}

template <typename T, typename Op> 
void WeightedCFLOBDDNodeHandleT<T,Op>::DisposeOfReduceCache()
{
	delete reduceCache<T,Op>;
	reduceCache<T,Op> = NULL;
}

// Canonicalization --------------------------------------------
template <typename T, typename Op> 
void WeightedCFLOBDDNodeHandleT<T,Op>::Canonicalize()
{
  WeightedCFLOBDDNode<T,Op> *answerContents;

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
template <typename T, typename Op> 
std::ostream& WeightedCFLOBDDNodeHandleT<T,Op>::print(std::ostream & out) const
{
  out << *handleContents << std::endl;
  return out;
}

// template <typename T, typename Op> 
// bool WeightedCFLOBDDNodeHandleT<T,Op>::IsValid()
// {
// 	return handleContents->IsValid();
// }

namespace CFL_OBDD {

template <typename T, typename Op> 
std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDNodeHandleT<T,Op> &d)
{
  d.print(out);
  return(out);
}

/*
template <typename T, typename Op> 
WeightedCFLOBDDNodeHandleT<T,Op> MkDistinction(unsigned int level, unsigned int i)
{
  if (level == 0) {
    assert(i == 0);
    return WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle;
  }
  else {  // Create an appropriate WeightedCFLOBDDInternalNode
    WeightedCFLOBDDInternalNode<T,Op> *n = new WeightedCFLOBDDInternalNode<T,Op>(level);
    CFLOBDDReturnMapHandle m1, m2, m3;
    if (i < (unsigned int)(1 << (level-1))) { // i falls in AConnection range
      m1.AddToEnd(0);
      m1.AddToEnd(1);
      m1.Canonicalize();
      WeightedCFLOBDDNodeHandleT temp;
	  temp = MkDistinction(level-1, i);
      n->AConnection = WConnection<T,Op>(temp, m1);

      n->numBConnections = 2;
      n->BConnection = new WConnection<T,Op>[n->numBConnections];
      m2.AddToEnd(0);
      m2.Canonicalize();
      n->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m2);
      m3.AddToEnd(1);
      m3.Canonicalize();
      n->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m3);
      n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
      n->InstallPathCounts();
#endif
    }
    else {         // i falls in BConnection range
      m1.AddToEnd(0);
      m1.Canonicalize();
      n->AConnection = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m1);

      i = i ^ (1 << (level-1));  // Mask off high-order bit for recursive call
      n->numBConnections = 1;
      n->BConnection = new WConnection<T,Op>[n->numBConnections];
      m2.AddToEnd(0);
      m2.AddToEnd(1);
      m2.Canonicalize();
      WeightedCFLOBDDNodeHandleT temp = MkDistinction(level-1, i);
      n->BConnection[0] = WConnection<T,Op>(temp, m2);
      n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
      n->InstallPathCounts();
#endif
    }
    return WeightedCFLOBDDNodeHandleT(n);
  }
} // MkDistinction


WeightedCFLOBDDNodeHandleT MkAdditionInterleavedRecursive(unsigned int level, bool carry)
{
  WeightedCFLOBDDInternalNode<T,Op> *n;

  //CFLOBDDReturnMapHandle m0, m1, m01, m10;
  CFLOBDDReturnMapHandle m02, m12, m21, m20, m012, m102;
  m0.AddToEnd(0);
  m0.Canonicalize();

  m1.AddToEnd(1);
  m1.Canonicalize();

  m01.AddToEnd(0);
  m01.AddToEnd(1);
  m01.Canonicalize();

  m02.AddToEnd(0);
  m02.AddToEnd(2);
  m02.Canonicalize();

  m10.AddToEnd(1);
  m10.AddToEnd(0);
  m10.Canonicalize();

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
    WeightedCFLOBDDInternalNode<T,Op> *n1 = new WeightedCFLOBDDInternalNode<T,Op>(1);
    n1->AConnection = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m01);//m01
    n1->numBConnections = 2;
    n1->BConnection = new WConnection<T,Op>[n1->numBConnections];
    n1->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m01);//m01
    n1->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m12);
    n1->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
    n1->InstallPathCounts();
#endif
    WeightedCFLOBDDNodeHandleT nh1(n1);

    WeightedCFLOBDDInternalNode<T,Op> *n2 = new WeightedCFLOBDDInternalNode<T,Op>(1);
    n2->AConnection =  WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m01);//m01
    n2->numBConnections = 2;
    n2->BConnection = new WConnection<T,Op>[n2->numBConnections];
    n2->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDDontCareNodeHandle, m01);//m0
    n2->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDDontCareNodeHandle, m01);//m1
    n2->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
    n2->InstallPathCounts();
#endif
    WeightedCFLOBDDNodeHandleT nh2(n2);

    n = new WeightedCFLOBDDInternalNode<T,Op>(2);
    n->AConnection = WConnection<T,Op>(nh1, m012);
    n->numBConnections = 3;
    n->BConnection = new WConnection<T,Op>[n->numBConnections];
    if (!carry) {
      n->BConnection[0] = WConnection<T,Op>(nh2, m01);//m01
      n->BConnection[1] = WConnection<T,Op>(nh2, m01);//m10
      n->BConnection[2] = WConnection<T,Op>(nh2, m21);
    }
    else {
      assert(carry == 1);
      n->BConnection[0] = WConnection<T,Op>(nh2, m01);//m01
      n->BConnection[1] = WConnection<T,Op>(nh2, m20);
      n->BConnection[2] = WConnection<T,Op>(nh2, m02);
    }
    n->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
    n->InstallPathCounts();
#endif
  }
  else {  // level > 2; Create an appropriate CFLOBDDInternalNode
    WeightedCFLOBDDNodeHandleT c0 = MkAdditionInterleavedRecursive(level-1, false);
    WeightedCFLOBDDNodeHandleT c1 = MkAdditionInterleavedRecursive(level-1, true);

    n = new WeightedCFLOBDDInternalNode<T,Op>(level);
    n->AConnection = WConnection<T,Op>((carry ? c1 : c0), m012);
    n->numBConnections = 3;
    n->BConnection = new WConnection<T,Op>[n->numBConnections];
    if (!carry) {
      n->BConnection[0] = WConnection<T,Op>(c0, m012);
      n->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m1);//m1
      n->BConnection[2] = WConnection<T,Op>(c1, m102);
    }
    else {
      assert(carry == 1);
      n->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m0);//m0
      n->BConnection[1] = WConnection<T,Op>(c0, m102);
      n->BConnection[2] = WConnection<T,Op>(c1, m012);
    }
    n->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
    n->InstallPathCounts();
#endif
  }
  return WeightedCFLOBDDNodeHandleT(n);
} // MkAdditionInterleavedRecursive

WeightedCFLOBDDNodeHandleT MkParity(unsigned int level)
{
  if (level == 0) {
    return WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle;
  }
  else {  // Create an appropriate WeightedCFLOBDDInternalNode<T,Op>
    WeightedCFLOBDDInternalNode<T,Op> *n = new WeightedCFLOBDDInternalNode<T,Op>(level);
    CFLOBDDReturnMapHandle m1, m2;
    m1.AddToEnd(0);
    m1.AddToEnd(1);
    m1.Canonicalize();
    WeightedCFLOBDDNodeHandleT temp = MkParity(level-1);
    n->AConnection = WConnection<T,Op>(temp, m01);//m01

    n->numBConnections = 2;
    n->BConnection = new WConnection<T,Op>[n->numBConnections];
    n->BConnection[0] = WConnection<T,Op>(temp, m01);//m01
    //m2.AddToEnd(1);
    //m2.AddToEnd(0);
    //m2.Canonicalize();
    n->BConnection[1] = WConnection<T,Op>(temp, m10);//m10

    n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
    n->InstallPathCounts();
#endif
    return WeightedCFLOBDDNodeHandleT(n);
  }
} // MkParity

WeightedCFLOBDDNodeHandleT MkStepOneFourth(unsigned int level)
{
  assert(level >= 1);

  WeightedCFLOBDDInternalNode<T,Op> *n = new WeightedCFLOBDDInternalNode<T,Op>(level);
  
  if (level == 1) {
    CFLOBDDReturnMapHandle m1, m2;
    m1.AddToEnd(0);
    m1.AddToEnd(1);
    m1.Canonicalize();
    n->AConnection = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m1);

    n->numBConnections = 2;
    n->BConnection = new WConnection<T,Op>[n->numBConnections];
    n->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m1);
    m2.AddToEnd(1);
    m2.Canonicalize();
    n->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDDontCareNodeHandle, m2);
  }
  else {  // Create an appropriate WeightedCFLOBDDInternalNode
    CFLOBDDReturnMapHandle m1, m2, m3;
    m1.AddToEnd(0);
    m1.AddToEnd(1);
    m1.Canonicalize();
    WeightedCFLOBDDNodeHandleT temp = MkStepOneFourth(level-1);
    n->AConnection = WConnection<T,Op>(temp, m1);

    n->numBConnections = 2;
    n->BConnection = new WConnection<T,Op>[n->numBConnections];
    m2.AddToEnd(0);
    m2.Canonicalize();
    n->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m2);
    m3.AddToEnd(1);
    m3.Canonicalize();
    n->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m3);
  }
  n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
  n->InstallPathCounts();
#endif
  return WeightedCFLOBDDNodeHandleT(n);
} // MkStepOneFourth

double ComputeProbabilityNode(WeightedCFLOBDDNodeHandleT g, std::vector<double>& var_probs, std::vector<double>& path_probs, int start, int end){
	if (g.handleContents->level == 0){
		if (g == WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle){
				return (1 - var_probs[start]) * path_probs[0] + var_probs[start] * path_probs[1];
		} else {
			return path_probs[0];
		}
	} else if (g == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g.handleContents->level]){
			return path_probs[0];
	} else {
		WeightedCFLOBDDInternalNode<T,Op>* gh = (WeightedCFLOBDDInternalNode<T,Op> *) g.handleContents;
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

std::vector<double> ComputeProbabilityOfListNode(WeightedCFLOBDDNodeHandleT g, std::vector<std::vector<double>>& var_probs, std::vector<std::vector<double>>& path_probs, int start, int end){
	if (g.handleContents->level == 0){
		if (g == WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle){
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
	} else if (g == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g.handleContents->level]){
		std::vector<double> ans (path_probs[0].size(), 0);
		for (int i = 0; i <  path_probs[0].size(); i++)
			ans[i] = path_probs[0][i];
		return ans;
	} else {
		WeightedCFLOBDDInternalNode<T,Op>* gh = (WeightedCFLOBDDInternalNode<T,Op> *) g.handleContents;
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
*/
// TODO
// std::vector<double> ComputeEntropyOfListNode(WeightedCFLOBDDNodeHandleT g, std::vector<std::vector<double>>& var_probs, std::vector<std::vector<double>>& path_probs, std::vector<std::vector<double>>& entropy, int start, int end){
// 	if (g.handleContents->level == 0){
// 		if (g == WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle){
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
// 	} else if (g == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g.handleContents->level]){
// 		std::vector<double> ans (path_probs[0].size(), 0);
// 		for (int i = 0; i <  path_probs[0].size(); i++)
// 			ans[i] = path_probs[0][i];
// 		return ans;
// 	} else {
// 		WeightedCFLOBDDInternalNode<T,Op>* gh = (WeightedCFLOBDDInternalNode<T,Op> *) g.handleContents;
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
inline unsigned int binQuotient(unsigned int v, unsigned int level) {
  assert(level < 5);
  //  return v >> (1 << level);
  return v/(1 << (1 << level));
}

//
// binRemainder
//
// Return v mod (2**(2**level))
//
inline unsigned int binRemainder(unsigned int v, unsigned int level) {
  assert(level < 5);
  //  return v & (~((~0) << (1 << level)));
  return v % (1 << (1 << level));
}

namespace CFL_OBDD {
#ifdef ARBITRARY_STEP_FUNCTIONS
	WeightedCFLOBDDNodeHandleT MkStepNode(unsigned int level, unsigned int left, unsigned int middle, unsigned int right)
	{
		assert(level <= 5);    // Need LONG ints
		assert(middle == 0 || middle == 1);

		if (level == 0) {
			if (left == 2 || right == 2) {
				return WeightedCFLOBDDNodeHandleT::CFLOBDDDontCareNodeHandle;
			}
			else {
				return WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle;
			}
		}
		else {  // Create an appropriate WeightedCFLOBDDInternalNode<T,Op>
			WeightedCFLOBDDInternalNode<T,Op> *n = new WeightedCFLOBDDInternalNode<T,Op>(level);
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
			WeightedCFLOBDDNodeHandleT temp = MkStepNode(level-1, a, b, c);
			n->AConnection = WConnection<T,Op>(temp, m1);

			// Create the BConnections ------------------------------
			n->numBConnections = numberOfAExits;
			n->BConnection = new WConnection<T,Op>[n->numBConnections];
			int curConnection = 0;
			if (a != 0) {
				CFLOBDDReturnMapHandle m2;
				m2.AddToEnd(0);  // Connect to first exit
				m2.Canonicalize();
				n->BConnection[curConnection++] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m2);
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
				WeightedCFLOBDDNodeHandleT temp = MkStepNode(level-1, aa, bb, cc);
				n->BConnection[curConnection++] = WConnection<T,Op>(temp, m3);
			}
			if (c != 0) {
				CFLOBDDReturnMapHandle m4;
				m4.AddToEnd(numberOfExits-1);  // Connect to last exit
				m4.Canonicalize();
				n->BConnection[curConnection++] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[level-1], m4);
			}
			n->numExits = numberOfExits;
#ifdef PATH_COUNTING_ENABLED
			n->InstallPathCounts();
#endif
			return WeightedCFLOBDDNodeHandleT(n);
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
    /*
	WeightedCFLOBDDNodeHandleT shiftAtoB(WeightedCFLOBDDNodeHandleT f, const unsigned int levelAtWhichToShift)
	{
		assert(f.handleContents->level > 0);
		WeightedCFLOBDDInternalNode<T,Op> *g = (WeightedCFLOBDDInternalNode<T,Op> *)(f.handleContents);

		// Check that each BConnection of *g is a NoDistinctionNode[level-1]
		for (unsigned int j = 0; j < g->numBConnections; j++) {
			if (*(g->BConnection[j].entryPointHandle) != WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level - 1]) {
				std::cout << "g->level = " << g->level << std::endl;
				std::cout << *(g->BConnection[j].entryPointHandle) << std::endl << std::endl;
				std::cout << WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level - 1] << std::endl << std::endl;
			}
			assert(*(g->BConnection[j].entryPointHandle) == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level - 1]);
		}

		// Create an appropriate WeightedCFLOBDDInternalNode
		WeightedCFLOBDDInternalNode<T,Op> *n = new WeightedCFLOBDDInternalNode<T,Op>(g->level);
		if (g->level == levelAtWhichToShift) {  // Base case: the current level is consistent with shift
			// Put a NoDistinctionNode in the AConnection
			CFLOBDDReturnMapHandle m1;
			m1.AddToEnd(0);
			m1.Canonicalize();
			n->AConnection = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level - 1], m1);

			// Transfer g's AConnection to the BConnection
			n->numBConnections = 1;
			n->BConnection = new WConnection<T,Op>[n->numBConnections];
			n->BConnection[0] = g->AConnection;
		}
		else {    // Haven't reached the correct level yet, so . . .
			// Apply shiftAtoB recursively to g's AConnection
			WeightedCFLOBDDNodeHandleT temp = shiftAtoB(*(g->AConnection.entryPointHandle), levelAtWhichToShift);
			n->AConnection = WConnection<T,Op>(temp, g->AConnection.returnMapHandle);

			// Copy over all BConnections from g
			n->numBConnections = g->numBConnections;
			n->BConnection = new WConnection<T,Op>[n->numBConnections];
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				n->BConnection[j] = g->BConnection[j];
			}
		}
		n->numExits = g->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return WeightedCFLOBDDNodeHandleT(n);
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
	WeightedCFLOBDDNodeHandleT shiftBtoA(WeightedCFLOBDDNodeHandleT f, const unsigned int levelAtWhichToShift)
	{
		assert(f.handleContents->level > 0);
		WeightedCFLOBDDInternalNode<T,Op> *g = (WeightedCFLOBDDInternalNode<T,Op> *)(f.handleContents);

		if (f == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level]) {
			return f;
		}

		// Check that the AConnection of *g is a NoDistinctionNode[level-1]
		if (*(g->AConnection.entryPointHandle) != WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level - 1]) {
			std::cout << "g->level = " << g->level << std::endl;
			std::cout << *(g->AConnection.entryPointHandle) << std::endl << std::endl;
			std::cout << WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level - 1] << std::endl << std::endl;
		}
		assert(*(g->AConnection.entryPointHandle) == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level - 1]);

		// Create an appropriate WeightedCFLOBDDInternalNode<T,Op>
		WeightedCFLOBDDInternalNode<T,Op> *n = new WeightedCFLOBDDInternalNode<T,Op>(g->level);
		if (g->level == levelAtWhichToShift) {  // Base case: the current level is consistent with shift
			// Transfer g's BConnection[0] to the AConnection
			n->AConnection = g->BConnection[0];

			// Put a NoDistinctionNode in each BConnection
			n->numBConnections = g->numExits;
			n->BConnection = new WConnection<T,Op>[n->numBConnections];
			for (unsigned int j = 0; j < n->numBConnections; j++) {
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(j);
				m1.Canonicalize();
				n->BConnection[j] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level - 1], m1);
			}
		}
		else {    // Haven't reached the correct level yet, so . . .
			// Copy over the AConnection from g
			n->AConnection = g->AConnection;

			// Apply shiftBtoA recursively to g's BConnection[0]
			WeightedCFLOBDDNodeHandleT temp = shiftBtoA(*(g->BConnection[0].entryPointHandle), levelAtWhichToShift);
			n->numBConnections = g->numBConnections;
			n->BConnection = new WConnection<T,Op>[n->numBConnections];
			n->BConnection[0] = WConnection<T,Op>(temp, g->BConnection[0].returnMapHandle);
		}
		n->numExits = g->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return WeightedCFLOBDDNodeHandleT(n);
	}

	//
	// shiftAtoBAtLevelOne
	//
	// Given a WeightedCFLOBDDNodeHandleT for which all BConnections at level 1 are
	// NoDistinctionNode[0], create a corresponding WeightedCFLOBDDNodeHandleT
	// in which
	//  1. The AConnection is moved to BConnection[0], and
	//  2. All AConnections at level 1 are NoDistinctionNode[0]
	//
	// Precondition: *this is a WeightedCFLOBDDNodeHandleT for which all BConnections at level 1 are
	//               NoDistinctionNode[0]
	// Postcondition: *ans is a WeightedCFLOBDDNodeHandleT that satisfies conditions 1 and 2
	//
	// TODO: Checks on g == NoDistinctionNode[g->level] and function caching
	//
	// Hashtable<WeightedCFLOBDDNodeHandleT, WeightedCFLOBDDNodeHandleT> *memoTable
	//
	WeightedCFLOBDDNodeHandleT shiftAtoBAtLevelOne(Hashset<WeightedCFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, WeightedCFLOBDDNodeHandleT f)
	{
		assert(f.handleContents->level > 0);
		WeightedCFLOBDDInternalNode<T,Op> *g = (WeightedCFLOBDDInternalNode<T,Op> *)(f.handleContents);

		if (f == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level]) {
			return f;
		}
		totalVisitCount++;
		if (visitedNodes->Lookup(g) == NULL) {
			visitedNodes->Insert(g);
		}
		else {
			redundantVisitCount++;
		}

		// Create an appropriate WeightedCFLOBDDInternalNode<T,Op>
		WeightedCFLOBDDInternalNode<T,Op> *n = new WeightedCFLOBDDInternalNode<T,Op>(g->level);
		if (g->level == 1) {  // Base case: the current level is consistent with shift
			// Check that each BConnection of *g is a NoDistinctionNode[0]
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				if (*(g->BConnection[j].entryPointHandle) != WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0]) {
					std::cout << "g->level = " << g->level << std::endl;
					std::cout << *(g->BConnection[j].entryPointHandle) << std::endl << std::endl;
					std::cout << WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0] << std::endl << std::endl;
				}
				assert(*(g->BConnection[j].entryPointHandle) == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0]);
			}

			// Put a NoDistinctionNode[0] in the AConnection
			CFLOBDDReturnMapHandle m1;
			m1.AddToEnd(0);
			m1.Canonicalize();
			n->AConnection = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0], m1);

			// Transfer g's AConnection to BConnection[0]
			n->numBConnections = 1;
			n->BConnection = new WConnection<T,Op>[1];
			n->BConnection[0] = g->AConnection;
		}
		else {    // Haven't reached the correct level yet, so apply shiftAtoBAtLevelOne recursively
			// Create the AConnection
			WeightedCFLOBDDNodeHandleT temp1 = shiftAtoBAtLevelOne(visitedNodes, totalVisitCount, redundantVisitCount, *(g->AConnection.entryPointHandle));
			n->AConnection = WConnection<T,Op>(temp1, g->AConnection.returnMapHandle);

			// Create the BConnections
			n->numBConnections = g->numBConnections;
			n->BConnection = new WConnection<T,Op>[n->numBConnections];
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				WeightedCFLOBDDNodeHandleT temp2 = shiftAtoBAtLevelOne(visitedNodes, totalVisitCount, redundantVisitCount, *(g->BConnection[j].entryPointHandle));
				n->BConnection[j] = WConnection<T,Op>(temp2, g->BConnection[j].returnMapHandle);
			}
		}
		n->numExits = g->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return WeightedCFLOBDDNodeHandleT(n);
	}

	//
	// shiftBtoAAtLevelOne
	//
	// Given a WeightedCFLOBDDNodeHandleT for which all AConnections at level 1 are
	// NoDistinctionNode[0], create a corresponding WeightedCFLOBDDNodeHandleT
	// in which
	//  1. The BConnection[0] is moved to AConnection, and
	//  2. All BConnections at level 1 are NoDistinctionNode[0]
	//
	// Precondition: *this is a WeightedCFLOBDDNodeHandleT for which all AConnections at level 1 are
	//               NoDistinctionNode[0]
	// Postcondition: *ans is a WeightedCFLOBDDNodeHandleT that satisfies conditions 1 and 2
	//
	// TODO: Checks on g == NoDistinctionNode[g->level] and function caching
	//
	WeightedCFLOBDDNodeHandleT shiftBtoAAtLevelOne(Hashset<WeightedCFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, WeightedCFLOBDDNodeHandleT f)
	{
		assert(f.handleContents->level > 0);
		WeightedCFLOBDDInternalNode<T,Op> *g = (WeightedCFLOBDDInternalNode<T,Op> *)(f.handleContents);

		if (f == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[g->level]) {
			return f;
		}
		totalVisitCount++;
		if (visitedNodes->Lookup(g) == NULL) {
			visitedNodes->Insert(g);
		}
		else {
			redundantVisitCount++;
		}

		// Create an appropriate WeightedCFLOBDDInternalNode<T,Op>
		WeightedCFLOBDDInternalNode<T,Op> *n = new WeightedCFLOBDDInternalNode<T,Op>(g->level);
		if (g->level == 1) {  // Base case: the current level is consistent with shift
			// Check that the AConnection of *g is a NoDistinctionNode[0]
			if (*(g->AConnection.entryPointHandle) != WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0]) {
				std::cout << "g->level = " << g->level << std::endl;
				std::cout << *(g->AConnection.entryPointHandle) << std::endl << std::endl;
				std::cout << WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0] << std::endl << std::endl;
			}
			assert(*(g->AConnection.entryPointHandle) == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0]);

			// Transfer g's BConnection[0] to the AConnection
			n->AConnection = g->BConnection[0];

			// Put a NoDistinctionNode in each BConnection
			n->numBConnections = g->numExits;
			n->BConnection = new WConnection<T,Op>[n->numBConnections];
			for (unsigned int j = 0; j < n->numBConnections; j++) {
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(j);
				m1.Canonicalize();
				n->BConnection[j] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0], m1);
			}
		}
		else {    // Haven't reached the correct level yet, so apply shiftBtoAAtLevelOne recursively
			// Create the AConnection
			WeightedCFLOBDDNodeHandleT temp1 = shiftBtoAAtLevelOne(visitedNodes, totalVisitCount, redundantVisitCount, *(g->AConnection.entryPointHandle));
			n->AConnection = WConnection<T,Op>(temp1, g->AConnection.returnMapHandle);

			// Create the BConnections
			n->numBConnections = g->numBConnections;
			n->BConnection = new WConnection<T,Op>[n->numBConnections];
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				WeightedCFLOBDDNodeHandleT temp2 = shiftBtoAAtLevelOne(visitedNodes, totalVisitCount, redundantVisitCount, *(g->BConnection[j].entryPointHandle));
				n->BConnection[j] = WConnection<T,Op>(temp2, g->BConnection[j].returnMapHandle);
			}
		}
		n->numExits = g->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return WeightedCFLOBDDNodeHandleT(n);
	}

	//
	// duplicateAinBAtLevelOne
	//
	// Given a WeightedCFLOBDDNodeHandleT for which all BConnections at level 1 are
	// NoDistinctionNode[0], create a corresponding WeightedCFLOBDDNodeHandleT
	// in which
	//  1. The AConnection is duplicated in BConnection[0], and
	//  2. The AConnection and BConnection[0] "match"
	// "Match" means that if the AConnection's handle is WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle
	// then the three assignments of the level-one Boolean variables (0,0), (0,1), and (1,0)
	// are routed to exit[0] and (1,1) is routed to exit[1].
	//
	// Precondition: *this is a WeightedCFLOBDDNodeHandleT for which all BConnections at level 1 are
	//               NoDistinctionNode[0]
	// Postcondition: *ans is a WeightedCFLOBDDNodeHandleT that satisfies conditions 1 and 2
	//
	// TODO: Checks on g == NoDistinctionNode[g->level] and function caching
	//
	WeightedCFLOBDDNodeHandleT duplicateAinBAtLevelOne(WeightedCFLOBDDNodeHandleT f)
	{
		assert(f.handleContents->level > 0);
		WeightedCFLOBDDInternalNode<T,Op> *g = (WeightedCFLOBDDInternalNode<T,Op> *)(f.handleContents);

		// Create an appropriate WeightedCFLOBDDInternalNode<T,Op>
		WeightedCFLOBDDInternalNode<T,Op> *n = new WeightedCFLOBDDInternalNode<T,Op>(g->level);
		if (g->level == 1) {  // Base case: the level at which to create a node with duplicated AConnection
			// Check that each BConnection of *g is a NoDistinctionNode[0]
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				if (*(g->BConnection[j].entryPointHandle) != WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0]) {
					std::cout << "g->level = " << g->level << std::endl;
					std::cout << g->BConnection[j].entryPointHandle << std::endl << std::endl;
					std::cout << WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0] << std::endl << std::endl;
				}
				assert(*(g->BConnection[j].entryPointHandle) == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0]);
			}

			// Transfer g's AConnection to n
			n->AConnection = g->AConnection;

			// Create an appropriate set of BConnections
			n->numBConnections = g->numBConnections;
			n->BConnection = new WConnection<T,Op>[n->numBConnections];
			if (g->numBConnections == 1) {
				assert(*(g->AConnection.entryPointHandle) == WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0]);
				// Transfer g's BConnection[0] to n; n should be equal to WeightedCFLOBDDNodeHandleT::NoDistinctionNode[1]
				n->BConnection[0] = g->AConnection;
			}
			else {
				assert(g->numBConnections == 2 && *(g->AConnection.entryPointHandle) == WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle);
				// Put a NoDistinctionNode[0] in BConnection[0]
				CFLOBDDReturnMapHandle m0;
				m0.AddToEnd(0);
				m0.Canonicalize();
				n->BConnection[0] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::NoDistinctionNode[0], m0);

				// Put a ForkNode in BConnection[1]
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m01.Canonicalize();
				n->BConnection[1] = WConnection<T,Op>(WeightedCFLOBDDNodeHandleT::CFLOBDDForkNodeHandle, m01);
				assert(n->BConnection[1] == g->AConnection);
			}
		}
		else {    // Haven't reached the correct level yet, so apply duplicateAinBAtLevelOne recursively
			// Create the AConnection
			WeightedCFLOBDDNodeHandleT temp1 = duplicateAinBAtLevelOne(*(g->AConnection.entryPointHandle));
			n->AConnection = WConnection<T,Op>(temp1, g->AConnection.returnMapHandle);

			// Create the BConnections
			n->numBConnections = g->numBConnections;
			n->BConnection = new WConnection<T,Op>[n->numBConnections];
			for (unsigned int j = 0; j < g->numBConnections; j++) {
				WeightedCFLOBDDNodeHandleT temp2 = duplicateAinBAtLevelOne(*(g->BConnection[j].entryPointHandle));
				n->BConnection[j] = WConnection<T,Op>(temp2, g->BConnection[j].returnMapHandle);
			}
		}
		n->numExits = g->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return WeightedCFLOBDDNodeHandleT(n);
	}
    */

}  // namespace CFL_OBDD

//********************************************************************
// CFLReduceKey
//********************************************************************

// Constructor
template <typename T, typename Op>
WeightedCFLReduceKey<T,Op>::WeightedCFLReduceKey(WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle, ReductionMapHandle redMapHandle, WeightedValuesListHandle<T> valList)
  :  nodeHandle(nodeHandle), redMapHandle(redMapHandle), valList(valList)
{
}

// Hash
template <typename T, typename Op>
unsigned int WeightedCFLReduceKey<T,Op>::Hash(unsigned int modsize)
{
  unsigned int hvalue = 0;
  hvalue = (997 * nodeHandle.Hash(modsize) + 97 * redMapHandle.Hash(modsize) + valList.Hash(modsize)) % modsize;
  return hvalue;
}

// print
template <typename T, typename Op>
std::ostream& WeightedCFLReduceKey<T,Op>::print(std::ostream & out) const
{
  out << "(" << nodeHandle << ", " << redMapHandle << ", " << valList << ")";
  return out;
}

template <typename T, typename Op>
std::ostream& operator<< (std::ostream & out, const WeightedCFLReduceKey<T,Op> &p)
{
  p.print(out);
  return(out);
}

template <typename T, typename Op>
WeightedCFLReduceKey<T,Op>& WeightedCFLReduceKey<T,Op>::operator= (const WeightedCFLReduceKey<T,Op>& i)
{
  if (this != &i)      // don't assign to self!
  {
    nodeHandle = i.nodeHandle;
    redMapHandle = i.redMapHandle;
    valList = i.valList;
  }
  return *this;        
}

// Overloaded !=
template <typename T, typename Op>
bool WeightedCFLReduceKey<T,Op>::operator!=(const WeightedCFLReduceKey<T,Op>& p)
{
  return (nodeHandle != p.nodeHandle) || (redMapHandle != p.redMapHandle) || (valList != p.valList);
}

// Overloaded ==
template <typename T, typename Op>
bool WeightedCFLReduceKey<T,Op>::operator==(const WeightedCFLReduceKey<T,Op>& p)
{
  return (nodeHandle == p.nodeHandle) && (redMapHandle == p.redMapHandle) && (valList == p.valList);
}

//********************************************************************
// WeightedCFLOBDDNode
//********************************************************************

// Initializations of static members ---------------------------------

template <typename T, typename Op>
unsigned int const WeightedCFLOBDDNode<T,Op>::maxLevel = CFLOBDDMaxlevel;

// Constructors/Destructor -------------------------------------------

// Default constructor
template <typename T, typename Op>
WeightedCFLOBDDNode<T,Op>::WeightedCFLOBDDNode()
	: level(maxLevel), refCount(0), isCanonical(false), isNumPathsMemAllocated(false), isWeightsOfPathsMemAllocated(false)
{
}

// Constructor
template <typename T, typename Op>
WeightedCFLOBDDNode<T,Op>::WeightedCFLOBDDNode(const unsigned int l)
	: level(l), refCount(0), isCanonical(false), isNumPathsMemAllocated(false), isWeightsOfPathsMemAllocated(false)
{
}

template <typename T, typename Op>
WeightedCFLOBDDNode<T,Op>::~WeightedCFLOBDDNode()
{
}

// print
namespace CFL_OBDD {
    template <typename T, typename Op>
std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDNode<T,Op> &n)
{
  n.print(out);
  return(out);
}
}

//********************************************************************
// WeightedCFLOBDDInternalNode
//********************************************************************

// Constructors/Destructor -------------------------------------------

template <typename T, typename Op>
WeightedCFLOBDDInternalNode<T,Op>::WeightedCFLOBDDInternalNode(const unsigned int l)
  :  WeightedCFLOBDDNode<T,Op>(l)
{
}

template <typename T, typename Op>
WeightedCFLOBDDInternalNode<T,Op>::~WeightedCFLOBDDInternalNode()
{
  delete [] BConnection;
//#ifdef PATH_COUNTING_ENABLED
  if (this->isNumPathsMemAllocated)
	delete [] this->numPathsToExit;
  if (this->isWeightsOfPathsMemAllocated)
    delete [] this->numWeightsOfPathsAsAmpsToExit;
//#endif
}

// print
// TODO: Check this
template <typename T, typename Op>
std::ostream& WeightedCFLOBDDInternalNode<T,Op>::print(std::ostream & out) const
{
    unsigned int i, j;
    unsigned int level = this->level;
    unsigned int maxLevel = this->maxLevel;
    for (i = level; i < maxLevel; i++) {  // Indentation
        out << "  ";
    }
    out << "A: " << std::endl;
    if (WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[level - 1] == *AConnection.entryPointHandle) {
	    for (i = level-1; i < maxLevel; i++) {  // Indentation
		    out << "  ";
	    }
	    out << "NoDistinctionNode[" << level - 1 << "] " << getIdentityValue<T,Op>() << std::endl;
    }
    else if (WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[level - 1] == *AConnection.entryPointHandle) {
        for (i = level-1; i < maxLevel; i++) {  // Indentation
		    out << "  ";
	    }
	    out << "NoDistinctionNode[" << level - 1 << "] " << getAnnhilatorValue<T,Op>() << std::endl;
    }
    else if (WeightedCFLOBDDNodeHandleT<T,Op>::IdentityNode[level - 1] == *AConnection.entryPointHandle) {
        for (i = level-1; i < maxLevel; i++) {  // Indentation
		    out << "  ";
	    }
	    out << "IdentityNode[" << level - 1 << "] " << std::endl;
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
	    if (WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[level - 1] == *BConnection[j].entryPointHandle) {
		    for (i = level - 1; i < maxLevel; i++) {  // Indentation
			    out << "  ";
		    }
			out << "NoDistinctionNode[" << level - 1 << "] " << getIdentityValue<T,Op>() << std::endl;
        }
        else if (WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[level - 1] == *BConnection[j].entryPointHandle) {
		    for (i = level - 1; i < maxLevel; i++) {  // Indentation
			    out << "  ";
		    }
			out << "NoDistinctionNode[" << level - 1 << "] " << getAnnhilatorValue<T,Op>() << std::endl;
        }
        else if (WeightedCFLOBDDNodeHandleT<T,Op>::IdentityNode[level - 1] == *BConnection[j].entryPointHandle) {
           for (i = level-1; i < maxLevel; i++) {  // Indentation
                out << "  ";
            }
            out << "IdentityNode[" << level - 1 << "] " << std::endl; 
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

// template <typename T, typename Op>
// bool WeightedCFLOBDDInternalNode<T,Op>::IsValid()
// {
// 	if (isValid)
// 		return true;
	
// 	if (!(level == (AConnection.entryPointHandle->handleContents->level + 1)))
// 		return false;
// 	for (unsigned int i = 0; i < numBConnections; i++){
// 		if (!(level == (BConnection[i].entryPointHandle->handleContents->level + 1)))
// 			return false;
// 	}

// 	for (unsigned int i = 0; i < AConnection.returnMapHandle.Size(); i++){
// 		if (AConnection.returnMapHandle[i] != i)
// 			return false;
// 	}

// 	unsigned int tempNumExits = 0;
// 	for (unsigned int i = 0; i < numBConnections; i++){
// 		std::unordered_set<int> s;
// 		for (unsigned int j = 0; j < BConnection[i].returnMapHandle.Size(); j++){
// 			if (s.find(BConnection[i].returnMapHandle[j]) != s.end())
// 				return false;
// 			s.insert(BConnection[i].returnMapHandle[j]);
// 			if (BConnection[i].returnMapHandle[j] > tempNumExits)
// 				return false;
// 			else if (BConnection[i].returnMapHandle[j] == tempNumExits)
// 				tempNumExits++;
// 		}
// 	}

// 	if (numExits != tempNumExits)
// 		return false;

// 	if (!AConnection.entryPointHandle->handleContents->IsValid())
// 		return false;
// 	for (unsigned int i = 0; i < numBConnections; i++){
// 		if (!BConnection[i].entryPointHandle->handleContents->IsValid())
// 			return false;
// 	}
	
// 	isValid = true;
// 	return true;
// }

template <typename T, typename Op>
void WeightedCFLOBDDInternalNode<T,Op>::FillSatisfyingAssignment(unsigned int exitNumber, SH_OBDD::Assignment &assignment, unsigned int &index)
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

template <typename T, typename Op>
int WeightedCFLOBDDInternalNode<T,Op>::Traverse(SH_OBDD::AssignmentIterator &ai)
{
  int ans;

  if (this == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[this->level].handleContents) {
    ai.Advance((unsigned int)(((unsigned int)1) << this->level));
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


inline CFLOBDDReturnMapHandle ComposeAndReduce(CFLOBDDReturnMapHandle& mapHandle, ReductionMapHandle& redMapHandle, ReductionMapHandle& inducedRedMapHandle)
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

template <typename T>
WeightedValuesListHandle<T> ComposeValueList(CFLOBDDReturnMapHandle& returnMap, WeightedValuesListHandle<T>& valList)
{
    WeightedValuesListHandle<T> answerList;
    for (unsigned int i = 0; i < returnMap.Size(); i++)
    {
        answerList.AddToEnd(valList[returnMap[i]]);
    }
    return answerList;
}

// Check this
template <typename T, typename Op>
std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> WeightedCFLOBDDInternalNode<T,Op>::Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valList, bool forceReduce)
{
  WeightedCFLOBDDInternalNode *n = new WeightedCFLOBDDInternalNode(this->level);
  ReductionMapHandle AReductionMapHandle;        // To record duplicate BConnections
  WeightedValuesListHandle<T> AValueList;

  // Reduce the B connections
     n->BConnection = new WConnection<T,Op>[numBConnections];   // May create shorter version later
     n->numBConnections = 0;
     for (unsigned int i = 0; i < numBConnections; i++) {
        ReductionMapHandle inducedReductionMapHandle(redMapHandle.Size());
        CFLOBDDReturnMapHandle inducedReturnMap;
		inducedReturnMap = ComposeAndReduce(BConnection[i].returnMapHandle, redMapHandle, inducedReductionMapHandle);
        WeightedValuesListHandle<T> inducedValueList = ComposeValueList(BConnection[i].returnMapHandle, valList);
        std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> temp = BConnection[i].entryPointHandle->Reduce(inducedReductionMapHandle, inducedReturnMap.Size(), inducedValueList, forceReduce);
        WConnection<T,Op> c(temp.first, inducedReturnMap);
        unsigned int position = n->InsertBConnection(n->numBConnections, c);
        AReductionMapHandle.AddToEnd(position);
        AValueList.AddToEnd(temp.second);
     }
     AReductionMapHandle.Canonicalize();
     AValueList.Canonicalize();
     if (n->numBConnections < numBConnections) {  // Shorten
       WConnection<T,Op> *temp = n->BConnection;
       n->BConnection = new WConnection<T,Op>[n->numBConnections];
       for (unsigned int j = 0; j < n->numBConnections; j++) {
         n->BConnection[j] = temp[j];
       }
       delete [] temp;
     }

  // Reduce the A connection
     ReductionMapHandle inducedAReductionMapHandle;
	 CFLOBDDReturnMapHandle inducedAReturnMap;
	 inducedAReturnMap = ComposeAndReduce(AConnection.returnMapHandle, AReductionMapHandle, inducedAReductionMapHandle);
     WeightedValuesListHandle<T> inducedAValueList = ComposeValueList(AConnection.returnMapHandle, AValueList);
     std::pair<WeightedCFLOBDDNodeHandleT<T,Op>, T> tempHandle = AConnection.entryPointHandle->Reduce(inducedAReductionMapHandle, inducedAReturnMap.Size(), inducedAValueList, forceReduce);
     n->AConnection = WConnection<T,Op>(tempHandle.first, inducedAReturnMap);

  // Other material that has to be filled in
     n->numExits = replacementNumExits;
#ifdef PATH_COUNTING_ENABLED
     n->InstallPathCounts();
#endif
    return std::make_pair(WeightedCFLOBDDNodeHandleT(n), tempHandle.second);
} // CFLOBDDInternalNode::Reduce

template <typename T, typename Op>
unsigned int WeightedCFLOBDDInternalNode<T,Op>::Hash(unsigned int modsize)
{
  unsigned int hvalue = AConnection.Hash(modsize);
  for (unsigned int j = 0; j < numBConnections; j++) {
    hvalue = (997 * hvalue + BConnection[j].Hash(modsize)) % modsize;
  }
  return hvalue;
}

template <typename T, typename Op>
void WeightedCFLOBDDInternalNode<T,Op>::DumpConnections(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visited, std::ostream & out /* = std::cout */)
{
	if (visited->Lookup(new WeightedCFLOBDDNodeHandleT(this)) == NULL) {
    unsigned int i;
	visited->Insert(new WeightedCFLOBDDNodeHandleT(this));
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

template <typename T, typename Op>
void WeightedCFLOBDDInternalNode<T,Op>::CountNodesAndEdges(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
	unsigned int &nodeCount, unsigned int &edgeCount, unsigned int &returnEdgesCount)
{
  if (visitedNodes->Lookup(new WeightedCFLOBDDNodeHandleT(this)) == NULL) {
    visitedNodes->Insert(new WeightedCFLOBDDNodeHandleT(this));
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

template <typename T, typename Op>
void WeightedCFLOBDDInternalNode<T,Op>::CountNodes(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, unsigned int &nodeCount)
{
	if (visitedNodes->Lookup(new WeightedCFLOBDDNodeHandleT(this)) == NULL) {
		visitedNodes->Insert(new WeightedCFLOBDDNodeHandleT(this));
		nodeCount++;
		AConnection.entryPointHandle->handleContents->CountNodes(visitedNodes, nodeCount);
		for (unsigned int i = 0; i < numBConnections; i++) {
			BConnection[i].entryPointHandle->handleContents->CountNodes(visitedNodes, nodeCount);
		}
	}
}

template <typename T, typename Op>
void WeightedCFLOBDDInternalNode<T,Op>::CountPaths(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes)
{
	WeightedCFLOBDDNodeHandleT<T,Op>* handle = new WeightedCFLOBDDNodeHandleT<T,Op>(this);
	if (visitedNodes->Lookup(handle) == NULL) {
		visitedNodes->Insert(handle);
		AConnection.entryPointHandle->handleContents->CountPaths(visitedNodes);
		for (unsigned int i = 0; i < numBConnections; i++) {
			BConnection[i].entryPointHandle->handleContents->CountPaths(visitedNodes);
		}
		InstallPathCounts();
	}
}

template <typename T, typename Op>
void WeightedCFLOBDDInternalNode<T,Op>::ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes)
{
	WeightedCFLOBDDNodeHandleT<T,Op>* handle = new WeightedCFLOBDDNodeHandleT<T,Op>(this);
	if (visitedNodes->Lookup(handle) == NULL) {
		visitedNodes->Insert(handle);
		AConnection.entryPointHandle->handleContents->ComputeWeightOfPathsAsAmpsToExits(visitedNodes);
		for (unsigned int i = 0; i < numBConnections; i++) {
			BConnection[i].entryPointHandle->handleContents->ComputeWeightOfPathsAsAmpsToExits(visitedNodes);
		}
		InstallWeightsOfPathsAsAmpsToExits();
	}
}


// Overloaded !=
template <typename T, typename Op>
bool WeightedCFLOBDDInternalNode<T,Op>::operator!= (const WeightedCFLOBDDNode<T,Op> & n)
{
  return !(*this == n);
}

// Overloaded ==
template <typename T, typename Op>
bool WeightedCFLOBDDInternalNode<T,Op>::operator== (const WeightedCFLOBDDNode<T,Op> & n)
{
  if (n.NodeKind() != W_CFLOBDD_INTERNAL)
    return false;
  WeightedCFLOBDDInternalNode &m = (WeightedCFLOBDDInternalNode &)n;
  if (this->level != m.level)
    return false;
  if (this->numExits != m.numExits)
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

template <typename T, typename Op>
void WeightedCFLOBDDInternalNode<T,Op>::IncrRef()
{
  this->refCount++;    // Warning: Saturation not checked
}

template <typename T, typename Op>
void WeightedCFLOBDDInternalNode<T,Op>::DecrRef()
{
  if (--this->refCount == 0) {    // Warning: Saturation not checked
    if (this->isCanonical) {
      WeightedCFLOBDDNodeHandleT<T,Op>::canonicalNodeTable->DeleteEq(this);
    }
    delete this;
  }
}

// Insert Connection c at BConnection position j, but only if c does not
// duplicate an existing BConnection.
// Return the position at which c was found or inserted.
template <typename T, typename Op>
unsigned int WeightedCFLOBDDInternalNode<T,Op>::InsertBConnection(unsigned int &j, WConnection<T,Op> &c)
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

inline long double addNumPathsToExit(long double path1, long double path2){
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

inline long double addNumPathsToExit(std::vector<long double>& logOfPaths){
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
template <typename T, typename Op>
void WeightedCFLOBDDInternalNode<T,Op>::InstallPathCounts()
{
  this->numPathsToExit = new long double[this->numExits];
  this->isNumPathsMemAllocated = true;
  for (unsigned int i = 0; i < this->numExits; i++) {
    this->numPathsToExit[i] = -1.0 * std::numeric_limits<long double>::infinity();
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
		this->numPathsToExit[it->first] = addNumPathsToExit(it->second);
	}
  }
}


// InstallWeightsOfPathsAsAmpsToExits
template <typename T, typename Op>
void WeightedCFLOBDDInternalNode<T,Op>::InstallWeightsOfPathsAsAmpsToExits()
{
  this->numWeightsOfPathsAsAmpsToExit = new long double[this->numExits];
  this->isWeightsOfPathsMemAllocated = true;
  for (unsigned int i = 0; i < this->numExits; i++) {
    this->numWeightsOfPathsAsAmpsToExit[i] = 0;
  }

//   std::map<unsigned int, std::vector<T>> storingNumPathsToExit;

  for (unsigned int i = 0; i < AConnection.entryPointHandle->handleContents->numExits; i++) {
    for (unsigned int j = 0; j < BConnection[i].entryPointHandle->handleContents->numExits; j++) {
      unsigned int k = BConnection[i].returnMapHandle.Lookup(j);
	  long double numPathsValue = AConnection.entryPointHandle->handleContents->numWeightsOfPathsAsAmpsToExit[i] * BConnection[i].entryPointHandle->handleContents->numWeightsOfPathsAsAmpsToExit[j];
	//   if (storingNumPathsToExit.find(k) == storingNumPathsToExit.end()){
	// 	  std::vector<long double> logOfPaths;
	// 	  logOfPaths.push_back(numPathsValue);
	// 	  storingNumPathsToExit[k] = logOfPaths;
	//   }
	//   else{
	// 	  storingNumPathsToExit[k].push_back(numPathsValue);
	//   }
	  this->numWeightsOfPathsAsAmpsToExit[k] = this->numWeightsOfPathsAsAmpsToExit[k] + numPathsValue;
    }

	// for (std::map<unsigned int, std::vector<long double>>::iterator it = storingNumPathsToExit.begin(); it != storingNumPathsToExit.end(); it++){
	// 	sort(it->second.begin(), it->second.end());
	// 	this->numPathsToExit[it->first] = addNumPathsToExit(it->second);
	// }
  }
}



//#endif

//********************************************************************
// CFLOBDDLeafNode
//********************************************************************

// Constructors/Destructor -------------------------------------------

// Default constructor
template <typename T, typename Op>
WeightedCFLOBDDLeafNode<T, Op>::WeightedCFLOBDDLeafNode(bool flag)
  :  WeightedCFLOBDDNode<T,Op>(0)
{
  this->refCount = 1;
  if (flag)
    lweight = rweight = getIdentityValue<T, Op>();
  else
    lweight = rweight = getAnnhilatorValue<T, Op>();
}

template <typename T, typename Op>
WeightedCFLOBDDLeafNode<T, Op>::WeightedCFLOBDDLeafNode(T lw, T rw)
  :  WeightedCFLOBDDNode<T,Op>(0)
{
  this->refCount = 1;
  lweight = lw;
  rweight = rw;
}

template <typename T, typename Op>
WeightedCFLOBDDLeafNode<T, Op>::~WeightedCFLOBDDLeafNode()
{
}

template <typename T, typename Op>
void WeightedCFLOBDDLeafNode<T, Op>::DumpConnections(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *, std::ostream & /* = std::cout */ )
{
}

template <typename T, typename Op>
void WeightedCFLOBDDLeafNode<T, Op>::CountNodesAndEdges(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *, unsigned int &nodeCount, 
	unsigned int &, unsigned int& )
{
  if (visitedNodes->Lookup(new WeightedCFLOBDDNodeHandleT<T,Op>(this)) == NULL) {
    visitedNodes->Insert(new WeightedCFLOBDDNodeHandleT<T,Op>(this));
    nodeCount++;
  }
}

template <typename T, typename Op>
void WeightedCFLOBDDLeafNode<T,Op>::CountNodes(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, unsigned int &nodeCount)
{
	if (visitedNodes->Lookup(new WeightedCFLOBDDNodeHandleT<T,Op>(this)) == NULL) {
		visitedNodes->Insert(new WeightedCFLOBDDNodeHandleT<T,Op>(this));
		nodeCount++;
	}
}

template <typename T, typename Op>
void WeightedCFLOBDDLeafNode<T,Op>::CountPaths(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes)
{
	WeightedCFLOBDDNodeHandleT<T,Op>* handle = new WeightedCFLOBDDNodeHandleT<T,Op>(this);
	if (visitedNodes->Lookup(handle) == NULL) {
		visitedNodes->Insert(handle);
	}
}

template <typename T, typename Op>
void WeightedCFLOBDDLeafNode<T,Op>::ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes)
{
	WeightedCFLOBDDNodeHandleT<T,Op>* handle = new WeightedCFLOBDDNodeHandleT<T,Op>(this);
	if (visitedNodes->Lookup(handle) == NULL) {
		visitedNodes->Insert(handle);
	}
    InstallWeightsOfPathsAsAmpsToExits();
}

template <typename T, typename Op>
void WeightedCFLOBDDLeafNode<T,Op>::IncrRef() 
{
    this->refCount++;
}
template <typename T, typename Op>
void WeightedCFLOBDDLeafNode<T,Op>::DecrRef() 
{ 
  if (--this->refCount == 0) {    // Warning: Saturation not checked
    if (this->isCanonical) {
      WeightedCFLOBDDNodeHandleT<T,Op>::canonicalNodeTable->DeleteEq(this);
    }
    delete this;
  }  
}

//********************************************************************
// CFLOBDDForkNode
//********************************************************************

// Constructors/Destructor -------------------------------------------

// Default constructor
template <typename T, typename Op>
WeightedCFLOBDDForkNode<T,Op>::WeightedCFLOBDDForkNode(bool flag)
  :  WeightedCFLOBDDLeafNode<T,Op>(flag)
{
  this->numExits = 2;
  this->numPathsToExit = new long double[2];
  this->numPathsToExit[0] = 0;
  this->numPathsToExit[1] = 0;
}

template <typename T, typename Op>
WeightedCFLOBDDForkNode<T,Op>::WeightedCFLOBDDForkNode(T lweight, T rweight)
  :  WeightedCFLOBDDLeafNode<T,Op>(lweight, rweight)
{
  this->numExits = 2;
  this->numPathsToExit = new long double[2];
  this->numPathsToExit[0] = 0;
  this->numPathsToExit[1] = 0;
}

template <typename T, typename Op>
WeightedCFLOBDDForkNode<T,Op>::~WeightedCFLOBDDForkNode()
{
    if (this->isWeightsOfPathsMemAllocated)
        delete [] this->numWeightsOfPathsAsAmpsToExit;
}

template <typename T, typename Op>
void WeightedCFLOBDDForkNode<T,Op>::InstallWeightsOfPathsAsAmpsToExits()
{
    this->isWeightsOfPathsMemAllocated = true;
    this->numWeightsOfPathsAsAmpsToExit = new long double[2];
    this->numWeightsOfPathsAsAmpsToExit[0] = computeProbabilityFromAmplitude<T,Op>(this->lweight);
    this->numWeightsOfPathsAsAmpsToExit[1] = computeProbabilityFromAmplitude<T,Op>(this->rweight);
}

// print
template <typename T, typename Op>
std::ostream& WeightedCFLOBDDForkNode<T, Op>::print(std::ostream & out) const
{
  for (unsigned int i = this->level; i < this->maxLevel; i++) {
    out << "  ";
  }
  out << "Fork (" << this->lweight << "," << this->rweight << ")";
  return out;
}

// template <typename T, typename Op>
// bool WeightedCFLOBDDForkNode<T,Op>::IsValid()
// {
// 	if (level == 0)
// 		isValid = true;
// 	return level == 0;
// }

template <typename T, typename Op>
void WeightedCFLOBDDForkNode<T,Op>::FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index)
{
  assert(i <= 1);
  index--;
  assignment[index] = i;
}

template <typename T, typename Op>
int WeightedCFLOBDDForkNode<T,Op>::Traverse(SH_OBDD::AssignmentIterator &ai)
{
  bool val = ai.Current();
  ai.Next();
  return (int)val;
}

template <typename T, typename Op>
std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> WeightedCFLOBDDForkNode<T,Op>::Reduce(ReductionMapHandle&, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valList, bool forceReduce)
{
    // std::cout << "Fork " << replacementNumExits << " " << valList << std::endl;
    if (replacementNumExits == 1) {
        if (valList[0] == valList[1] && valList[0] == getIdentityValue<T,Op>())
            return std::make_pair(WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDDontCareNodeHandle, getIdentityValue<T,Op>());
        else if (valList[0] == valList[1] && valList[0] == getAnnhilatorValue<T,Op>())
            return std::make_pair(WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[0], valList[0]);
        else{
            T lw = valList[0] * this->lweight;
            T rw = valList[1] * this->rweight;
            std::tuple<T,T,T> w = computeInverseValue<T, Op>(lw, rw);
            if (std::get<1>(w) == getIdentityValue<T,Op>() && std::get<2>(w) == getIdentityValue<T,Op>())
                return std::make_pair(WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDDontCareNodeHandle, std::get<0>(w));
            WeightedCFLOBDDDontCareNode<T,Op>* n = new WeightedCFLOBDDDontCareNode<T,Op>(std::get<1>(w), std::get<2>(w));
            return std::make_pair(WeightedCFLOBDDNodeHandleT(n), std::get<0>(w));
        }
    }
    else {
        if (valList[0] == valList[1] && valList[0] == getIdentityValue<T,Op>())
            return std::make_pair(WeightedCFLOBDDNodeHandleT(this), getIdentityValue<T,Op>());
        else{
            T lw = valList[0] * this->lweight;
            T rw = valList[1] * this->rweight;
            std::tuple<T,T,T> w = computeInverseValue<T, Op>(lw, rw);
            if (std::get<1>(w) == getIdentityValue<T,Op>() && std::get<2>(w) == getIdentityValue<T,Op>())
                return std::make_pair(WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle, std::get<0>(w));
            if (std::get<1>(w) == getIdentityValue<T,Op>() && std::get<2>(w) == getAnnhilatorValue<T,Op>())
                return std::make_pair(WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle10, std::get<0>(w));
            if (std::get<1>(w) == getAnnhilatorValue<T,Op>() && std::get<2>(w) == getIdentityValue<T,Op>())
                return std::make_pair(WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle01, std::get<0>(w));
            WeightedCFLOBDDForkNode<T,Op>* n = new WeightedCFLOBDDForkNode<T,Op>(std::get<1>(w), std::get<2>(w));
            return std::make_pair(WeightedCFLOBDDNodeHandleT(n), std::get<0>(w));
        }
    }
}

template <typename T, typename Op>
unsigned int WeightedCFLOBDDForkNode<T,Op>::Hash(unsigned int modsize)
{
    boost::hash<T> boost_hash;
    return (997 * boost_hash(this->lweight) + 97 * boost_hash(this->rweight) + 2) % modsize; 
//   return (((unsigned int)reinterpret_cast<uintptr_t>(this) >> 2) + 997 * boost_hash(this->lweight) + 97 * boost_hash(this->rweight)) % modsize;
}

// Overloaded !=
template <typename T, typename Op>
bool WeightedCFLOBDDForkNode<T,Op>::operator!= (const WeightedCFLOBDDNode<T,Op> & n)
{
  return n.NodeKind() != W_CFLOBDD_FORK;
}

// Overloaded ==
template <typename T, typename Op>
bool WeightedCFLOBDDForkNode<T,Op>::operator== (const WeightedCFLOBDDNode<T,Op> & n)
{
  return n.NodeKind() == W_CFLOBDD_FORK;
}

//********************************************************************
// WeightedCFLOBDDDontCareNode
//********************************************************************

// Constructors/Destructor -------------------------------------------

// Default constructor
template <typename T, typename Op>
WeightedCFLOBDDDontCareNode<T, Op>::WeightedCFLOBDDDontCareNode(bool flag)
  :  WeightedCFLOBDDLeafNode<T, Op>(flag)
{
  this->numExits = 1;
  this->numPathsToExit = new long double[1];
  this->numPathsToExit[0] = 1;
}

template <typename T, typename Op>
WeightedCFLOBDDDontCareNode<T, Op>::WeightedCFLOBDDDontCareNode(T lweight, T rweight)
  :  WeightedCFLOBDDLeafNode<T, Op>(lweight, rweight)
{
  this->numExits = 1;
  this->numPathsToExit = new long double[1];
  this->numPathsToExit[0] = 1;
}

template <typename T, typename Op>
WeightedCFLOBDDDontCareNode<T,Op>::~WeightedCFLOBDDDontCareNode()
{
    if (this->isWeightsOfPathsMemAllocated)
        delete [] this->numWeightsOfPathsAsAmpsToExit;
}

template <typename T, typename Op>
void WeightedCFLOBDDDontCareNode<T,Op>::InstallWeightsOfPathsAsAmpsToExits()
{
    this->isWeightsOfPathsMemAllocated = true;
    this->numWeightsOfPathsAsAmpsToExit = new long double[1];
    this->numWeightsOfPathsAsAmpsToExit[0] = computeProbabilityFromAmplitude<T,Op>(this->lweight) + computeProbabilityFromAmplitude<T,Op>(this->rweight);
}

// print
template <typename T, typename Op>
std::ostream& WeightedCFLOBDDDontCareNode<T,Op>::print(std::ostream & out) const
{
  for (unsigned int i = this->level; i < this->maxLevel; i++) {
    out << "  ";
  }
  out << "Don't care (" << this->lweight << "," << this->rweight << ")";
  return out;
}

// template <typename T, typename Op>
// bool WeightedCFLOBDDDontCareNode<T,Op>::IsValid()
// {
// 	if (level == 0)
// 		isValid = true;
// 	return (level == 0);
// }

template <typename T, typename Op>
void WeightedCFLOBDDDontCareNode<T,Op>::FillSatisfyingAssignment(unsigned int, SH_OBDD::Assignment &assignment, unsigned int &index)
{
  index--;
  assignment[index] = 0;
}

template <typename T, typename Op>
int WeightedCFLOBDDDontCareNode<T,Op>::Traverse(SH_OBDD::AssignmentIterator &ai)
{
  ai.Next();
  return 0;
}

template <typename T, typename Op>
std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> WeightedCFLOBDDDontCareNode<T,Op>::Reduce(ReductionMapHandle&, unsigned int, WeightedValuesListHandle<T>& valList, bool)
{
    if (valList[0] == valList[1] && valList[0] == getIdentityValue<T,Op>())
        return std::make_pair(WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDDontCareNodeHandle, getIdentityValue<T,Op>());
    else if (valList[0] == valList[1] && valList[0] == getAnnhilatorValue<T,Op>())
        return std::make_pair(WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[0], valList[0]);
    else if (valList[0] == valList[1])
        return std::make_pair(this, valList[0]);
    else{
        T lw = valList[0] * this->lweight;
        T rw = valList[1] * this->rweight;
        std::tuple<T,T,T> w = computeInverseValue<T, Op>(lw, rw);
        if (std::get<1>(w) == getIdentityValue<T,Op>() && std::get<2>(w) == getIdentityValue<T,Op>())
            return std::make_pair(WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDDontCareNodeHandle, std::get<0>(w));
        WeightedCFLOBDDDontCareNode<T,Op>* n = new WeightedCFLOBDDDontCareNode<T,Op>(std::get<1>(w), std::get<2>(w));
        return std::make_pair(WeightedCFLOBDDNodeHandleT(n), std::get<0>(w));
    }
}

template <typename T, typename Op>
unsigned int WeightedCFLOBDDDontCareNode<T,Op>::Hash(unsigned int modsize)
{
    boost::hash<T> boost_hash;
    return (997 * boost_hash(this->lweight) + 97 * boost_hash(this->rweight) + 1) % modsize;
//   return (((unsigned int) reinterpret_cast<uintptr_t>(this) >> 2) + 997 * boost_hash(this->lweight) + 97 * boost_hash(this->rweight)) % modsize;
}

// Overloaded !=
template <typename T, typename Op>
bool WeightedCFLOBDDDontCareNode<T,Op>::operator!= (const WeightedCFLOBDDNode<T,Op> & n)
{
  return n.NodeKind() != W_CFLOBDD_DONTCARE;
}

// Overloaded ==
template <typename T, typename Op>
bool WeightedCFLOBDDDontCareNode<T,Op>::operator== (const WeightedCFLOBDDNode<T,Op> & n)
{
  return n.NodeKind() == W_CFLOBDD_DONTCARE;
}

// Restrict -----------------------------------------------------------
namespace CFL_OBDD {
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T,Op> Restrict(WeightedCFLOBDDNodeHandleT<T,Op> g, unsigned int i, bool val,
                           CFLOBDDReturnMapHandle &MapHandle
                          )
{
    WeightedCFLOBDDNodeHandleT<T,Op> answer;
  
    if (g.handleContents->NodeKind() == W_CFLOBDD_INTERNAL) {
      answer = Restrict((WeightedCFLOBDDInternalNode<T,Op> *)g.handleContents, i, val,
                        MapHandle
                       );
    }
    else if (g.handleContents->NodeKind() == W_CFLOBDD_FORK) {
      if (val == false) {
        MapHandle.AddToEnd(0);
        MapHandle.Canonicalize();
        answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDDontCareNodeHandle;
      }
      else { /* val == true */
        MapHandle.AddToEnd(1);
        MapHandle.Canonicalize();
        answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDDontCareNodeHandle;
      }
    }
    else { /* g.handleContents->NodeKind() == W_CFLOBDD_DONTCARE */
      MapHandle.AddToEnd(0);
      MapHandle.Canonicalize();
      answer = g;
    }
    return answer;
}

// TODO: Check this
template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T,Op> Restrict(WeightedCFLOBDDInternalNode<T,Op> *g, unsigned int i, bool val,
                           CFLOBDDReturnMapHandle &MapHandle
                          )
{
  if (g == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[g->level].handleContents) {
    MapHandle.AddToEnd(0);
    MapHandle.Canonicalize();
    return WeightedCFLOBDDNodeHandleT<T,Op>(g);
  }

  CFLOBDDReturnMapHandle AMap;
  WeightedCFLOBDDInternalNode<T,Op> *n = new WeightedCFLOBDDInternalNode<T,Op>(g->level);
  unsigned int j;
  unsigned int curExit;
  int b;

  if (i < (unsigned int)(1 << (g->level-1))) { // i falls in AConnection range
  	WeightedCFLOBDDNodeHandleT<T,Op> aHandle = Restrict(*(g->AConnection.entryPointHandle), i, val, AMap);
    n->AConnection.entryPointHandle = &aHandle;
    for (unsigned int k = 0; k < AMap.Size(); k++) {
      n->AConnection.returnMapHandle.AddToEnd(k);
    }
    n->AConnection.returnMapHandle.Canonicalize();
    j = 0;
    curExit = 0;
    n->numBConnections = AMap.Size();
    n->BConnection = new WConnection<T,Op>[n->numBConnections];
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
    n->BConnection = new WConnection<T,Op>[g->numBConnections];   // May create shorter version later
    n->numBConnections = 0;
    for (j = 0; j < g->numBConnections; j++) { // Perform a Restrict for each middle vertex
      CFLOBDDReturnMapHandle BMap;
      WeightedCFLOBDDNodeHandleT<T,Op> m = Restrict(*(g->BConnection[j].entryPointHandle), i-(unsigned int)(1 << (g->level-1)), val, BMap);

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
      WConnection<T,Op> candidate(m, inducedReturnMapHandleB);
      unsigned int position = n->InsertBConnection(n->numBConnections, candidate);
      AReductionMapHandle.AddToEnd(position);
    }
    AReductionMapHandle.Canonicalize();

    if (n->numBConnections < g->numBConnections) {  // Shorten
      WConnection<T,Op> *temp = n->BConnection;
      n->BConnection = new WConnection<T,Op>[n->numBConnections];
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
       WeightedCFLOBDDNodeHandleT<T,Op> tempHandle = g->AConnection.entryPointHandle->Reduce(inducedAReductionMapHandle, inducedAReturnMap.Size());
       n->AConnection = WConnection<T,Op>(tempHandle, inducedAReturnMap);
  }

  n->numExits = curExit;
#ifdef PATH_COUNTING_ENABLED
  n->InstallPathCounts();
#endif
  MapHandle.Canonicalize();
  return WeightedCFLOBDDNodeHandleT<T,Op>(n);
}

template <typename T, typename Op>
WeightedBDDTopNode<T,Op>::WeightedBDDTopNode(unsigned int numberOfVars) : numberOfVars(numberOfVars)
{
}

template <typename T, typename Op>
WeightedBDDTopNode<T,Op>::~WeightedBDDTopNode()
{
}

template <typename T, typename Op>
bool WeightedBDDTopNode<T,Op>::operator!= (const WeightedCFLOBDDNode<T,Op> &n)
{
  if (n.NodeKind() != W_BDD_TOPNODE)
    return true;
  WeightedBDDTopNode& nh = (WeightedBDDTopNode&)n;
  return bddContents != nh.bddContents;
}

template <typename T, typename Op>
bool WeightedBDDTopNode<T,Op>::operator== (const WeightedCFLOBDDNode<T,Op> &n)
{
  if (n.NodeKind() != W_BDD_TOPNODE)
    return false;
  WeightedBDDTopNode& nh = (WeightedBDDTopNode&)n;
  return bddContents == nh.bddContents;
}
template <typename T, typename Op>
void WeightedBDDTopNode<T,Op>::IncrRef()
{
  this->refCount++;
}
template <typename T, typename Op>
void WeightedBDDTopNode<T,Op>::DecrRef()
{
  if (--this->refCount == 0) {    // Warning: Saturation not checked
    if (this->isCanonical) {
      WeightedCFLOBDDNodeHandleT<T,Op>::canonicalNodeTable->DeleteEq(this);
    }
    delete this;
  }  
}

template <typename T, typename Op>
unsigned int WeightedBDDTopNode<T,Op>::Hash(unsigned int modsize)
{
  return (117 * bddContents.Hash(modsize) + 1) % modsize;
}

template <typename T, typename Op>
std::ostream& WeightedBDDTopNode<T,Op>::print(std::ostream &out) const
{
  bddContents.print(out);
  return (out);
}

template <typename T, typename Op>
void WeightedBDDTopNode<T,Op>::FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index) {}
template <typename T, typename Op>
int WeightedBDDTopNode<T,Op>::Traverse(SH_OBDD::AssignmentIterator &ai) {}
template <typename T, typename Op>
std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> WeightedBDDTopNode<T,Op>::Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valList, bool forceReduce) {}
template <typename T, typename Op>
void WeightedBDDTopNode<T,Op>::DumpConnections(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visited, std::ostream & out) {}
template <typename T, typename Op>
void WeightedBDDTopNode<T,Op>::CountNodesAndEdges(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
unsigned int &nodeCount, unsigned int &edgeCount, unsigned int& returnEdgesCount) {}
template <typename T, typename Op>
void WeightedBDDTopNode<T,Op>::CountNodes(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, unsigned int &nodeCount) {}
template <typename T, typename Op>
void WeightedBDDTopNode<T,Op>::CountPaths(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes) {}
template <typename T, typename Op>
void WeightedBDDTopNode<T,Op>::ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes) {
  bddContents.ComputeWeightofPathsAsAmpsToExits();
}

}
