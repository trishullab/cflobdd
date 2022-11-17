#ifndef CFLOBDD_NODE_GUARD
#define CFLOBDD_NODE_GUARD

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


// Configuration flags --------------------------------------------
//#define PATH_COUNTING_ENABLED 1

namespace CFL_OBDD {

// Node classes declared in this file --------------------------------
class CFLOBDDNode;
class CFLOBDDInternalNode;   //  : public CFLOBDDNode
class CFLOBDDLeafNode;       //  : public CFLOBDDNode
class CFLOBDDForkNode;       //  : public CFLOBDDLeafNode
class CFLOBDDDontCareNode;   //  : public CFLOBDDLeafNode
class CFLOBDDNodeHandle;
}

#include <iostream>
#include <fstream>
//#include <mpirxx.h>
//#include <boost/multiprecision/cpp_int.hpp>

#include "list_T.h"
#include "list_TPtr.h"
#include "intpair.h"
#include "conscell.h"
#include "traverse_state_cfl.h"
#include "assignment.h"
#include "bool_op.h"
#include "return_map_T.h"

//using namespace boost::multiprecision;

namespace CFL_OBDD {
	typedef ReturnMapBody<int> CFLOBDDReturnMapBody;
	typedef ReturnMapHandle<int> CFLOBDDReturnMapHandle;
}
#include "connectionT.h"
namespace CFL_OBDD {
	typedef ConnectionT<CFLOBDDReturnMapHandle> Connection;
}
#include "reduction_map.h"
#include "hash.h"
#include "hashset.h"
#include "ref_ptr.h"



namespace CFL_OBDD {

//********************************************************************
// CFLOBDDNodeHandle
//********************************************************************

class CFLOBDDNodeHandle {
#define CFLOBDD_NODE_HANDLE_GUARD
	// friend CFLOBDDNodeHandle ApplyAndReduce(CFLOBDDNodeHandle n1, CFLOBDDNodeHandle n2, BoolOp op);
	// friend CFLOBDDNodeHandle ApplyAndReduce(CFLOBDDNodeHandle n1, CFLOBDDNodeHandle n2, CFLOBDDNodeHandle n3, BoolOp3 op);
 public:
  CFLOBDDNodeHandle();                                        // Default constructor
  CFLOBDDNodeHandle(CFLOBDDNode *n);                          // Constructor
  CFLOBDDNodeHandle(const CFLOBDDNodeHandle &nh);              // Copy constructor
  ~CFLOBDDNodeHandle();                                       // Destructor
  unsigned int Hash(unsigned int modsize);
  bool operator!= (const CFLOBDDNodeHandle &nh);              // Overloaded !=
  bool operator== (const CFLOBDDNodeHandle &nh) const;              // Overloaded ==
  CFLOBDDNodeHandle & operator= (const CFLOBDDNodeHandle &nh); // assignment

  // Distinguished CFLOBDDNodeHandles -----------------------
     public:
      static CFLOBDDNodeHandle *NoDistinctionNode;        // CFLOBDDNodeHandle NoDistinctionNode[maxLevel+1] for levels 0, ..., CFLOBDDMaxLevel
      static CFLOBDDNodeHandle CFLOBDDForkNodeHandle;
      static CFLOBDDNodeHandle CFLOBDDDontCareNodeHandle;
      static void InitNoDistinctionTable();

	  static CFLOBDDNodeHandle *AdditionInterleavedNoCarryNode;  // CFLOBDDNodeHandle AdditionInterleavedNoCarryNode[maxLevel-1] for levels 2, ..., CFLOBDDMaxLevel
	  static CFLOBDDNodeHandle *AdditionInterleavedCarryNode;    // CFLOBDDNodeHandle AdditionInterleavedCarryNode[maxLevel-1] for levels 2, ..., CFLOBDDMaxLevel
	  static void InitAdditionInterleavedTable();
	  
  // The data member
     CFLOBDDNode *handleContents;

 // Table of canonical nodes -------------------------
    public:
     static Hashset<CFLOBDDNode> *canonicalNodeTable;
     void Canonicalize();

 // Reduce and its associated cache ---------------
    public:
     CFLOBDDNodeHandle Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, bool forceReduce = false);
     static void InitReduceCache();
     static void DisposeOfReduceCache();

	public:
	 std::ostream& print(std::ostream & out = std::cout) const;
	 bool IsValid();
	 struct CFLOBDDNodeHandle_Hash {
	 public:
		 size_t operator()(const CFLOBDDNodeHandle& c) const {
			 return ((reinterpret_cast<std::uintptr_t>(c.handleContents) >> 2) % 997);
		 }
	 };
};

std::ostream& operator<< (std::ostream & out, const CFLOBDDNodeHandle &d);


typedef Hashtable<CFLOBDDNodeHandle, CFLOBDDNodeHandle> CFLOBDDNodeMemoTable;
typedef ref_ptr<CFLOBDDNodeMemoTable> CFLOBDDNodeMemoTableRefPtr;

extern CFLOBDDNodeHandle MkDistinction(unsigned int level, unsigned int i);
extern CFLOBDDNodeHandle Restrict(CFLOBDDNodeHandle g, unsigned int i, bool val, CFLOBDDReturnMapHandle &MapHandle);
extern CFLOBDDNodeHandle Restrict(CFLOBDDInternalNode *g, unsigned int i, bool val, CFLOBDDReturnMapHandle &MapHandle);
extern CFLOBDDNodeHandle duplicateAinBAtLevelOne(CFLOBDDNodeHandle f);
extern CFLOBDDNodeHandle MkAdditionInterleavedRecursive(unsigned int level, bool carry);
extern CFLOBDDNodeHandle MkParity(unsigned int level);
extern CFLOBDDNodeHandle MkStepOneFourth(unsigned int level);
extern CFLOBDDNodeHandle shiftAtoB(CFLOBDDNodeHandle f, const unsigned int levelAtWhichToShift);
extern CFLOBDDNodeHandle shiftAtoBAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, CFLOBDDNodeHandle f);
extern CFLOBDDNodeHandle shiftBtoA(CFLOBDDNodeHandle f, const unsigned int levelAtWhichToShift);
extern CFLOBDDNodeHandle shiftBtoAAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, CFLOBDDNodeHandle f);
#ifdef ARBITRARY_STEP_FUNCTIONS
extern CFLOBDDNodeHandle MkStepNode(unsigned int level, unsigned int left, unsigned int middle, unsigned int right);
#endif

extern double ComputeProbabilityNode(CFLOBDDNodeHandle g, std::vector<double>& var_probs, std::vector<double>& path_probs, int start, int end); 
extern std::vector<double> ComputeProbabilityOfListNode(CFLOBDDNodeHandle g, std::vector<std::vector<double>>& var_probs, std::vector<std::vector<double>>& path_probs, int start, int end); 
// extern std::vector<double> ComputeEntropyOfListNode(CFLOBDDNodeHandle g, std::vector<std::vector<double>>& var_probs, std::vector<std::vector<double>>& path_probs, std::vector<std::vector<double>>& entropy, int start, int end); 


}

#include "cross_product.h"

namespace CFL_OBDD {
	// Other classes and types declared in this file ---------------------
	typedef List<intpair> RestrictMap;
	typedef ListIterator<intpair> RestrictMapIterator;



	// Auxiliary data -----------------------------------------------

	// Note: If CFLOBDDMaxLevel >= 27, allocating an Assignment may cause
	//       virtual memory to be exceeded.
#define CFLOBDD_MAX_LEVEL 30
	unsigned int const CFLOBDDMaxLevel = CFLOBDD_MAX_LEVEL;

	//********************************************************************
	// CFLReduceKey
	//********************************************************************

	class CFLReduceKey {

	public:
		CFLReduceKey(CFLOBDDNodeHandle nodeHandle, ReductionMapHandle redMap); // Constructor
		unsigned int Hash(unsigned int modsize);
		CFLReduceKey& operator= (const CFLReduceKey& p);  // Overloaded assignment
		bool operator!= (const CFLReduceKey& p);        // Overloaded !=
		bool operator== (const CFLReduceKey& p);        // Overloaded ==
		CFLOBDDNodeHandle NodeHandle() const { return nodeHandle; }      // Access function
		ReductionMapHandle RedMapHandle() const { return redMapHandle; } // Access function
		std::ostream& print(std::ostream & out) const;

	private:
		CFLOBDDNodeHandle nodeHandle;
		ReductionMapHandle redMapHandle;
		CFLReduceKey();                                 // Default constructor (hidden)
	};

	std::ostream& operator<< (std::ostream & out, const CFLReduceKey &p);

}

#include "cflobdd_top_node_t.h"

//********************************************************************
// CFLOBDDNode
//********************************************************************

namespace CFL_OBDD {

enum CFLOBDD_NODEKIND { CFLOBDD_INTERNAL, CFLOBDD_FORK, CFLOBDD_DONTCARE };

class CFLOBDDNode {
  friend void CFLOBDDNodeHandle::InitNoDistinctionTable();
  friend void CFLOBDDNodeHandle::InitAdditionInterleavedTable();
 public:
  CFLOBDDNode();                       // Constructor
  CFLOBDDNode(const unsigned int l);   // Constructor
  virtual ~CFLOBDDNode();              // Destructor
  virtual CFLOBDD_NODEKIND NodeKind() const = 0;
  unsigned int numExits;
  static unsigned int const maxLevel;
//#ifdef PATH_COUNTING_ENABLED
  //unsigned long long int *numPathsToExit;       // unsigned int numPathsToExit[numExits]
  //cpp_int *numPathsToExit;
  long double *numPathsToExit;
  bool isNumPathsMemAllocated;
//#endif
  virtual void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index) = 0;
  virtual int Traverse(SH_OBDD::AssignmentIterator &ai) = 0;
  virtual CFLOBDDNodeHandle Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, bool forceReduce = false) = 0;
  virtual unsigned int Hash(unsigned int modsize) = 0;
  virtual void DumpConnections(Hashset<CFLOBDDNodeHandle> *visited, std::ostream & out = std::cout) = 0;
  virtual void CountNodesAndEdges(Hashset<CFLOBDDNodeHandle> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
	  unsigned int &nodeCount, unsigned int &edgeCount, unsigned int &returnEdgesCount) = 0;
  virtual void CountNodes(Hashset<CFLOBDDNodeHandle> *visitedNodes, unsigned int &nodeCount) = 0;
  virtual void CountPaths(Hashset<CFLOBDDNodeHandle> *visitedNodes) = 0;
  virtual bool operator!= (const CFLOBDDNode & n) = 0;  // Overloaded !=
  virtual bool operator== (const CFLOBDDNode & n) = 0;  // Overloaded ==
  virtual void IncrRef() = 0;
  virtual void DecrRef() = 0;
  const unsigned int Level() const { return level; }
  const bool IsCanonical() const { return isCanonical; }
  void SetCanonical() { isCanonical = true;  }
  unsigned int GetRefCount(){ return refCount; }
 public:
  virtual std::ostream& print(std::ostream & out = std::cout) const = 0;
  virtual bool IsValid() = 0;
  const unsigned int level;
 protected:
  unsigned int refCount;
  bool isCanonical;              // Is this CFLOBDDNode in canonicalNodeTable?
  bool isValid;
};

std::ostream& operator<< (std::ostream & out, const CFLOBDDNode &n);

//********************************************************************
// CFLOBDDInternalNode
//********************************************************************

class CFLOBDDInternalNode : public CFLOBDDNode {
  friend void CFLOBDDNodeHandle::InitNoDistinctionTable();
  friend void CFLOBDDNodeHandle::InitAdditionInterleavedTable();
  friend CFLOBDDNodeHandle Restrict(CFLOBDDInternalNode *g, unsigned int i, bool val, CFLOBDDReturnMapHandle &MapHandle);
  friend CFLOBDDNodeHandle PairProduct(CFLOBDDInternalNode *n1, CFLOBDDInternalNode *n2, PairProductMapHandle &pairProductMap);
  friend CFLOBDDNodeHandle TripleProduct(CFLOBDDInternalNode *n1, CFLOBDDInternalNode *n2, CFLOBDDInternalNode *n3, TripleProductMapHandle &tripleProductMap);
//  friend bool CFLOBDDTopNode::EvaluateIteratively(SH_OBDD::Assignment &assignment);
//  friend void CFLOBDDTopNode::PrintYieldAux(std::ostream * out, List<ConsCell<TraverseState> *> &T, ConsCell<TraverseState> *S);
#ifdef ARBITRARY_STEP_FUNCTIONS
  friend CFLOBDDNodeHandle MkStepNode(unsigned int level, unsigned int left, unsigned int middle, unsigned int right);
#endif

 public:
  CFLOBDDInternalNode(const unsigned int l);   // Constructor
  ~CFLOBDDInternalNode();                      // Destructor
  CFLOBDD_NODEKIND NodeKind() const { return CFLOBDD_INTERNAL; }
  void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index);
  int Traverse(SH_OBDD::AssignmentIterator &ai);
  CFLOBDDNodeHandle Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, bool forceReduce = false);
  unsigned int Hash(unsigned int modsize);
  void DumpConnections(Hashset<CFLOBDDNodeHandle> *visited, std::ostream & out = std::cout);
  void CountNodesAndEdges(Hashset<CFLOBDDNodeHandle> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
	  unsigned int &nodeCount, unsigned int &edgeCount, unsigned int& returnEdgesCount);
  void CountNodes(Hashset<CFLOBDDNodeHandle> *visitedNodes, unsigned int &nodeCount);
  void CountPaths(Hashset<CFLOBDDNodeHandle> *visitedNodes);
  bool operator!= (const CFLOBDDNode & n);        // Overloaded !=
  bool operator== (const CFLOBDDNode & n);        // Overloaded ==
  void IncrRef();
  void DecrRef();

 public:
  std::ostream& print(std::ostream & out = std::cout) const;
  bool IsValid();
  Connection AConnection;              // A single Connection
  unsigned int numBConnections;
  Connection *BConnection;             // Connection BConnection[numBConnections];
  unsigned int InsertBConnection(unsigned int &j, Connection &c);
//#ifdef PATH_COUNTING_ENABLED
  void InstallPathCounts();
//#endif

 private:
  CFLOBDDInternalNode();                                         // Default constructor (hidden)
  CFLOBDDInternalNode(const CFLOBDDInternalNode &n);             // Copy constructor (hidden)
  CFLOBDDInternalNode& operator= (const CFLOBDDInternalNode &n); // Overloaded = (hidden)
};

//********************************************************************

// std::ostream& operator<< (std::ostream & out, const CFLOBDDTopNode &d);

//********************************************************************
// CFLOBDDLeafNode
//********************************************************************

class CFLOBDDLeafNode : public CFLOBDDNode {
 public:
  CFLOBDDLeafNode();                   // Constructor
  virtual ~CFLOBDDLeafNode();          // Destructor
  virtual CFLOBDD_NODEKIND NodeKind() const = 0;
  virtual void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index) = 0;
  virtual int Traverse(SH_OBDD::AssignmentIterator &ai) = 0;
  virtual CFLOBDDNodeHandle Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, bool forceReduce = false) = 0;
  virtual unsigned int Hash(unsigned int modsize) = 0;
  void DumpConnections(Hashset<CFLOBDDNodeHandle> *visited, std::ostream & out = std::cout);
  void CountNodesAndEdges(Hashset<CFLOBDDNodeHandle> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
	  unsigned int &nodeCount, unsigned int &edgeCount, unsigned int& returnEdgesCount);
  void CountNodes(Hashset<CFLOBDDNodeHandle> *visitedNodes, unsigned int &nodeCount);
  void CountPaths(Hashset<CFLOBDDNodeHandle> *visitedNodes);
  virtual bool operator!= (const CFLOBDDNode & n) = 0;  // Overloaded !=
  virtual bool operator== (const CFLOBDDNode & n) = 0;  // Overloaded ==
  void IncrRef();
  void DecrRef();

 public:
	 virtual std::ostream& print(std::ostream & out = std::cout) const = 0;
	 virtual bool IsValid() = 0;
};

//********************************************************************
// CFLOBDDForkNode
//********************************************************************

class CFLOBDDForkNode : public CFLOBDDLeafNode {
 public:
  CFLOBDDForkNode();                   // Constructor
  ~CFLOBDDForkNode();                  // Destructor
  CFLOBDD_NODEKIND NodeKind() const { return CFLOBDD_FORK; }
  void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index);
  int Traverse(SH_OBDD::AssignmentIterator &ai);
  CFLOBDDNodeHandle Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, bool forceReduce = false);
  unsigned int Hash(unsigned int modsize);
  bool operator!= (const CFLOBDDNode & n);        // Overloaded !=
  bool operator== (const CFLOBDDNode & n);        // Overloaded ==

 public:
	 std::ostream& print(std::ostream & out = std::cout) const;
	 bool IsValid();

 private:
  CFLOBDDForkNode(const CFLOBDDForkNode &n);   // Copy constructor (hidden)
  CFLOBDDForkNode& operator= (const CFLOBDDForkNode &n); // Overloaded = (hidden)
};

//********************************************************************
// CFLOBDDDontCareNode
//********************************************************************

class CFLOBDDDontCareNode : public CFLOBDDLeafNode {
 public:
  CFLOBDDDontCareNode();                   // Constructor
  ~CFLOBDDDontCareNode();                  // Destructor
  CFLOBDD_NODEKIND NodeKind() const { return CFLOBDD_DONTCARE; }
  void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index);
  int Traverse(SH_OBDD::AssignmentIterator &ai);
  CFLOBDDNodeHandle Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, bool forceReduce = false);
  unsigned int Hash(unsigned int modsize);
  bool operator!= (const CFLOBDDNode & n);        // Overloaded !=
  bool operator== (const CFLOBDDNode & n);        // Overloaded ==

 public:
	 std::ostream& print(std::ostream & out = std::cout) const;
	 bool IsValid();

 private:
  CFLOBDDDontCareNode(const CFLOBDDDontCareNode &n);   // Copy constructor (hidden)
  CFLOBDDDontCareNode& operator= (const CFLOBDDDontCareNode &n); // Overloaded = (hidden)
};

 } // namespace CFL_OBDD

#endif
