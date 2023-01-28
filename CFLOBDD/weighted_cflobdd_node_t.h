#ifndef WEIGHTED_CFLOBDD_NODE_T_GUARD
#define WEIGHTED_CFLOBDD_NODE_T_GUARD



// Configuration flags --------------------------------------------
//#define PATH_COUNTING_ENABLED 1

namespace CFL_OBDD {

// Node classes declared in this file --------------------------------
template <typename, typename>
class WeightedCFLOBDDNode;
template <typename, typename>
class WeightedCFLOBDDInternalNode;   //  : public CFLOBDDNode
template <typename, typename>
class WeightedCFLOBDDLeafNode;       //  : public CFLOBDDNode
template <typename, typename>
class WeightedCFLOBDDForkNode;       //  : public CFLOBDDLeafNode
template <typename, typename>
class WeightedCFLOBDDDontCareNode;   //  : public CFLOBDDLeafNode
template <typename, typename>
class WeightedCFLOBDDNodeHandleT;
template <typename, typename>
class WeightedBDDTopNode;
}

#include <iostream>
#include <fstream>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "list_T.h"
#include "list_TPtr.h"
#include "intpair.h"
// #include "conscell.h"
#include "assignment.h"
#include "bool_op.h"
#include "return_map_T.h"

//using namespace boost::multiprecision;

namespace CFL_OBDD {
	typedef ReturnMapBody<int> CFLOBDDReturnMapBody;
	typedef ReturnMapHandle<int> CFLOBDDReturnMapHandle;
}
#include "weighted_connectionT.h"
namespace CFL_OBDD {
	// template <typename T, typename Op> 
  // using WConnection = WeightedConnectionT<T, Op> ;
  // typedef WeightedConnectionT WConnection
}
#include "reduction_map.h"
#include "hash.h"
#include "hashset.h"
#include "ref_ptr.h"
#include "weighted_values_list.h"
#include "weighted_bdd_node_t.h"



namespace CFL_OBDD {

//********************************************************************
// WeightedCFLOBDDNodeHandleT
//********************************************************************
template<typename T, typename Op>
class WeightedCFLOBDDNodeHandleT {
#define CFLOBDD_NODE_HANDLE_GUARD
	// friend WeightedCFLOBDDNodeHandleT ApplyAndReduce(WeightedCFLOBDDNodeHandleT n1, WeightedCFLOBDDNodeHandleT n2, BoolOp op);
	// friend WeightedCFLOBDDNodeHandleT ApplyAndReduce(WeightedCFLOBDDNodeHandleT n1, WeightedCFLOBDDNodeHandleT n2, WeightedCFLOBDDNodeHandleT n3, BoolOp3 op);
 public:
  WeightedCFLOBDDNodeHandleT();                                        // Default constructor
  WeightedCFLOBDDNodeHandleT(WeightedCFLOBDDNode<T,Op> *n);                          // Constructor
  WeightedCFLOBDDNodeHandleT(const WeightedCFLOBDDNodeHandleT<T,Op> &nh);              // Copy constructor
  ~WeightedCFLOBDDNodeHandleT();                                       // Destructor
  unsigned int Hash(unsigned int modsize);
  bool operator!= (const WeightedCFLOBDDNodeHandleT &nh);              // Overloaded !=
  bool operator== (const WeightedCFLOBDDNodeHandleT &nh) const;              // Overloaded ==
  WeightedCFLOBDDNodeHandleT<T,Op> & operator= (const WeightedCFLOBDDNodeHandleT<T,Op> &nh); // assignment

  // Distinguished WeightedCFLOBDDNodeHandleTs -----------------------
     public:
      static WeightedCFLOBDDNodeHandleT *NoDistinctionNode;        // WeightedCFLOBDDNodeHandleT NoDistinctionNode[maxLevel+1] for levels 0, ..., CFLOBDDMaxLevel
      static WeightedCFLOBDDNodeHandleT *NoDistinctionNode_Ann;        // WeightedCFLOBDDNodeHandleT NoDistinctionNode_Ann[maxLevel+1] for levels 0, ..., CFLOBDDMaxLevel
      static WeightedCFLOBDDNodeHandleT *IdentityNode;        // WeightedCFLOBDDNodeHandleT IdentityNode[maxLevel+1] for levels 1, ..., CFLOBDDMaxLevel
      static WeightedCFLOBDDNodeHandleT CFLOBDDForkNodeHandle;
      static WeightedCFLOBDDNodeHandleT CFLOBDDForkNodeHandle01;
      static WeightedCFLOBDDNodeHandleT CFLOBDDForkNodeHandle10;
      static WeightedCFLOBDDNodeHandleT CFLOBDDDontCareNodeHandle;
      static void InitNoDistinctionTable();
      static void InitNoDistinctionTable_Ann();
      static void InitIdentityNodeTable();

	  static WeightedCFLOBDDNodeHandleT *AdditionInterleavedNoCarryNode;  // WeightedCFLOBDDNodeHandleT AdditionInterleavedNoCarryNode[maxLevel-1] for levels 2, ..., CFLOBDDMaxLevel
	  static WeightedCFLOBDDNodeHandleT *AdditionInterleavedCarryNode;    // WeightedCFLOBDDNodeHandleT AdditionInterleavedCarryNode[maxLevel-1] for levels 2, ..., CFLOBDDMaxLevel
	  static void InitAdditionInterleavedTable();
	  
  // The data member
    WeightedCFLOBDDNode<T, Op> *handleContents;

 // Table of canonical nodes -------------------------
    public:
     static Hashset<WeightedCFLOBDDNode<T, Op>> *canonicalNodeTable;
     void Canonicalize();

 // Reduce and its associated cache ---------------
    public:
     std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valueTuple, bool forceReduce = false);
     static void InitReduceCache();
     static void DisposeOfReduceCache();

	public:
	 std::ostream& print(std::ostream & out = std::cout) const;
	//  bool IsValid();
	 struct WeightedCFLOBDDNodeHandleT_Hash {
	 public:
		 size_t operator()(const WeightedCFLOBDDNodeHandleT<T,Op>& c) const {
			 return ((reinterpret_cast<std::uintptr_t>(c.handleContents) >> 2) % 997);
		 }
	 };
};

template <typename T, typename Op>
std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDNodeHandleT<T, Op> &d);



template <typename T, typename Op> 
using WeightedCFLOBDDNodeMemoTable = Hashtable<WeightedCFLOBDDNodeHandleT<T, Op>, WeightedCFLOBDDNodeHandleT<T, Op>>;
template <typename T, typename Op> 
using WeightedCFLOBDDNodeMemoTableRefPtr = ref_ptr<WeightedCFLOBDDNodeMemoTable<T,Op>>;

// extern WeightedCFLOBDDNodeHandleT MkDistinction(unsigned int level, unsigned int i);
// extern WeightedCFLOBDDNodeHandleT Restrict(WeightedCFLOBDDNodeHandleT g, unsigned int i, bool val, CFLOBDDReturnMapHandle &MapHandle);
// extern WeightedCFLOBDDNodeHandleT Restrict(CFLOBDDInternalNode *g, unsigned int i, bool val, CFLOBDDReturnMapHandle &MapHandle);
// extern WeightedCFLOBDDNodeHandleT duplicateAinBAtLevelOne(WeightedCFLOBDDNodeHandleT f);
// extern WeightedCFLOBDDNodeHandleT MkAdditionInterleavedRecursive(unsigned int level, bool carry);
// extern WeightedCFLOBDDNodeHandleT MkParity(unsigned int level);
// extern WeightedCFLOBDDNodeHandleT MkStepOneFourth(unsigned int level);
// extern WeightedCFLOBDDNodeHandleT shiftAtoB(WeightedCFLOBDDNodeHandleT f, const unsigned int levelAtWhichToShift);
// extern WeightedCFLOBDDNodeHandleT shiftAtoBAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, WeightedCFLOBDDNodeHandleT f);
// extern WeightedCFLOBDDNodeHandleT shiftBtoA(WeightedCFLOBDDNodeHandleT f, const unsigned int levelAtWhichToShift);
// extern WeightedCFLOBDDNodeHandleT shiftBtoAAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, WeightedCFLOBDDNodeHandleT f);
// #ifdef ARBITRARY_STEP_FUNCTIONS
// extern WeightedCFLOBDDNodeHandleT MkStepNode(unsigned int level, unsigned int left, unsigned int middle, unsigned int right);
// #endif

// extern double ComputeProbabilityNode(WeightedCFLOBDDNodeHandleT g, std::vector<double>& var_probs, std::vector<double>& path_probs, int start, int end); 
// extern std::vector<double> ComputeProbabilityOfListNode(WeightedCFLOBDDNodeHandleT g, std::vector<std::vector<double>>& var_probs, std::vector<std::vector<double>>& path_probs, int start, int end); 
// extern std::vector<double> ComputeEntropyOfListNode(WeightedCFLOBDDNodeHandleT g, std::vector<std::vector<double>>& var_probs, std::vector<std::vector<double>>& path_probs, std::vector<std::vector<double>>& entropy, int start, int end); 


}

// #include "cross_product.h"

namespace CFL_OBDD {
	// Other classes and types declared in this file ---------------------

  #define CFLOBDD_MAX_LEVEL 30
	  unsigned int const CFLOBDDMaxlevel = CFLOBDD_MAX_LEVEL;

	//********************************************************************
	// WeightedCFLReduceKey
	//********************************************************************

    template <typename T, typename Op>
	class WeightedCFLReduceKey {

	public:
		WeightedCFLReduceKey(WeightedCFLOBDDNodeHandleT<T, Op> nodeHandle, ReductionMapHandle redMap, WeightedValuesListHandle<T> valList ); // Constructor
		unsigned int Hash(unsigned int modsize);
		WeightedCFLReduceKey& operator= (const WeightedCFLReduceKey& p);  // Overloaded assignment
		bool operator!= (const WeightedCFLReduceKey& p);        // Overloaded !=
		bool operator== (const WeightedCFLReduceKey& p);        // Overloaded ==
		WeightedCFLOBDDNodeHandleT<T,Op> NodeHandle() const { return nodeHandle; }      // Access function
		ReductionMapHandle RedMapHandle() const { return redMapHandle; } // Access function
    WeightedValuesListHandle<T> ValList() const { return valList; }
		std::ostream& print(std::ostream & out) const;

	private:
		WeightedCFLOBDDNodeHandleT<T, Op> nodeHandle;
		ReductionMapHandle redMapHandle;
    WeightedValuesListHandle<T> valList;
		WeightedCFLReduceKey();                                 // Default constructor (hidden)
	};

    template <typename T, typename Op>
	std::ostream& operator<< (std::ostream & out, const WeightedCFLReduceKey<T, Op> &p);

}

#include "weighted_cflobdd_top_node_t.h"

//********************************************************************
// CFLOBDDNode
//********************************************************************

namespace CFL_OBDD {

// TODO: Change this during integration
enum WCFLOBDD_NODEKIND { W_CFLOBDD_INTERNAL, W_CFLOBDD_FORK, W_CFLOBDD_DONTCARE, W_BDD_TOPNODE };

template <typename T, typename Op>
class WeightedCFLOBDDNode {
  friend void WeightedCFLOBDDNodeHandleT<T,Op>::InitNoDistinctionTable();
  friend void WeightedCFLOBDDNodeHandleT<T,Op>::InitNoDistinctionTable_Ann();
  friend void WeightedCFLOBDDNodeHandleT<T,Op>::InitAdditionInterleavedTable();
 public:
  WeightedCFLOBDDNode();                       // Constructor
  WeightedCFLOBDDNode(const unsigned int l);   // Constructor
  virtual ~WeightedCFLOBDDNode();              // Destructor
  virtual WCFLOBDD_NODEKIND NodeKind() const = 0;
  unsigned int numExits;
  static unsigned int const maxLevel;
//#ifdef PATH_COUNTING_ENABLED
  //unsigned long long int *numPathsToExit;       // unsigned int numPathsToExit[numExits]
  //cpp_int *numPathsToExit;
  long double *numPathsToExit;
  long double *numWeightsOfPathsAsAmpsToExit;
  bool isNumPathsMemAllocated;
  bool isWeightsOfPathsMemAllocated;
//#endif
  virtual void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index) = 0;
  virtual int Traverse(SH_OBDD::AssignmentIterator &ai) = 0;
  virtual std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valList, bool forceReduce = false) = 0;
  virtual unsigned int Hash(unsigned int modsize) = 0;
  virtual void DumpConnections(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visited, std::ostream & out = std::cout) = 0;
  virtual void CountNodesAndEdges(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
        unsigned int &nodeCount, unsigned int &edgeCount, unsigned int &returnEdgesCount) = 0;
  virtual void CountNodes(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, unsigned int &nodeCount) = 0;
  virtual void CountPaths(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes) = 0;
  virtual void ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes) = 0;
  virtual bool operator!= (const WeightedCFLOBDDNode<T,Op> & n) = 0;  // Overloaded !=
  virtual bool operator== (const WeightedCFLOBDDNode<T,Op> & n) = 0;  // Overloaded ==
  virtual void IncrRef() = 0;
  virtual void DecrRef() = 0;
  const unsigned int Level() const { return level; }
  const bool IsCanonical() const { return isCanonical; }
  void SetCanonical() { isCanonical = true;  }
  unsigned int GetRefCount(){ return refCount; }
 public:
  virtual std::ostream& print(std::ostream & out = std::cout) const = 0;
  // virtual bool IsValid() = 0;
  const unsigned int level;
 protected:
  unsigned int refCount;
  bool isCanonical;              // Is this CFLOBDDNode in canonicalNodeTable?
  // bool isValid;
};

template <typename T, typename Op>
std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDNode<T, Op> &n);

//********************************************************************
// WeightedCFLOBDDInternalNode
//********************************************************************

template <typename T, typename Op>
class WeightedCFLOBDDInternalNode : public WeightedCFLOBDDNode<T,Op> {
  friend void WeightedCFLOBDDNodeHandleT<T, Op>::InitNoDistinctionTable();
  friend void WeightedCFLOBDDNodeHandleT<T, Op>::InitNoDistinctionTable_Ann();
  friend void WeightedCFLOBDDNodeHandleT<T, Op>::InitAdditionInterleavedTable();
  friend WeightedCFLOBDDNodeHandleT<T,Op> Restrict(WeightedCFLOBDDInternalNode *g, unsigned int i, bool val, CFLOBDDReturnMapHandle &MapHandle);
  // friend WeightedCFLOBDDNodeHandleT<T,Op> PairProduct(WeightedCFLOBDDInternalNode *n1, WeightedCFLOBDDInternalNode *n2, PairProductMapHandle &pairProductMap);
  // friend WeightedCFLOBDDNodeHandleT<T,Op> TripleProduct(WeightedCFLOBDDInternalNode *n1, WeightedCFLOBDDInternalNode *n2, WeightedCFLOBDDInternalNode *n3, TripleProductMapHandle &tripleProductMap);
//  friend bool CFLOBDDTopNode::EvaluateIteratively(SH_OBDD::Assignment &assignment);
//  friend void CFLOBDDTopNode::PrintYieldAux(std::ostream * out, List<ConsCell<TraverseState> *> &T, ConsCell<TraverseState> *S);
#ifdef ARBITRARY_STEP_FUNCTIONS
  friend WeightedCFLOBDDNodeHandleT MkStepNode(unsigned int level, unsigned int left, unsigned int middle, unsigned int right);
#endif

 public:
  WeightedCFLOBDDInternalNode(const unsigned int l);   // Constructor
  ~WeightedCFLOBDDInternalNode();                      // Destructor
  WCFLOBDD_NODEKIND NodeKind() const { return W_CFLOBDD_INTERNAL; }
  void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index);
  int Traverse(SH_OBDD::AssignmentIterator &ai);
  std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valList, bool forceReduce = false);
  unsigned int Hash(unsigned int modsize);
  void DumpConnections(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visited, std::ostream & out = std::cout);
  void CountNodesAndEdges(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
	  unsigned int &nodeCount, unsigned int &edgeCount, unsigned int& returnEdgesCount);
  void CountNodes(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, unsigned int &nodeCount);
  void CountPaths(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes);
  void ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes);
  bool operator!= (const WeightedCFLOBDDNode<T,Op> & n);        // Overloaded !=
  bool operator== (const WeightedCFLOBDDNode<T,Op> & n);        // Overloaded ==
  void IncrRef();
  void DecrRef();

 public:
  std::ostream& print(std::ostream & out = std::cout) const;
  // bool IsValid();
  WConnection<T,Op> AConnection;              // A single Connection
  unsigned int numBConnections;
  WConnection<T,Op> *BConnection;             // Connection BConnection[numBConnections];
  unsigned int InsertBConnection(unsigned int &j, WConnection<T,Op> &c);
//#ifdef PATH_COUNTING_ENABLED
  void InstallPathCounts();
  void InstallWeightsOfPathsAsAmpsToExits();
//#endif

 private:
  WeightedCFLOBDDInternalNode();                                         // Default constructor (hidden)
  WeightedCFLOBDDInternalNode(const WeightedCFLOBDDInternalNode &n);             // Copy constructor (hidden)
  WeightedCFLOBDDInternalNode& operator= (const WeightedCFLOBDDInternalNode &n); // Overloaded = (hidden)
};

//********************************************************************

// std::ostream& operator<< (std::ostream & out, const CFLOBDDTopNode &d);

//********************************************************************
// CFLOBDDLeafNode
//********************************************************************
template <typename T, typename Op>
class WeightedCFLOBDDLeafNode : public WeightedCFLOBDDNode<T,Op> {
 public:
  WeightedCFLOBDDLeafNode(bool flag = true);                   // Constructor
  WeightedCFLOBDDLeafNode(T lweight, T rweight);                   // Constructor
  virtual ~WeightedCFLOBDDLeafNode();          // Destructor
  virtual WCFLOBDD_NODEKIND NodeKind() const = 0;
  virtual void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index) = 0;
  virtual int Traverse(SH_OBDD::AssignmentIterator &ai) = 0;
  virtual std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valList, bool forceReduce = false) = 0;
  virtual unsigned int Hash(unsigned int modsize) = 0;
  void DumpConnections(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visited, std::ostream & out = std::cout);
  void CountNodesAndEdges(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
	  unsigned int &nodeCount, unsigned int &edgeCount, unsigned int& returnEdgesCount);
  void CountNodes(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, unsigned int &nodeCount);
  void CountPaths(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes);
  void ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes);
  virtual bool operator!= (const WeightedCFLOBDDNode<T,Op> & n) = 0;  // Overloaded !=
  virtual bool operator== (const WeightedCFLOBDDNode<T,Op> & n) = 0;  // Overloaded ==
  virtual void InstallWeightsOfPathsAsAmpsToExits() = 0;
  void IncrRef();
  void DecrRef();
  T lweight;
  T rweight;

 public:
	 virtual std::ostream& print(std::ostream & out = std::cout) const = 0;
	//  virtual bool IsValid() = 0;
};

//********************************************************************
// CFLOBDDForkNode
//********************************************************************

template <typename T, typename Op>
class WeightedCFLOBDDForkNode : public WeightedCFLOBDDLeafNode<T,Op> {
 public:
  WeightedCFLOBDDForkNode(bool flag = true);                   // Constructor
  WeightedCFLOBDDForkNode(T lweight, T rweight);                   // Constructor
  ~WeightedCFLOBDDForkNode();                  // Destructor
  WCFLOBDD_NODEKIND NodeKind() const { return W_CFLOBDD_FORK; }
  void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index);
  int Traverse(SH_OBDD::AssignmentIterator &ai);
  std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valList, bool forceReduce = false);
  unsigned int Hash(unsigned int modsize);
  bool operator!= (const WeightedCFLOBDDNode<T,Op> & n);        // Overloaded !=
  bool operator== (const WeightedCFLOBDDNode<T,Op> & n);        // Overloaded ==
  void InstallWeightsOfPathsAsAmpsToExits(); 

 public:
	 std::ostream& print(std::ostream & out = std::cout) const;
	//  bool IsValid();

 private:
  WeightedCFLOBDDForkNode(const WeightedCFLOBDDForkNode &n);   // Copy constructor (hidden)
  WeightedCFLOBDDForkNode& operator= (const WeightedCFLOBDDForkNode &n); // Overloaded = (hidden)
};

//********************************************************************
// CFLOBDDDontCareNode
//********************************************************************

template <typename T, typename Op>
class WeightedCFLOBDDDontCareNode : public WeightedCFLOBDDLeafNode<T,Op> {
 public:
  WeightedCFLOBDDDontCareNode(bool flag = true);                   // Constructor
  WeightedCFLOBDDDontCareNode(T lweight, T rweight);                   // Constructor
  ~WeightedCFLOBDDDontCareNode();                  // Destructor
  WCFLOBDD_NODEKIND NodeKind() const { return W_CFLOBDD_DONTCARE; }
  void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index);
  int Traverse(SH_OBDD::AssignmentIterator &ai);
  std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valList, bool forceReduce = false);
  unsigned int Hash(unsigned int modsize);
  bool operator!= (const WeightedCFLOBDDNode<T,Op> & n);        // Overloaded !=
  bool operator== (const WeightedCFLOBDDNode<T,Op> & n);        // Overloaded ==
  void InstallWeightsOfPathsAsAmpsToExits();

 public:
	 std::ostream& print(std::ostream & out = std::cout) const;
	//  bool IsValid();

 private:
  WeightedCFLOBDDDontCareNode(const WeightedCFLOBDDDontCareNode &n);   // Copy constructor (hidden)
  WeightedCFLOBDDDontCareNode& operator= (const WeightedCFLOBDDDontCareNode &n); // Overloaded = (hidden)
};

template <typename T, typename Op>
class WeightedBDDTopNode : public WeightedCFLOBDDNode<T,Op> {
  public:
    WeightedBDDTopNode(unsigned int numberOfVars);
    ~WeightedBDDTopNode();
    WCFLOBDD_NODEKIND NodeKind() const { return W_BDD_TOPNODE; }
    bool operator!= (const WeightedCFLOBDDNode<T,Op> &n);
    bool operator== (const WeightedCFLOBDDNode<T,Op> &n);
    unsigned int Hash(unsigned int modsize);
    void IncrRef();
    void DecrRef();
    WeightedBDDNodeHandle<T, Op> bddContents;
    unsigned int numberOfVars;

    void FillSatisfyingAssignment(unsigned int i, SH_OBDD::Assignment &assignment, unsigned int &index);
    int Traverse(SH_OBDD::AssignmentIterator &ai);
    std::pair<WeightedCFLOBDDNodeHandleT<T,Op>,T> Reduce(ReductionMapHandle& redMapHandle, unsigned int replacementNumExits, WeightedValuesListHandle<T>& valList, bool forceReduce = false);
    void DumpConnections(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visited, std::ostream & out = std::cout);
    void CountNodesAndEdges(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
	  unsigned int &nodeCount, unsigned int &edgeCount, unsigned int& returnEdgesCount);
    void CountNodes(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes, unsigned int &nodeCount);
    void CountPaths(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes);
    void ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes);
    
    public:
      std::ostream& print(std::ostream &out = std::cout) const;
};

 } // namespace CFL_OBDD

#endif
