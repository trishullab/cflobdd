#ifndef WEIGHTED_CFLOBDD_TOP_NODE_T_GUARD
#define WEIGHTED_CFLOBDD_TOP_NODE_T_GUARD


#include <iostream>
#include <fstream>
#include <boost/unordered_map.hpp>
#include "assignment.h"
#include "bool_op.h"
#include "hash.h"
#include "hashset.h"
#include "ref_ptr.h"
#include "weighted_cflobdd_node_t.h"
#include "weighted_rootConnectionT.h"
// #include "traverse_state_cfl.h"

namespace CFL_OBDD {

//********************************************************************
// CFLOBDDTopNode
//********************************************************************

template <typename T, typename Op>
class WeightedCFLOBDDTopNodeT {
 public:
  WeightedCFLOBDDTopNodeT(WeightedCFLOBDDNode<T, Op> *n, ReturnMapHandle<T> &mapHandle, T factor = 1);                // Constructor
  WeightedCFLOBDDTopNodeT(WeightedCFLOBDDNodeHandleT<T, Op> &nodeHandle, ReturnMapHandle<T> &mapHandle, T factor = 1); // Constructor
  ~WeightedCFLOBDDTopNodeT();                                   // Destructor

  T Evaluate(SH_OBDD::Assignment &assignment);             // Evaluate a Boolean function (recursive)
//   T EvaluateIteratively(SH_OBDD::Assignment &assignment);  // Evaluate a Boolean function (iterative)
//   void PrintYield(std::ostream * out);
//   void PrintYieldAux(std::ostream * out, List<ConsCell<TraverseState<WeightedCFLOBDDNode<T,Op>>> *> &T1, ConsCell<TraverseState<WeightedCFLOBDDNode<T,Op>>> *S);
//   void PrintYieldSemantic(std::ostream & out);       // print the yield of the "tree"
//   bool IsValid();				// check if the CFLOBDD is valid (satisfies structural invariants)
#ifdef PATH_COUNTING_ENABLED
  unsigned int NumSatisfyingAssignments();
  //mpz_class NumSatisfyingAssignments();
#endif
  bool FindOneSatisfyingAssignment(SH_OBDD::Assignment * &assignment);
  unsigned int Hash(unsigned int modsize);
  void DumpConnections(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visited, std::ostream & out = std::cout);
  void CountNodesAndEdges(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
	  unsigned int &nodeCount, unsigned int &edgeCount, unsigned int& returnEdgesCount, unsigned int &returnEdgesObjCount);
  void CountNodes(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visitedNodes, unsigned int &nodeCount);
  void CountPaths(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visitedNodes);
  void ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visitedNodes);
  bool operator!= (const WeightedCFLOBDDTopNodeT<T, Op> & C);          // Overloaded !=
  bool operator== (const WeightedCFLOBDDTopNodeT<T, Op> & C);          // Overloaded ==
  static unsigned int const maxLevel;
  unsigned int level;
  WRootConnection<ReturnMapHandle<T>, T, Op> rootConnection;                           // A single Connection
  RefCounter count;
  void DeallocateMemory();


 private:
  static Hashset<WeightedCFLOBDDTopNodeT<T, Op> > *computedCache;       // TEMPORARY: should be HashCache
  WeightedCFLOBDDTopNodeT();                                    // Default constructor (hidden)
  WeightedCFLOBDDTopNodeT(const WeightedCFLOBDDTopNodeT<T, Op> &n);             // Copy constructor (hidden)
  WeightedCFLOBDDTopNodeT& operator= (const WeightedCFLOBDDTopNodeT<T, Op> &n); // Overloaded = (hidden)
 public:
  std::ostream& print(std::ostream & out = std::cout) const;

  // Usage: CFLOBDDTopNodeT<T, Op>::CFLOBDDTopNodeTRefPtr
  typedef ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>> WeightedCFLOBDDTopNodeTRefPtr;
};

template <typename T, typename Op>
std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDTopNodeT<T, Op> &d);

template <typename T, typename Op> 
typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    ApplyAndReduce(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n1,
		typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTTRefPtr n2,
		BoolOp op
	);
template <typename T, typename Op> 
  typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    ApplyAndReduce(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n1,
        typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n2,
        T(*func)(T, T),
        bool flag
	);
template <typename T, typename Op> 
  typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    ApplyAndReduce(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n1,
        typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n2,
        typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n3,
        BoolOp3 op
	);
template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    MkPlusTopNode(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f,
                typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g
	);
template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    operator+(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f, 
	typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g);

template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    MkExorTopNode(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f,
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g
    );

template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    operator^(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f, typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g);

template <typename T, typename T1, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    MkLeftScalarTimesTopNode(T1 c, typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g);

template<typename T, typename T1, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    operator*(T1 c, typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g);

template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    MkRightScalarTimesTopNode(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f, int c);

template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    operator*(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f, int c);

template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    MkTimesTopNode(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f,
                typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g);

template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    operator*(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f, typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g);

} // namespace CFL_OBDD

#endif
