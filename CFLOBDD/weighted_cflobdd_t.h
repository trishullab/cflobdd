#ifndef WEIGHTED_CFLOBDD_T_GUARD
#define WEIGHTED_CFLOBDD_T_GUARD


#include <iostream>
#include <fstream>
#include "ref_ptr.h"
#include "weighted_cflobdd_node_t.h"
#include "assignment.h"
#include "bool_op.h"
#include "weighted_cflobdd_top_node_t.h"

// TODO: Ideally, merge this class and CFLOBDD_T by making one a subclass of another.

namespace CFL_OBDD {

//********************************************************************
// WEIGHTED_CFLOBDD_T (use of Op as the composition operator)
//********************************************************************

// Representation of a Boolean function
template <typename T, typename Op>
class WEIGHTED_CFLOBDD_T {
  friend void GroupCountNodesAndEdgesStart(unsigned int &nodeCount, unsigned int &edgeCount);
  friend void GroupCountNodesAndEdgesEnd();
  friend void GroupDumpConnectionsStart();
  friend void GroupDumpConnectionsEnd();
 public:
  WEIGHTED_CFLOBDD_T();      // Default constructor (rep. of \a.true)
  WEIGHTED_CFLOBDD_T(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n); // Constructor
  WEIGHTED_CFLOBDD_T(const WEIGHTED_CFLOBDD_T &d);                            // Copy constructor
  ~WEIGHTED_CFLOBDD_T();                                             // Destructor
  static unsigned int const maxLevel;
  bool Evaluate(SH_OBDD::Assignment &assignment);        // Evaluate a Boolean function
//   void PrintYield(std::ostream * out);               // print the yield of the "tree"
//   bool IsValid();									// check if the CFLOBDD is valid (satisying all structural invariants)
#ifdef PATH_COUNTING_ENABLED
  unsigned int NumSatisfyingAssignments();      // Return number of satisfying assignments
 // mpz_class NumSatisfyingAssignments();      // Return number of satisfying assignments
#endif
  bool FindOneSatisfyingAssignment(SH_OBDD::Assignment * &assignment);
                                                // Find a satisfying assignment
  unsigned int Hash(unsigned int modsize);
  bool operator!= (const WEIGHTED_CFLOBDD_T & C) const;          // Overloaded !=
  bool operator== (const WEIGHTED_CFLOBDD_T & C) const;          // Overloaded ==
  WEIGHTED_CFLOBDD_T& operator= (const WEIGHTED_CFLOBDD_T &c);       // assignment
  ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>> root;

  void DumpConnections(std::ostream & out = std::cout);
  void GroupDumpConnections(std::ostream & out = std::cout);
  void CountNodesAndEdges(unsigned int &nodeCount, unsigned int &edgeCount, unsigned int &returnEdgesCount, unsigned int &returnEdgesObjCount);
  void CountNodes(unsigned int &nodeCount);
  void GroupCountNodesAndEdges(unsigned int &nodeCount, unsigned int &edgeCount);
  void CountPaths();
  void ComputeWeightOfPathsAsAmpsToExits();

 public:
	 std::ostream& print(std::ostream & out = std::cout) const;
  static Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodesDuringGroupCountNodesAndEdges;
  static Hashset<CFLOBDDReturnMapBody> *visitedEdgesDuringGroupCountNodesAndEdges;
  static Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedDuringGroupDumpConnections;
};

template <typename T, typename Op>
std::ostream& operator<< (std::ostream & out, const WEIGHTED_CFLOBDD_T<T, Op> &d);

void GroupCountNodesAndEdgesStart(unsigned int &nodeCount, unsigned int &edgeCount);
void GroupCountNodesAndEdgesEnd();

template <typename T, typename Op>
unsigned int const WEIGHTED_CFLOBDD_T<T, Op>::maxLevel = CFLOBDDMaxlevel;

// Default constructor must be provided by each specialization XXX
// template<>
// CFLOBDD_T<XXX>::CFLOBDD_T()
// {
// 	root = default XXX value, e.g., for int, we use MkTrueTop();
// }

template <typename T, typename Op>
WEIGHTED_CFLOBDD_T<T, Op>::WEIGHTED_CFLOBDD_T(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n)
{
	root = n;
}

// Copy constructor
template <typename T, typename Op>
WEIGHTED_CFLOBDD_T<T, Op>::WEIGHTED_CFLOBDD_T(const WEIGHTED_CFLOBDD_T<T, Op> &c)
{
	root = c.root;
}

template <typename T, typename Op>
WEIGHTED_CFLOBDD_T<T, Op>::~WEIGHTED_CFLOBDD_T()
{
}

// Operations -----------------------------------------------------

// Evaluate
//    Return the value of the Boolean function under the given assignment
template <typename T, typename Op>
bool WEIGHTED_CFLOBDD_T<T, Op>::Evaluate(SH_OBDD::Assignment &assignment)
{
	// TODO: uncomment and change to correct output
	return true;// root->EvaluateIteratively(assignment);
}


// PrintYield

// print the yield of the CFLOBDD (i.e., the leaves of 0's and 1's
// in "left-to-right order").

// template <typename T, typename Op>
// void WEIGHTED_CFLOBDD_T<T, Op>::PrintYield(std::ostream * out)
// {
// 	root->PrintYield(out);
// }

// template <typename T, typename Op>
// bool WEIGHTED_CFLOBDD_T<T, Op>::IsValid()
// {
// 	return root->IsValid();
// }

// Satisfaction Operations ------------------------------------

#ifdef PATH_COUNTING_ENABLED
// NumSatisfyingAssignments
//
// Return the number of satisfying assignments
//
// Running time: Linear in the size of the CFLOBDD
//
template <typename T, typename Op>
unsigned int CFLOBDD_T<T, Op>::NumSatisfyingAssignments()
{
	return root->NumSatisfyingAssignments();
}
#endif

// FindOneSatisfyingAssignment
//
// If a satisfying assignment exists, allocate and place such an
//    assignment in variable "assignment" and return true.
// Otherwise return false.
//
// Running time: Linear in the number of variables
//
template <typename T, typename Op>
bool WEIGHTED_CFLOBDD_T<T, Op>::FindOneSatisfyingAssignment(SH_OBDD::Assignment * &assignment)
{
	return root->FindOneSatisfyingAssignment(assignment);
}

// Hash
template <typename T, typename Op>
unsigned int WEIGHTED_CFLOBDD_T<T, Op>::Hash(unsigned int modsize)
{
	return root->Hash(modsize);
}

// Count Nodes
template <typename T, typename Op>
void WEIGHTED_CFLOBDD_T<T, Op>::CountNodes(unsigned int &nodeCount)
{
	Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes = new Hashset<WeightedCFLOBDDNodeHandleT<T,Op>>;
	nodeCount = 0;
	root->CountNodes(visitedNodes, nodeCount);
	delete visitedNodes;
}

// Count Paths
template <typename T, typename Op>
void WEIGHTED_CFLOBDD_T<T, Op>::CountPaths()
{
	Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes = new Hashset<WeightedCFLOBDDNodeHandleT<T,Op>>;
	root->CountPaths(visitedNodes);
	delete visitedNodes;
}

// Compute Weights of Paths
template <typename T, typename Op>
void WEIGHTED_CFLOBDD_T<T, Op>::ComputeWeightOfPathsAsAmpsToExits()
{
	Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes = new Hashset<WeightedCFLOBDDNodeHandleT<T,Op>>;
	root->ComputeWeightOfPathsAsAmpsToExits(visitedNodes);
	delete visitedNodes;
}

template <typename T, typename Op>
void WEIGHTED_CFLOBDD_T<T, Op>::CountNodesAndEdges(unsigned int &nodeCount, unsigned int &edgeCount, unsigned int &returnEdgesCount, unsigned int &returnEdgesObjCount)
{
	Hashset<WeightedCFLOBDDNodeHandleT<T,Op>> *visitedNodes = new Hashset<WeightedCFLOBDDNodeHandleT<T,Op>>;
	Hashset<CFLOBDDReturnMapBody> *visitedEdges = new Hashset<CFLOBDDReturnMapBody>;
	nodeCount = 0;
	edgeCount = 0;
	returnEdgesCount = 0;
	root->CountNodesAndEdges(visitedNodes, visitedEdges, nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
	returnEdgesObjCount = visitedEdges->Size();
	delete visitedNodes;
	//delete visitedEdges;
}

// Overloaded !=
template <typename T, typename Op>
bool WEIGHTED_CFLOBDD_T<T, Op>::operator!= (const WEIGHTED_CFLOBDD_T<T, Op> & C) const
{
	return *root != *C.root;
}

// Overloaded ==
template <typename T, typename Op>
bool WEIGHTED_CFLOBDD_T<T, Op>::operator== (const WEIGHTED_CFLOBDD_T<T, Op> & C) const
{
	return *root == *C.root;
}

// Overloaded assignment
template <typename T, typename Op>
WEIGHTED_CFLOBDD_T<T, Op> & WEIGHTED_CFLOBDD_T<T, Op>::operator= (const WEIGHTED_CFLOBDD_T<T, Op> &c)
{
	if (this != &c)      // don't assign to self!
	{
		root = c.root;
	}
	return *this;
}


// print
template <typename T, typename Op>
std::ostream& WEIGHTED_CFLOBDD_T<T, Op>::print(std::ostream & out) const
{
	out << *root << std::endl;
	return out;
}

template <typename T, typename Op>
std::ostream& operator<< (std::ostream & out, const WEIGHTED_CFLOBDD_T<T, Op> &d)
{
	d.print(out);
	return(out);
}

// Linear operations -----------------------------------------------

// Pointwise addition: \f.\g.(f + g)
template <typename T, typename Op>
WEIGHTED_CFLOBDD_T<T, Op> operator+(WEIGHTED_CFLOBDD_T<T, Op> f, WEIGHTED_CFLOBDD_T<T, Op> g)
{
	return WEIGHTED_CFLOBDD_T<T, Op>(MkPlusTopNode<T, Op>(f.root, g.root));
}

// Pointwise addition: \f.\g.(f ^ g)
template <typename T, typename Op>
WEIGHTED_CFLOBDD_T<T, Op> operator^(WEIGHTED_CFLOBDD_T<T, Op> f, WEIGHTED_CFLOBDD_T<T, Op> g)
{
	return WEIGHTED_CFLOBDD_T<T, Op>(MkExorTopNode<T, Op>(f.root, g.root));
}

//// Left scalar-multiplication: \c:unsigned int.\g.(c * g)
//template <typename T, typename Op>
//WEIGHTED_CFLOBDD_T<T, Op> operator*(unsigned int c, WEIGHTED_CFLOBDD_T<T, Op> g)
//{
//	return WEIGHTED_CFLOBDD_T<T, Op>(MkLeftScalarTimesTopNode(c, g.root));
//}
//
//// Left scalar-multiplication: \c:int.\g.(c * g)
//template <typename T, typename Op>
//WEIGHTED_CFLOBDD_T<T, Op> operator*(int c, WEIGHTED_CFLOBDD_T<T, Op> g)
//{
//	return WEIGHTED_CFLOBDD_T<T, Op>(MkLeftScalarTimesTopNode(c, g.root));
//}
//
//// Left scalar-multiplication: \c:double.\g.(c * g)
//template <typename T, typename Op>
//WEIGHTED_CFLOBDD_T<T, Op> operator*(double c, WEIGHTED_CFLOBDD_T<T, Op> g)
//{
//	return WEIGHTED_CFLOBDD_T<T, Op>(MkLeftScalarTimesTopNode(c, g.root));
//}

// Left scalar-multiplication: \c:unsigned int.\g.(c * g)
template<typename T, typename T1, typename Op>
WEIGHTED_CFLOBDD_T<T, Op> operator*(T1 c, WEIGHTED_CFLOBDD_T<T, Op> g)
{
	return WEIGHTED_CFLOBDD_T<T, Op>(MkLeftScalarTimesTopNode<T, T1, Op>(c, g.root));
}

// Right scalar-multiplication: \f.\c:int.(f * c)
template <typename T, typename Op>
WEIGHTED_CFLOBDD_T<T, Op> operator*(WEIGHTED_CFLOBDD_T<T, Op> f, int c)
{
	return WEIGHTED_CFLOBDD_T<T, Op>(MkRightScalarTimesTopNode<T, Op>(f.root, c));
}

// Pointwise multiplication: \f.\g.(f * g) -------------------------------------------
template <typename T, typename Op>
WEIGHTED_CFLOBDD_T<T, Op> operator*(WEIGHTED_CFLOBDD_T<T, Op> f, WEIGHTED_CFLOBDD_T<T, Op> g)
{
	return WEIGHTED_CFLOBDD_T<T, Op>(MkTimesTopNode<T, Op>(f.root, g.root));
}


} // namespace CFL_OBDD

#endif
