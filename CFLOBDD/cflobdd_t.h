#ifndef CFLOBDD_T_GUARD
#define CFLOBDD_T_GUARD

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

#include <iostream>
#include <fstream>
#include "ref_ptr.h"
#include "cflobdd_node.h"
#include "assignment.h"
#include "bool_op.h"
#include "cflobdd_top_node_t.h"

namespace CFL_OBDD {

//********************************************************************
// CFLOBDD_T
//********************************************************************

// Representation of a Boolean function
template <typename T>
class CFLOBDD_T {
  friend void GroupCountNodesAndEdgesStart(unsigned int &nodeCount, unsigned int &edgeCount);
  friend void GroupCountNodesAndEdgesEnd();
  friend void GroupDumpConnectionsStart();
  friend void GroupDumpConnectionsEnd();
 public:
  CFLOBDD_T();      // Default constructor (rep. of \a.true)
  CFLOBDD_T(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n); // Constructor
  CFLOBDD_T(const CFLOBDD_T &d);                            // Copy constructor
  ~CFLOBDD_T();                                             // Destructor
  static unsigned int const maxLevel;
  bool Evaluate(SH_OBDD::Assignment &assignment);        // Evaluate a Boolean function
  void PrintYield(std::ostream * out);               // print the yield of the "tree"
  bool IsValid();									// check if the CFLOBDD is valid (satisying all structural invariants)
#ifdef PATH_COUNTING_ENABLED
  unsigned int NumSatisfyingAssignments();      // Return number of satisfying assignments
 // mpz_class NumSatisfyingAssignments();      // Return number of satisfying assignments
#endif
  bool FindOneSatisfyingAssignment(SH_OBDD::Assignment * &assignment);
                                                // Find a satisfying assignment
  unsigned int Hash(unsigned int modsize);
  bool operator!= (const CFLOBDD_T & C) const;          // Overloaded !=
  bool operator== (const CFLOBDD_T & C) const;          // Overloaded ==
  CFLOBDD_T& operator= (const CFLOBDD_T &c);       // assignment
  ref_ptr<CFLOBDDTopNodeT<T>> root;

  void DumpConnections(std::ostream & out = std::cout);
  void GroupDumpConnections(std::ostream & out = std::cout);
  void CountNodesAndEdges(unsigned int &nodeCount, unsigned int &edgeCount, unsigned int &returnEdgesCount, unsigned int &returnEdgesObjCount);
  void CountNodes(unsigned int &nodeCount);
  void GroupCountNodesAndEdges(unsigned int &nodeCount, unsigned int &edgeCount);
  void CountPaths();

 public:
	 std::ostream& print(std::ostream & out = std::cout) const;
  static Hashset<CFLOBDDNodeHandle> *visitedNodesDuringGroupCountNodesAndEdges;
  static Hashset<CFLOBDDReturnMapBody> *visitedEdgesDuringGroupCountNodesAndEdges;
  static Hashset<CFLOBDDNodeHandle> *visitedDuringGroupDumpConnections;
};

template <typename T>
std::ostream& operator<< (std::ostream & out, const CFLOBDD_T<T> &d);

void GroupCountNodesAndEdgesStart(unsigned int &nodeCount, unsigned int &edgeCount);
void GroupCountNodesAndEdgesEnd();

template <typename T>
unsigned int const CFLOBDD_T<T>::maxLevel = CFLOBDDMaxLevel;

// Default constructor must be provided by each specialization XXX
// template<>
// CFLOBDD_T<XXX>::CFLOBDD_T()
// {
// 	root = default XXX value, e.g., for int, we use MkTrueTop();
// }

template<typename T>
CFLOBDD_T<T>::CFLOBDD_T(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n)
{
	root = n;
}

// Copy constructor
template<typename T>
CFLOBDD_T<T>::CFLOBDD_T(const CFLOBDD_T<T> &c)
{
	root = c.root;
}

template<typename T>
CFLOBDD_T<T>::~CFLOBDD_T()
{
}

// Operations -----------------------------------------------------

// Evaluate
//    Return the value of the Boolean function under the given assignment
template<typename T>
bool CFLOBDD_T<T>::Evaluate(SH_OBDD::Assignment &assignment)
{
	// TODO: uncomment and change to correct output
	return true;// root->EvaluateIteratively(assignment);
}


// PrintYield
//
// print the yield of the CFLOBDD (i.e., the leaves of 0's and 1's
// in "left-to-right order").
//
template<typename T>
void CFLOBDD_T<T>::PrintYield(std::ostream * out)
{
	root->PrintYield(out);
}

template<typename T>
bool CFLOBDD_T<T>::IsValid()
{
	return root->IsValid();
}

// Satisfaction Operations ------------------------------------

#ifdef PATH_COUNTING_ENABLED
// NumSatisfyingAssignments
//
// Return the number of satisfying assignments
//
// Running time: Linear in the size of the CFLOBDD
//
template<typename T>
unsigned int CFLOBDD_T<T>::NumSatisfyingAssignments()
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
template<typename T>
bool CFLOBDD_T<T>::FindOneSatisfyingAssignment(SH_OBDD::Assignment * &assignment)
{
	return root->FindOneSatisfyingAssignment(assignment);
}

// Hash
template<typename T>
unsigned int CFLOBDD_T<T>::Hash(unsigned int modsize)
{
	return root->Hash(modsize);
}

// Count Nodes
template<typename T>
void CFLOBDD_T<T>::CountNodes(unsigned int &nodeCount)
{
	Hashset<CFLOBDDNodeHandle> *visitedNodes = new Hashset<CFLOBDDNodeHandle>;
	nodeCount = 0;
	root->CountNodes(visitedNodes, nodeCount);
	delete visitedNodes;
}

// Count Paths
template<typename T>
void CFLOBDD_T<T>::CountPaths()
{
	Hashset<CFLOBDDNodeHandle> *visitedNodes = new Hashset<CFLOBDDNodeHandle>;
	root->CountPaths(visitedNodes);
	delete visitedNodes;
}

template<typename T>
void CFLOBDD_T<T>::CountNodesAndEdges(unsigned int &nodeCount, unsigned int &edgeCount, unsigned int &returnEdgesCount, unsigned int &returnEdgesObjCount)
{
	Hashset<CFLOBDDNodeHandle> *visitedNodes = new Hashset<CFLOBDDNodeHandle>;
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
template<typename T>
bool CFLOBDD_T<T>::operator!= (const CFLOBDD_T<T> & C) const
{
	return *root != *C.root;
}

// Overloaded ==
template<typename T>
bool CFLOBDD_T<T>::operator== (const CFLOBDD_T<T> & C) const
{
	return *root == *C.root;
}

// Overloaded assignment
template<typename T>
CFLOBDD_T<T> & CFLOBDD_T<T>::operator= (const CFLOBDD_T<T> &c)
{
	if (this != &c)      // don't assign to self!
	{
		root = c.root;
	}
	return *this;
}


// print
template<typename T>
std::ostream& CFLOBDD_T<T>::print(std::ostream & out) const
{
	out << *root << std::endl;
	return out;
}

template<typename T>
std::ostream& operator<< (std::ostream & out, const CFLOBDD_T<T> &d)
{
	d.print(out);
	return(out);
}

// Linear operations -----------------------------------------------

// Pointwise addition: \f.\g.(f + g)
template<typename T>
CFLOBDD_T<T> operator+(CFLOBDD_T<T> f, CFLOBDD_T<T> g)
{
	return CFLOBDD_T<T>(MkPlusTopNode<T>(f.root, g.root));
}

// Pointwise addition: \f.\g.(f ^ g)
template<typename T>
CFLOBDD_T<T> operator^(CFLOBDD_T<T> f, CFLOBDD_T<T> g)
{
	return CFLOBDD_T<T>(MkExorTopNode<T>(f.root, g.root));
}

//// Left scalar-multiplication: \c:unsigned int.\g.(c * g)
//template<typename T>
//CFLOBDD_T<T> operator*(unsigned int c, CFLOBDD_T<T> g)
//{
//	return CFLOBDD_T<T>(MkLeftScalarTimesTopNode(c, g.root));
//}
//
//// Left scalar-multiplication: \c:int.\g.(c * g)
//template<typename T>
//CFLOBDD_T<T> operator*(int c, CFLOBDD_T<T> g)
//{
//	return CFLOBDD_T<T>(MkLeftScalarTimesTopNode(c, g.root));
//}
//
//// Left scalar-multiplication: \c:double.\g.(c * g)
//template<typename T>
//CFLOBDD_T<T> operator*(double c, CFLOBDD_T<T> g)
//{
//	return CFLOBDD_T<T>(MkLeftScalarTimesTopNode(c, g.root));
//}

// Left scalar-multiplication: \c:unsigned int.\g.(c * g)
template<typename T, typename T1>
CFLOBDD_T<T> operator*(T1 c, CFLOBDD_T<T> g)
{
	return CFLOBDD_T<T>(MkLeftScalarTimesTopNode<T>(c, g.root));
}

// Right scalar-multiplication: \f.\c:int.(f * c)
template<typename T>
CFLOBDD_T<T> operator*(CFLOBDD_T<T> f, int c)
{
	return CFLOBDD_T<T>(MkRightScalarTimesTopNode<T>(f.root, c));
}

// Pointwise multiplication: \f.\g.(f * g) -------------------------------------------
template<typename T>
CFLOBDD_T<T> operator*(CFLOBDD_T<T> f, CFLOBDD_T<T> g)
{
	return CFLOBDD_T<T>(MkTimesTopNode<T>(f.root, g.root));
}


} // namespace CFL_OBDD

#endif
