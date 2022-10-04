#ifndef CFLOBDD_TOP_NODE_T_GUARD
#define CFLOBDD_TOP_NODE_T_GUARD

//
//    Copyright (c) 2017 Thomas W. Reps
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
#include <boost/unordered_map.hpp>
#include "traverse_state_cfl.h"
#include "assignment.h"
#include "bool_op.h"
#include "hash.h"
#include "hashset.h"
#include "ref_ptr.h"
#include "connectionT.h"
#include "cflobdd_node.h"
#include "matmult_map.h"

namespace CFL_OBDD {

//********************************************************************
// CFLOBDDTopNode
//********************************************************************

template <typename T>
class CFLOBDDTopNodeT {
 public:
  CFLOBDDTopNodeT(CFLOBDDNode *n, ReturnMapHandle<T> &mapHandle);                // Constructor
  CFLOBDDTopNodeT(CFLOBDDNodeHandle &nodeHandle, ReturnMapHandle<T> &mapHandle); // Constructor
  ~CFLOBDDTopNodeT();                                   // Destructor

  T Evaluate(SH_OBDD::Assignment &assignment);             // Evaluate a Boolean function (recursive)
  T EvaluateIteratively(SH_OBDD::Assignment &assignment);  // Evaluate a Boolean function (iterative)
  void PrintYield(std::ostream * out);
  void PrintYieldAux(std::ostream * out, List<ConsCell<TraverseState> *> &T1, ConsCell<TraverseState> *S);
  void PrintYieldSemantic(std::ostream & out);       // print the yield of the "tree"
  bool IsValid();				// check if the CFLOBDD is valid (satisfies structural invariants)
#ifdef PATH_COUNTING_ENABLED
  unsigned int NumSatisfyingAssignments();
  //mpz_class NumSatisfyingAssignments();
#endif
  bool FindOneSatisfyingAssignment(SH_OBDD::Assignment * &assignment);
  unsigned int Hash(unsigned int modsize);
  void DumpConnections(Hashset<CFLOBDDNodeHandle> *visited, std::ostream & out = std::cout);
  void CountNodesAndEdges(Hashset<CFLOBDDNodeHandle> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
	  unsigned int &nodeCount, unsigned int &edgeCount, unsigned int& returnEdgesCount, unsigned int &returnEdgesObjCount);
  void CountNodes(Hashset<CFLOBDDNodeHandle> *visitedNodes, unsigned int &nodeCount);
  void CountPaths(Hashset<CFLOBDDNodeHandle> *visitedNodes);
  bool operator!= (const CFLOBDDTopNodeT<T> & C);          // Overloaded !=
  bool operator== (const CFLOBDDTopNodeT<T> & C);          // Overloaded ==
  static unsigned int const maxLevel;
  unsigned int level;
  ConnectionT<ReturnMapHandle<T>> rootConnection;                           // A single Connection
  RefCounter count;
  void DeallocateMemory();


 private:
  static Hashset<CFLOBDDTopNodeT<T> > *computedCache;       // TEMPORARY: should be HashCache
  CFLOBDDTopNodeT();                                    // Default constructor (hidden)
  CFLOBDDTopNodeT(const CFLOBDDTopNodeT<T> &n);             // Copy constructor (hidden)
  CFLOBDDTopNodeT& operator= (const CFLOBDDTopNodeT<T> &n); // Overloaded = (hidden)
 public:
  std::ostream& print(std::ostream & out = std::cout) const;

  // Usage: CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
  typedef ref_ptr<CFLOBDDTopNodeT<T>> CFLOBDDTopNodeTRefPtr;
};

template <typename T>
std::ostream& operator<< (std::ostream & out, const CFLOBDDTopNodeT<T> &d);

template <typename T> 
typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    ApplyAndReduce(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n1,
		typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n2,
		BoolOp op
	);
template <typename T> 
  typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    ApplyAndReduce(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n1,
        typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n2,
        T(*func)(T, T)
	);
template <typename T> 
  typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    ApplyAndReduce(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n1,
        typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n2,
        typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n3,
        BoolOp3 op
	);
template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    MkPlusTopNode(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f,
                typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g
	);
template<typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    operator+(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f, 
	typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g);

template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    MkExorTopNode(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f,
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g
    );

template<typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    operator^(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f, typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g);

template <typename T, typename T1>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    MkLeftScalarTimesTopNode(T1 c, typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g);

template<typename T, typename T1>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    operator*(T1 c, typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g);

template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    MkRightScalarTimesTopNode(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f, int c);

template<typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    operator*(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f, int c);

template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    MkTimesTopNode(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f,
                typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g);

template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    operator*(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f, typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g);

} // namespace CFL_OBDD

#endif
