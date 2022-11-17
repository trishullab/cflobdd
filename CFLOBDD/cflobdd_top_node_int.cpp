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

#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdarg>
//#include <mpirxx.h>

#include "cflobdd_top_node_int.h"
#include "cflobdd_node.h"

using namespace CFL_OBDD;

// FindOneSatisfyingAssignment
//
// If a satisfying assignment exists, allocate and place such an
//    assignment in variable "assignment" and return true.
// Otherwise return false.
//
// Running time: Linear in the number of variables
//
template <>
bool CFLOBDDTopNodeT<int>::FindOneSatisfyingAssignment(SH_OBDD::Assignment * &assignment)
{
	for (unsigned int i = 0; i < rootConnection.entryPointHandle->handleContents->numExits; i++) {
		unsigned int k = rootConnection.returnMapHandle.Lookup(i);
		if (k == 1) {  // A satisfying assignment must exist
			unsigned int size = 1 << CFLOBDDTopNode::maxLevel;
			assignment = new SH_OBDD::Assignment(size);
			rootConnection.entryPointHandle->handleContents->FillSatisfyingAssignment(i, *assignment, size);
			return true;
		}
	}
	return false;
}

#ifdef PATH_COUNTING_ENABLED
// NumSatisfyingAssignments
//
// Return the number of satisfying assignments
//
// Running time: Linear in the size of the CFLOBDDTopNode
//
template <>
unsigned int CFLOBDDTopNodeT<int>::NumSatisfyingAssignments()
{
	unsigned int ans = 0;

	for (unsigned int i = 0; i < rootConnection.entryPointHandle->handleContents->numExits; i++) {
		unsigned int k = rootConnection.returnMapHandle.Lookup(i);
		//unsigned int k = 1;
		if (k == 1) {
			ans += (int)rootConnection.entryPointHandle->handleContents->numPathsToExit[i];
		}
	}
	return ans;
}
#endif

/*
template <>
void CFLOBDDTopNodeT<int>::CountNodesAndEdges(Hashset<CFLOBDDNodeHandle> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, unsigned int &nodeCount, unsigned int &edgeCount)
{
	rootConnection.entryPointHandle->handleContents->CountNodesAndEdges(visitedNodes, visitedEdges, nodeCount, edgeCount);
	if (visitedEdges->Lookup(rootConnection.returnMapHandle.mapContents) == NULL) {
		visitedEdges->Insert(rootConnection.returnMapHandle.mapContents);
		edgeCount += rootConnection.returnMapHandle.Size();
	}
}
*/


//********************************************************************
// CFLOBDDTopNode
//********************************************************************

namespace CFL_OBDD {

template class CFLOBDDTopNodeT<int>;

// CFLOBDDTopNode-creation operations --------------------------------------

// Create representation of \x.true
CFLOBDDTopNodeIntRefPtr MkTrueTop(int level)
{
  CFLOBDDTopNodeIntRefPtr v;
  CFLOBDDReturnMapHandle m;

  m.AddToEnd(1);  // Map the only exit of the body to 1 (i.e., T)
  m.Canonicalize();
  if (level == -1)
    v = new CFLOBDDTopNode(CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDTopNode::maxLevel], m);
  else
    v = new CFLOBDDTopNode(CFLOBDDNodeHandle::NoDistinctionNode[level], m); 
  return v;
}

// Create representation of \x.false
CFLOBDDTopNodeIntRefPtr MkFalseTop(int level)
{
  CFLOBDDTopNodeIntRefPtr v;
  CFLOBDDReturnMapHandle m;

  m.AddToEnd(0);  // Map the only exit of the body to 0 (i.e., F)
  m.Canonicalize();
  if (level == -1)
    v = new CFLOBDDTopNode(CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDTopNode::maxLevel], m);
  else
    v = new CFLOBDDTopNode(CFLOBDDNodeHandle::NoDistinctionNode[level], m);
  return v;
}

// Create representation of \x.x_i
CFLOBDDTopNodeIntRefPtr MkDistinction(unsigned int i, int level)
{
  CFLOBDDTopNodeIntRefPtr v;
  CFLOBDDNodeHandle tempHandle;
  CFLOBDDReturnMapHandle m;

  assert(i < (1 << CFLOBDDTopNode::maxLevel));   // i.e., i < 2**maxLevel
  if (level == -1)
	tempHandle = MkDistinction(CFLOBDDTopNode::maxLevel, i);
  else
	  tempHandle = MkDistinction(level, i);
  m.AddToEnd(0);
  m.AddToEnd(1);
  m.Canonicalize();
  v = new CFLOBDDTopNode(tempHandle, m);
  return v;
}

// Create representation of addition relation with interleaved variables
// { (xi yi zi _)* | vec{x} + vec{y} = vec{z} }
CFLOBDDTopNodeIntRefPtr MkAdditionInterleavedRecursiveTop()
{
  CFLOBDDTopNodeIntRefPtr v;
  CFLOBDDNodeHandle n;
  CFLOBDDReturnMapHandle m10;

  n = MkAdditionInterleavedRecursive(CFLOBDDTopNode::maxLevel, false);

  // Reduce n by mapping the "carry=0" and "carry=1" exits to accept
     CFLOBDDReturnMapHandle retMapHandle;
     m10.AddToEnd(1);
     m10.AddToEnd(0);
     m10.Canonicalize();
     ReductionMapHandle reductionMapHandle;
     reductionMapHandle.AddToEnd(0);
     reductionMapHandle.AddToEnd(1);
     reductionMapHandle.AddToEnd(0);
     //CFLOBDDNodeHandle::InitReduceCache();
     CFLOBDDNodeHandle reduced_n = n.Reduce(reductionMapHandle, m10.Size());
     //CFLOBDDNodeHandle::DisposeOfReduceCache();

  // Create and return CFLOBDDTopNode
     v = new CFLOBDDTopNode(reduced_n, m10);
     return(v);
}

// Create representation of addition relation with interleaved variables
// { (xi yi zi _)* | vec{x} + vec{y} = vec{z} }
CFLOBDDTopNodeIntRefPtr MkAdditionInterleavedTop()
{
	CFLOBDDTopNodeIntRefPtr v;
	CFLOBDDNodeHandle n;
	CFLOBDDReturnMapHandle m10;

	n = CFLOBDDNodeHandle::AdditionInterleavedNoCarryNode[CFLOBDDTopNode::maxLevel-2];

	// Reduce n by mapping the "carry=0" and "carry=1" exits to accept
	CFLOBDDReturnMapHandle retMapHandle;
	m10.AddToEnd(1);
	m10.AddToEnd(0);
	m10.Canonicalize();
	ReductionMapHandle reductionMapHandle;
	reductionMapHandle.AddToEnd(0);
	reductionMapHandle.AddToEnd(1);
	reductionMapHandle.AddToEnd(0);
	//CFLOBDDNodeHandle::InitReduceCache();
	CFLOBDDNodeHandle reduced_n = n.Reduce(reductionMapHandle, m10.Size());
	//CFLOBDDNodeHandle::DisposeOfReduceCache();

	// Create and return CFLOBDDTopNode
	v = new CFLOBDDTopNode(reduced_n, m10);
	return(v);
}

// Create representation of parity function
CFLOBDDTopNodeIntRefPtr MkParityTop()
{
  CFLOBDDTopNodeIntRefPtr v;
  CFLOBDDNodeHandle tempHandle;
  CFLOBDDReturnMapHandle m;

  tempHandle = MkParity(CFLOBDDTopNode::maxLevel);
  m.AddToEnd(0);
  m.AddToEnd(1);
  m.Canonicalize();
  v = new CFLOBDDTopNode(tempHandle, m);
  return v;
}

// Create representation of step function
CFLOBDDTopNodeIntRefPtr MkStepUpOneFourthTop()
{
  CFLOBDDTopNodeIntRefPtr v;
  CFLOBDDNodeHandle tempHandle;
  CFLOBDDReturnMapHandle m;

  assert(CFLOBDDTopNode::maxLevel >= 1);
  tempHandle = MkStepOneFourth(CFLOBDDTopNode::maxLevel);
  m.AddToEnd(0);
  m.AddToEnd(1);
  m.Canonicalize();
  v = new CFLOBDDTopNode(tempHandle, m);
  return v;
}

// Create representation of step function
CFLOBDDTopNodeIntRefPtr MkStepDownOneFourthTop()
{
  return MkNot(MkStepUpOneFourthTop());
}

#ifdef ARBITRARY_STEP_FUNCTIONS
// Create representation of step function
CFLOBDDTopNodeIntRefPtr MkStepUpTop(unsigned int i)
{
  CFLOBDDTopNodeIntRefPtr v;
  CFLOBDDNodeHandle tempHandle;
  CFLOBDDReturnMapHandle m;

  assert(CFLOBDDTopNode::maxLevel <= 5);
  tempHandle = MkStepNode(CFLOBDDTopNode::maxLevel, i, 0, (1 << (1 << CFLOBDDTopNode::maxLevel)) - i);
  if (i != 0) {
    m.AddToEnd(0);
  }
  if (i != (1 << (1 << CFLOBDDTopNode::maxLevel))) {
    m.AddToEnd(1);
  }
  m.Canonicalize();
  v = new CFLOBDDTopNode(tempHandle, m);
  return v;
}

// Create representation of step function
CFLOBDDTopNodeIntRefPtr MkStepDownTop(unsigned int i)
{
  return MkNot(MkStepUpTop(i));
}
#endif

CFLOBDDTopNodeIntRefPtr shiftAtoBTop(CFLOBDDTopNodeIntRefPtr f, const unsigned int levelAtWhichToShift)
{
	CFLOBDDNodeHandle tempHandle = shiftAtoB(*(f->rootConnection.entryPointHandle), levelAtWhichToShift);
	CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, f->rootConnection.returnMapHandle);

	return v;
}

CFLOBDDTopNodeIntRefPtr shiftBtoATop(CFLOBDDTopNodeIntRefPtr f, const unsigned int levelAtWhichToShift)
{
	CFLOBDDNodeHandle tempHandle = shiftBtoA(*(f->rootConnection.entryPointHandle), levelAtWhichToShift);
	CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, f->rootConnection.returnMapHandle);
	return v;
}

CFLOBDDTopNodeIntRefPtr shiftAtoBAtLevelOneTop(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, CFLOBDDTopNodeIntRefPtr f)
{
	CFLOBDDNodeHandle tempHandle = shiftAtoBAtLevelOne(visitedNodes, totalVisitCount, redundantVisitCount, *(f->rootConnection.entryPointHandle));
	CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, f->rootConnection.returnMapHandle);
	return v;
}

CFLOBDDTopNodeIntRefPtr shiftBtoAAtLevelOneTop(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, CFLOBDDTopNodeIntRefPtr f)
{
	CFLOBDDNodeHandle tempHandle = shiftBtoAAtLevelOne(visitedNodes, totalVisitCount, redundantVisitCount, *(f->rootConnection.entryPointHandle));
	CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, f->rootConnection.returnMapHandle);
	return v;
}

CFLOBDDTopNodeIntRefPtr duplicateAinBAtLevelOneTop(CFLOBDDTopNodeIntRefPtr f)
{
	CFLOBDDNodeHandle tempHandle = duplicateAinBAtLevelOne(*(f->rootConnection.entryPointHandle));
	CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, f->rootConnection.returnMapHandle);
	return v;
}


// Unary operations on CFLOBDDTopNodes --------------------------------------

// Implements \f.!f
CFLOBDDTopNodeIntRefPtr MkNot(CFLOBDDTopNodeIntRefPtr f)
{
  CFLOBDDTopNodeIntRefPtr answer;
  CFLOBDDReturnMapHandle m = f->rootConnection.returnMapHandle.Complement();
  answer = new CFLOBDDTopNode(*(f->rootConnection.entryPointHandle), m);
  return answer;
}

// Binary operations on CFLOBDDTopNodes ------------------------------------

// \f.\g.(f && g)
CFLOBDDTopNodeIntRefPtr MkAnd(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, andOp);
}

// \f.\g.!(f && g)
CFLOBDDTopNodeIntRefPtr MkNand(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, nandOp);
}

// \f.\g.(f || g)
CFLOBDDTopNodeIntRefPtr MkOr(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, orOp);
}

// \f.\g.!(f || g)
CFLOBDDTopNodeIntRefPtr MkNor(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, norOp);
}

// \f.\g.(f == g)
CFLOBDDTopNodeIntRefPtr MkIff(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, iffOp);
}

// \f.\g.(f != g)
CFLOBDDTopNodeIntRefPtr MkExclusiveOr(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, exclusiveOrOp);
}

// \f.\g.(!f || g)
CFLOBDDTopNodeIntRefPtr MkImplies(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, impliesOp);
}

// \f.\g.(f && !g)
CFLOBDDTopNodeIntRefPtr MkMinus(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, minusOp);
}

// \f.\g.(!g || f)
CFLOBDDTopNodeIntRefPtr MkQuotient(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, quotientOp);
}

// \f.\g.(g && !f)
CFLOBDDTopNodeIntRefPtr MkNotQuotient(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, notQuotientOp);
}

// \f.\g.f
CFLOBDDTopNodeIntRefPtr MkFirst(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr )
{
  return f;
}

// \f.\g.!f
CFLOBDDTopNodeIntRefPtr MkNotFirst(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr )
{
  return MkNot(f);
}

// \f.\g.g
CFLOBDDTopNodeIntRefPtr MkSecond(CFLOBDDTopNodeIntRefPtr , CFLOBDDTopNodeIntRefPtr g)
{
  return g;
}

// \f.\g.!g
CFLOBDDTopNodeIntRefPtr MkNotSecond(CFLOBDDTopNodeIntRefPtr , CFLOBDDTopNodeIntRefPtr g)
{
  return MkNot(g);
}

/*
// \f.\g.(f + g)
CFLOBDDTopNodeIntRefPtr MkPlus(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, PlusFunc);
}

// \f.\g.(f * g)
CFLOBDDTopNodeIntRefPtr MkTimes(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g)
{
	return ApplyAndReduce<int>(f, g, TimesFunc);
}
*/

// Ternary operations on CFLOBDDTopNodes ------------------------------------

// \a.\b.\c.(a && b) || (!a && c)
CFLOBDDTopNodeIntRefPtr MkIfThenElse(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g, CFLOBDDTopNodeIntRefPtr h)
{
	return ApplyAndReduce<int>(f, g, h, ifThenElseOp);
}

// \a.\b.\c.(b && !a) || (c && !a) || (b && c)
CFLOBDDTopNodeIntRefPtr MkNegMajority(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g, CFLOBDDTopNodeIntRefPtr h)
{
	return ApplyAndReduce<int>(f, g, h, negMajorityOp);
}

CFLOBDDTopNodeIntRefPtr MkRestrict(CFLOBDDTopNodeIntRefPtr n, unsigned int i, bool val)
{
	CFLOBDDReturnMapHandle MapHandle;
	CFLOBDDNodeHandle g = Restrict(*(n->rootConnection.entryPointHandle), i, val,
		MapHandle);

	// Create returnMapHandle from MapHandle
	CFLOBDDReturnMapHandle returnMapHandle;
	unsigned MapSize = MapHandle.mapContents->mapArray.size();
	for (unsigned sBI = 0; sBI < MapSize; sBI++)
	{
		int d = MapHandle.mapContents->mapArray[sBI];
		int c = n->rootConnection.returnMapHandle.Lookup(d);
		returnMapHandle.AddToEnd(c);
	}
	returnMapHandle.Canonicalize();

	// Create and return CFLOBDDTopNode
	return(new CFLOBDDTopNode(g, returnMapHandle));
}

// Create representation of \f . exists x_i : f
CFLOBDDTopNodeIntRefPtr MkExists(CFLOBDDTopNodeIntRefPtr f, unsigned int i)
{
  CFLOBDDTopNodeIntRefPtr tempTrue = MkRestrict(f, i, true);
  CFLOBDDTopNodeIntRefPtr tempFalse = MkRestrict(f, i, false);
  return MkOr(tempTrue, tempFalse);
}

// Create representation of \f . forall x_i : f
CFLOBDDTopNodeIntRefPtr MkForall(CFLOBDDTopNodeIntRefPtr f, unsigned int i)
{
  CFLOBDDTopNodeIntRefPtr tempTrue = MkRestrict(f, i, true);
  CFLOBDDTopNodeIntRefPtr tempFalse = MkRestrict(f, i, false);
  return MkAnd(tempTrue, tempFalse);
}

CFLOBDDTopNodeIntRefPtr MkComposeTop(CFLOBDDTopNodeIntRefPtr f, int i, CFLOBDDTopNodeIntRefPtr g)              // \f. f | x_i = g
{
	  // DLC inefficient

  // Simple but slow method, see Bryant 1986 _GBAfBFM_:
  
  // f|x_i=g is g*(f|x_i=1) + (!g)*(f|x_i=0)
  // Better would be MkITE(g, f__x_i_1, f__x_i_0);

  CFLOBDDTopNodeIntRefPtr f__x_i_1 = MkRestrict(f, i, true);
  CFLOBDDTopNodeIntRefPtr f__x_i_0 = MkRestrict(f, i, false);
  CFLOBDDTopNodeIntRefPtr not_g = MkNot(g);

  CFLOBDDTopNodeIntRefPtr left_product = MkAnd(g, f__x_i_1);

  CFLOBDDTopNodeIntRefPtr right_product = MkAnd(not_g, f__x_i_0);

  CFLOBDDTopNodeIntRefPtr composition = MkOr(left_product, right_product);

  return composition;
}

double ComputeProbabilityTop(CFLOBDDTopNodeIntRefPtr f, std::vector<double>& var_probs)
{
  // assumption that the output contains only 0 and 1 values.
  std::vector<double> path_probs (f->rootConnection.returnMapHandle.Size(), 1);
  int index = f->rootConnection.returnMapHandle.LookupInv(0);
  if (index == -1)
    return 1;
  if (index == 0 && f->rootConnection.returnMapHandle.Size() == 1)
    return 0;
  path_probs[index] = 0;
  return ComputeProbabilityNode(*(f->rootConnection.entryPointHandle), var_probs, path_probs, 0, pow(2, f->level)-1);
}

std::vector<double> ComputeProbabilityOfListTop(CFLOBDDTopNodeIntRefPtr f, std::vector<std::vector<double>>& var_probs)
{
  // assumption that the output contains only 0 and 1 values.
  std::vector<std::vector<double>> path_probs; 
 
  for (int i = 0; i < f->rootConnection.returnMapHandle.Size(); i++){
     std::vector<double> tmp(var_probs[0].size(), 1);
     path_probs.push_back(tmp);
  }
  int index = f->rootConnection.returnMapHandle.LookupInv(0);
  if (index == -1){
    return std::vector<double>(var_probs[0].size(), 1);
  }
  if (index == 0 && f->rootConnection.returnMapHandle.Size() == 1)
    return std::vector<double>(var_probs[0].size(), 0);
  for (int i = 0; i < var_probs[0].size(); i++)
    path_probs[index][i] = 0;
  return ComputeProbabilityOfListNode(*(f->rootConnection.entryPointHandle), var_probs, path_probs, 0, pow(2, f->level)-1);
}

// std::vector<double> ComputeEntropyOfListTop(CFLOBDDTopNodeIntRefPtr f, std::vector<std::vector<double>>& var_probs)
// {
//   // assumption that the output contains only 0 and 1 values.
//   std::vector<std::vector<double>> path_probs; 
//   std::vector<std::vector<double>> entropy;
//   for (int i = 0; i < f->rootConnection.returnMapHandle.Size(); i++){
//     std::vector<double> tmp(var_probs[0].size(), 1);
//     path_probs.push_back(tmp);
//     entropy.push_back(tmp);
//   }
//   int index = f->rootConnection.returnMapHandle.LookupInv(0);
//   if (index == -1){
//     return std::vector<double>(var_probs[0].size(), 1);
//   }
//   // if (index == 0 && f->rootConnection.returnMapHandle.Size() == 1)
//   //   return std::vector<double>(var_probs[0].size(), 0);
//   for (int i = 0; i < var_probs[0].size(); i++){
//     path_probs[index][i] = 0;
//     entropy[index][i] = 0;
//   }
//   return ComputeEntropyOfListNode(*(f->rootConnection.entryPointHandle), var_probs, path_probs, entropy, 0, pow(2, f->level)-1);
// }


}
