#ifndef CFLOBDD_TOP_NODE_INT_GUARD
#define CFLOBDD_TOP_NODE_INT_GUARD

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
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.s
//


// #include <iostream>
// #include <fstream>
#include "cflobdd_top_node_t.h"

namespace CFL_OBDD {

typedef CFLOBDDTopNodeT<int> CFLOBDDTopNode;
typedef CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeIntRefPtr;

// CFLOBDDTopNode-creation operations --------------------------------------
extern CFLOBDDTopNodeIntRefPtr MkTrueTop(int level=-1);                    // Representation of \x.true
extern CFLOBDDTopNodeIntRefPtr MkFalseTop(int level=-1);                   // Representation of \x.false
extern CFLOBDDTopNodeIntRefPtr MkDistinction(unsigned int i, int level = -1);  // Representation of \x.x_i
extern CFLOBDDTopNodeIntRefPtr MkAdditionInterleavedRecursiveTop();     // Representation of addition relation, created recursively
extern CFLOBDDTopNodeIntRefPtr MkAdditionInterleavedTop();     // Representation of addition relation
extern CFLOBDDTopNodeIntRefPtr MkParityTop();                  // Representation of parity function
extern CFLOBDDTopNodeIntRefPtr MkStepUpOneFourthTop();         // Representation of step function
extern CFLOBDDTopNodeIntRefPtr MkStepDownOneFourthTop();       // Representation of step function
#ifdef ARBITRARY_STEP_FUNCTIONS
  extern CFLOBDDTopNodeIntRefPtr MkStepUpTop(unsigned int i);    // Representation of step function
  extern CFLOBDDTopNodeIntRefPtr MkStepDownTop(unsigned int i);  // Representation of step function
#endif

extern CFLOBDDTopNodeIntRefPtr shiftAtoBTop(CFLOBDDTopNodeIntRefPtr f, const unsigned int levelAtWhichToShift);
extern CFLOBDDTopNodeIntRefPtr shiftBtoATop(CFLOBDDTopNodeIntRefPtr f, const unsigned int levelAtWhichToShift);

extern CFLOBDDTopNodeIntRefPtr shiftAtoBAtLevelOneTop(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, CFLOBDDTopNodeIntRefPtr f);
extern CFLOBDDTopNodeIntRefPtr shiftBtoAAtLevelOneTop(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, CFLOBDDTopNodeIntRefPtr f);
extern CFLOBDDTopNodeIntRefPtr duplicateAinBAtLevelOneTop(CFLOBDDTopNodeIntRefPtr f);

// Unary operations on CFLOBDDTopNodes --------------------------------------
extern CFLOBDDTopNodeIntRefPtr MkNot(CFLOBDDTopNodeIntRefPtr f);               // \f.!f

// Binary operations on CFLOBDDTopNodes --------------------------------------
extern CFLOBDDTopNodeIntRefPtr MkAnd(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);          // \f.\g.(f && g)
extern CFLOBDDTopNodeIntRefPtr MkNand(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);         // \f.\g.!(f && g)
extern CFLOBDDTopNodeIntRefPtr MkOr(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);           // \f.\g.(f || g)
extern CFLOBDDTopNodeIntRefPtr MkNor(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);          // \f.\g.!(f || g)
extern CFLOBDDTopNodeIntRefPtr MkIff(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);          // \f.\g.(f == g)
extern CFLOBDDTopNodeIntRefPtr MkExclusiveOr(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);  // \f.\g.(f != g)
extern CFLOBDDTopNodeIntRefPtr MkImplies(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);      // \f.\g.(!f || g)
extern CFLOBDDTopNodeIntRefPtr MkMinus(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);        // \f.\g.(f && !g)
extern CFLOBDDTopNodeIntRefPtr MkQuotient(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);     // \f.\g.(!g || f)
extern CFLOBDDTopNodeIntRefPtr MkNotQuotient(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);  // \f.\g.(g && !f)
extern CFLOBDDTopNodeIntRefPtr MkFirst(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);        // \f.\g.f
extern CFLOBDDTopNodeIntRefPtr MkNotFirst(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);     // \f.\g.!f
extern CFLOBDDTopNodeIntRefPtr MkSecond(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);       // \f.\g.g
extern CFLOBDDTopNodeIntRefPtr MkNotSecond(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);    // \f.\g.!g

// extern CFLOBDDTopNodeIntRefPtr MkPlus(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);       // \f.\g.(f + g)
// extern CFLOBDDTopNodeIntRefPtr MkTimes(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g);       // \f.\g.(f * g)

// Ternary operations on CFLOBDDTopNodes ------------------------------------

// \a.\b.\c.(a && b) || (!a && c)
extern CFLOBDDTopNodeIntRefPtr MkIfThenElse(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g, CFLOBDDTopNodeIntRefPtr h);

extern CFLOBDDTopNodeIntRefPtr MkRestrict(CFLOBDDTopNodeIntRefPtr n, unsigned int i, bool val);

// \a.\b.\c.(b && !a) || (c && !a) || (b && c)
extern CFLOBDDTopNodeIntRefPtr MkNegMajority(CFLOBDDTopNodeIntRefPtr f, CFLOBDDTopNodeIntRefPtr g, CFLOBDDTopNodeIntRefPtr h);

extern CFLOBDDTopNodeIntRefPtr MkExists(CFLOBDDTopNodeIntRefPtr f, unsigned int i);              // \f. exists x_i : f
extern CFLOBDDTopNodeIntRefPtr MkForall(CFLOBDDTopNodeIntRefPtr f, unsigned int i);              // \f. forall x_i : f
extern CFLOBDDTopNodeIntRefPtr MkComposeTop(CFLOBDDTopNodeIntRefPtr f, int i, CFLOBDDTopNodeIntRefPtr g);              // \f. f | x_i = g

extern double ComputeProbabilityTop(CFLOBDDTopNodeIntRefPtr f, std::vector<double>& probs);
extern std::vector<double> ComputeProbabilityOfListTop(CFLOBDDTopNodeIntRefPtr f, std::vector<std::vector<double>>& probs);
// extern std::vector<double> ComputeEntropyOfListTop(CFLOBDDTopNodeIntRefPtr f, std::vector<std::vector<double>>& probs);

} // namespace CFL_OBDD

#endif
