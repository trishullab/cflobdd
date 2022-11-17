#ifndef CFLOBDD_INT_GUARD
#define CFLOBDD_INT_GUARD

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
#include "cflobdd_t.h"

namespace CFL_OBDD {

typedef CFLOBDD_T<int> CFLOBDD;

// CFLOBDD-creation operations --------------------------------------
extern CFLOBDD MkTrue(int level = -1);                             // Representation of \x.true
extern CFLOBDD MkFalse(int level = -1);                            // Representation of \x.false
extern CFLOBDD MkProjection(unsigned int i, int level = -1);         // Representation of \x.x_i
extern CFLOBDD MkAdditionInterleavedRecursive();     // Representation of addition relation { (xi yi zi _)* | vec{x} + vec{y} = vec{z} }, created recursively (expensive)
extern CFLOBDD MkAdditionInterleavedBruteForce();    // Representation of addition relation { (xi yi zi _)* | vec{x} + vec{y} = vec{z} }, created by brute force (expensive)
extern CFLOBDD MkAdditionInterleaved();              // Representation of addition relation { (xi yi zi _)* | vec{x} + vec{y} = vec{z} }, created by tabulation
extern CFLOBDD MkX(unsigned int i);                  // Bits of i -> vec{x} for addition relation
extern CFLOBDD MkY(unsigned int i);                  // Bits of i -> vec{y} for addition relation
extern CFLOBDD MkZ(unsigned int i);                  // Bits of i -> vec{z} for addition relation
extern CFLOBDD MkXBit(unsigned int j);
extern CFLOBDD MkYBit(unsigned int j);
extern CFLOBDD MkZBit(unsigned int j);
extern CFLOBDD MkMultiplicationInterleavedBruteForce(); // Representation of multiplication relation { (xi yi zi _)* | vec{x} * vec{y} = vec{z} }
extern bool FactorZ(CFLOBDD R, unsigned int z, CFLOBDD &f1, CFLOBDD &f2, unsigned int &v1, unsigned int &v2);   // Return true if z has a non-trivial factorization
extern CFLOBDD MkParity();                           // Representation of parity function
extern CFLOBDD MkStepUpOneFourth();                  // Representation of step function
extern CFLOBDD MkStepDownOneFourth();                // Representation of step function
#ifdef ARBITRARY_STEP_FUNCTIONS
  extern CFLOBDD MkStepUp(unsigned int i);           // Representation of step function
  extern CFLOBDD MkStepDown(unsigned int i);         // Representation of step function
  extern CFLOBDD MkPulseUp(unsigned int i, unsigned int j);
                                              // Representation of pulse function
  extern CFLOBDD MkPulseDown(unsigned int i, unsigned int j);
                                              // Representation of step function
#endif

  extern CFLOBDD shiftAtoB(const CFLOBDD f, const unsigned int levelAtWhichToShift = CFLOBDD_MAX_LEVEL);
											  // Special operation used for Karatsuba multiplication
											  // with ASCENDING_INTERLEAVED variable ordering
extern CFLOBDD shiftBtoA(const CFLOBDD f, const unsigned int levelAtWhichToShift = CFLOBDD_MAX_LEVEL);
											  // Special operation used for Karatsuba multiplication
											  // with DESCENDING_INTERLEAVED variable ordering

extern CFLOBDD shiftAtoBAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, const CFLOBDD f);
extern CFLOBDD shiftBtoAAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, const CFLOBDD f);
extern CFLOBDD duplicateAinBAtLevelOne(const CFLOBDD f);

// Unary operations on CFLOBDDs --------------------------------------
extern CFLOBDD MkNot(CFLOBDD f);                     // \f.!f

// Binary operations on CFLOBDDs --------------------------------------
extern CFLOBDD MkAnd(CFLOBDD f, CFLOBDD g);          // \f.\g.(f && g)
extern CFLOBDD MkNand(CFLOBDD f, CFLOBDD g);         // \f.\g.!(f && g)
extern CFLOBDD MkOr(CFLOBDD f, CFLOBDD g);           // \f.\g.(f || g)
extern CFLOBDD MkNor(CFLOBDD f, CFLOBDD g);          // \f.\g.!(f || g)
extern CFLOBDD MkIff(CFLOBDD f, CFLOBDD g);          // \f.\g.(f == g)
extern CFLOBDD MkExclusiveOr(CFLOBDD f, CFLOBDD g);  // \f.\g.(f != g)
extern CFLOBDD MkImplies(CFLOBDD f, CFLOBDD g);      // \f.\g.(!f || g)
extern CFLOBDD MkMinus(CFLOBDD f, CFLOBDD g);        // \f.\g.(f && !g)
extern CFLOBDD MkQuotient(CFLOBDD f, CFLOBDD g);     // \f.\g.(!g || f)
extern CFLOBDD MkNotQuotient(CFLOBDD f, CFLOBDD g);  // \f.\g.(g && !f)
extern CFLOBDD MkFirst(CFLOBDD f, CFLOBDD g);        // \f.\g.f
extern CFLOBDD MkNotFirst(CFLOBDD f, CFLOBDD g);     // \f.\g.!f
extern CFLOBDD MkSecond(CFLOBDD f, CFLOBDD g);       // \f.\g.g
extern CFLOBDD MkNotSecond(CFLOBDD f, CFLOBDD g);    // \f.\g.!g

// extern CFLOBDD MkPlus(CFLOBDD f, CFLOBDD g);      // \f.\g.(f + g)
// extern CFLOBDD MkTimes(CFLOBDD f, CFLOBDD g);     // \f.\g.(f * g)

// N-ary operations on CFLOBDDs --------------------------------------
extern CFLOBDD MkAnd(int N, ...);                    // \f1. ... \fk.(f1 && ... && fk)
extern CFLOBDD MkNand(int N, ...);                   // \f1. ... \fk.!(f1 && ... && fk)
extern CFLOBDD MkOr(int N, ...);                     // \f1. ... \fk.(f1 || ... || fk)
extern CFLOBDD MkNor(int N, ...);                    // \f1. ... \fk.!(f1 || ... || fk)

// Ternary operations on CFLOBDDs --------------------------------------
extern CFLOBDD MkIfThenElse(CFLOBDD f, CFLOBDD g, CFLOBDD h);  // \a.\b.\c.(a && b) || (!a && c)
extern CFLOBDD MkNegMajority(CFLOBDD f, CFLOBDD g, CFLOBDD h); // \a.\b.\c.(b && !a) || (c && !a) || (b && c)

extern CFLOBDD MkRestrict(CFLOBDD f, unsigned int i, bool val);  // \f. f | (x_i = val)
extern CFLOBDD MkExists(CFLOBDD f, unsigned int i);              // \f. exists x_i : f
extern CFLOBDD MkForall(CFLOBDD f, unsigned int i);              // \f. forall x_i : f
extern CFLOBDD MkCompose(CFLOBDD f, int i, CFLOBDD g);  // \f. f | (x_i = g)

// Other operations on CFLOBDDs ---------------------------------------------------
extern bool DependsOn(CFLOBDD f, int i);
extern bool IsPositiveCube(CFLOBDD f);
extern bool IsPositiveCubeInt(CFLOBDD f, int least);
extern bool SupportSetIs(CFLOBDD f, const apvector<int> &ss);
extern apvector<int> *GetSupportSet(CFLOBDD f);

extern double ComputeProbability(CFLOBDD f, std::vector<double>& probs); // probs is a list of probs of vars in the same order.
extern std::vector<double> ComputeProbabilityOfList(CFLOBDD f, std::vector<std::vector<double>>& probs); // probs is a list of probs of vars in the same order.
// extern std::vector<double> ComputeEntropyOfList(CFLOBDD f, std::vector<std::vector<double>>& probs); // probs is a list of probs of vars in the same order.
extern void PrintCFLOBDD(CFLOBDD f);

} // namespace CFL_OBDD


#endif
