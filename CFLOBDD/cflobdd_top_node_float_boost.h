#ifndef CFLOBDD_TOP_NODE_FLOAT_BOOST_GUARD
#define CFLOBDD_TOP_NODE_FLOAT_BOOST_GUARD

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


// #include <iostream>
// #include <fstream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "cflobdd_top_node_t.h"

namespace mp = boost::multiprecision;

namespace CFL_OBDD {

	typedef mp::cpp_dec_float_100 BIG_FLOAT;
	//typedef mp::number<mp::cpp_dec_float<1000> > BIG_FLOAT;
	typedef CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeFloatBoostRefPtr;
	typedef ReturnMapHandle<BIG_FLOAT> FloatBoostReturnMapHandle;
	typedef CFLOBDDTopNodeT<BIG_FLOAT> CFLOBDDFloatBoostTopNode;

	// CFLOBDDTopNode-creation operations --------------------------------------
	extern CFLOBDDTopNodeFloatBoostRefPtr MkTrueFloatTop();                    // Representation of \x.true
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkFalseTop();                   // Representation of \x.false
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkDistinction(unsigned int i);  // Representation of \x.x_i
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkAdditionInterleavedRecursiveTop();     // Representation of addition relation, created recursively
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkAdditionInterleavedTop();     // Representation of addition relation
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkParityTop();                  // Representation of parity function
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkStepUpOneFourthTop();         // Representation of step function
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkStepDownOneFourthTop();       // Representation of step function
// #ifdef ARBITRARY_STEP_FUNCTIONS
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkStepUpTop(unsigned int i);    // Representation of step function
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkStepDownTop(unsigned int i);  // Representation of step function
// #endif

// 	extern CFLOBDDTopNodeFloatBoostRefPtr shiftAtoBTop(CFLOBDDTopNodeFloatBoostRefPtr f, const unsigned int levelAtWhichToShift);
// 	extern CFLOBDDTopNodeFloatBoostRefPtr shiftBtoATop(CFLOBDDTopNodeFloatBoostRefPtr f, const unsigned int levelAtWhichToShift);

// 	extern CFLOBDDTopNodeFloatBoostRefPtr shiftAtoBAtLevelOneTop(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, CFLOBDDTopNodeFloatBoostRefPtr f);
// 	extern CFLOBDDTopNodeFloatBoostRefPtr shiftBtoAAtLevelOneTop(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, CFLOBDDTopNodeFloatBoostRefPtr f);
// 	extern CFLOBDDTopNodeFloatBoostRefPtr duplicateAinBAtLevelOneTop(CFLOBDDTopNodeFloatBoostRefPtr f);

// 	// Unary operations on CFLOBDDTopNodes --------------------------------------
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkNot(CFLOBDDTopNodeFloatBoostRefPtr f);               // \f.!f

// 	// Binary operations on CFLOBDDTopNodes --------------------------------------
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkAnd(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);          // \f.\g.(f && g)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkNand(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);         // \f.\g.!(f && g)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkOr(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);           // \f.\g.(f || g)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkNor(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);          // \f.\g.!(f || g)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkIff(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);          // \f.\g.(f == g)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkExclusiveOr(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);  // \f.\g.(f != g)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkImplies(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);      // \f.\g.(!f || g)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkMinus(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);        // \f.\g.(f && !g)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkQuotient(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);     // \f.\g.(!g || f)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkNotQuotient(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);  // \f.\g.(g && !f)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkFirst(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);        // \f.\g.f
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkNotFirst(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);     // \f.\g.!f
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkSecond(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);       // \f.\g.g
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkNotSecond(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);    // \f.\g.!g

// 	// extern CFLOBDDTopNodeFloatBoostRefPtr MkPlus(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);       // \f.\g.(f + g)
// 	// extern CFLOBDDTopNodeFloatBoostRefPtr MkTimes(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g);       // \f.\g.(f * g)

// 	// Ternary operations on CFLOBDDTopNodes ------------------------------------

// 	// \a.\b.\c.(a && b) || (!a && c)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkIfThenElse(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g, CFLOBDDTopNodeFloatBoostRefPtr h);

// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkRestrict(CFLOBDDTopNodeFloatBoostRefPtr n, unsigned int i, bool val);

// 	// \a.\b.\c.(b && !a) || (c && !a) || (b && c)
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkNegMajority(CFLOBDDTopNodeFloatBoostRefPtr f, CFLOBDDTopNodeFloatBoostRefPtr g, CFLOBDDTopNodeFloatBoostRefPtr h);

// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkExists(CFLOBDDTopNodeFloatBoostRefPtr f, unsigned int i);              // \f. exists x_i : f
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkForall(CFLOBDDTopNodeFloatBoostRefPtr f, unsigned int i);              // \f. forall x_i : f
// 	extern CFLOBDDTopNodeFloatBoostRefPtr MkComposeTop(CFLOBDDTopNodeFloatBoostRefPtr f, int i, CFLOBDDTopNodeFloatBoostRefPtr g);              // \f. f | x_i = g

} // namespace CFL_OBDD

#endif
