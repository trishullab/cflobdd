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

#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "assignment.h"
#include "bool_op.h"
#include "cflobdd_top_node_int.h"
#include "cflobdd_top_node_t.h"

using namespace CFL_OBDD;

//********************************************************************
// CFLOBDD_T specializations
//********************************************************************


#include "cflobdd_t.h"

// Constructors/Destructor -------------------------------------------

// Default constructor
template<>
CFLOBDD_T<int>::CFLOBDD_T()
{
	root = MkTrueTop();
}


namespace CFL_OBDD {

template class CFLOBDD_T<int>;    // aka CFLOBDD

// CFLOBDD-creation operations --------------------------------------

// Create representation of \x.true
CFLOBDD MkTrue(int level)
{
  return CFLOBDD(MkTrueTop(level));
}

// Create representation of \x.false
CFLOBDD MkFalse(int level)
{
  return CFLOBDD(MkFalseTop(level));
}

// Create representation of \x.x_i
CFLOBDD MkProjection(unsigned int i, int level)
{
  assert(i < (1 << CFLOBDD::maxLevel));   // i.e., i < 2**maxLevel
  return CFLOBDD(MkDistinction(i, level));
}

CFLOBDD MkCompose(CFLOBDD f, int i, CFLOBDD g)
{
#ifdef NDEBUG
  return  MkComposeTop(f.root, i, g.root);
#else
  CFLOBDD comp = MkComposeTop(f.root, i, g.root);
  assert(comp.root->count == 1);
  return comp;
#endif
}

// MkXBit, MkYBit, MkZBit -----------------------------------------------
//
// MISSING: introduce a level of indirection so that the bits can be
// in a different order than interleaved.
CFLOBDD MkXBit(unsigned int j) {
    return MkProjection(j * 4);
}

CFLOBDD MkYBit(unsigned int j) {
    return MkProjection(j * 4 + 1);
}

#define EVEN(x) (((x)/2) * 2 == (x))

CFLOBDD MkZBit(unsigned int j) {
//    return MkProjection((j/2) * 4 + 2 + (EVEN(j) ? 0 : 1));
    return MkProjection(j * 4 + 2);
}

// Create a representation of the addition relation for natural numbers
// { (xi yi zi _)* | vec{x} * vec{y} = vec{z} }.
// We create entries only for sums where the highest bit of vec{x} and vec{y}
// are 0 -- TEMPORARILY REMOVED
CFLOBDD MkAdditionInterleavedRecursive()
{
    CFLOBDD quasiAns = CFLOBDD(MkAdditionInterleavedRecursiveTop());
	
	return quasiAns;

    // Restrict the highest bits of x and y to 0
    // CFLOBDD xRestriction = MkNot(MkXBit(1 << ((CFLOBDD::maxLevel - 2) - 1)));
	// CFLOBDD yRestriction = MkNot(MkYBit(1 << ((CFLOBDD::maxLevel - 2) - 1)));
	// return MkAnd(quasiAns, MkAnd(xRestriction, yRestriction));
}

// Create a representation of the addition relation for natural numbers
// { (xi yi zi _)* | vec{x} * vec{y} = vec{z} }.
// We create entries only for sums where the highest bit of vec{x} and vec{y}
// are 0 -- TEMPORARILY REMOVED
CFLOBDD MkAdditionInterleaved()
{
	CFLOBDD quasiAns = CFLOBDD(MkAdditionInterleavedTop());

	return quasiAns;

	// Restrict the highest bits of x and y to 0
	// CFLOBDD xRestriction = MkNot(MkXBit(1 << ((CFLOBDD::maxLevel - 2) - 1)));
	// CFLOBDD yRestriction = MkNot(MkYBit(1 << ((CFLOBDD::maxLevel - 2) - 1)));
	// return MkAnd(quasiAns, MkAnd(xRestriction, yRestriction));
}

// Create a representation of the addition relation for natural numbers
// { (xi yi zi _)* | vec{x} * vec{y} = vec{z} } by brute force.
// We create entries only for sums where the highest bit of vec{x} and vec{y}
// are 0 -- TEMPORARILY REMOVED
CFLOBDD MkAdditionInterleavedBruteForce() {
    CFLOBDD ans = MkFalse();
    if (2 <= CFLOBDD::maxLevel && CFLOBDD::maxLevel <= 7) {
        unsigned int numBits = 1 << (CFLOBDD::maxLevel - 2);   // i.e., 2**(maxLevel-2)
        for (unsigned int i = 0; i < ((unsigned int)(1 << (numBits))); i++) {  // TWR: (numBits - 1)
            CFLOBDD X = MkX(i);
            for (unsigned int j = 0; j < ((unsigned int)(1 << (numBits))); j++) {  // TWR: (numBits - 1)
                CFLOBDD Y = MkY(j);
                CFLOBDD Z = MkZ(i+j);
                ans = MkOr(ans, MkAnd(X,MkAnd(Y,Z)));
				// std::cout << i << " + " << j << " = " << i + j << std::endl;
            }
        }
    }
    else {
        std::cerr << "Cannot use MkAdditionInterleavedBruteForce: maxLevel must be in the range [2..7]" << std::endl;
    }
    return ans;
}

// Bits of i -> vec{x} for addition relation
CFLOBDD MkX(unsigned int i) {
     CFLOBDD ans = MkTrue();
    if (2 <= CFLOBDD::maxLevel && CFLOBDD::maxLevel <= 7) {
        unsigned int bits = 1 << (CFLOBDD::maxLevel - 2);   // i.e., 2**(maxLevel-2)
        unsigned int mask = 1;
        for (unsigned int j = 0; j < bits; j++) {
            if (mask & i) ans = MkAnd(ans, MkXBit(j));
            else          ans = MkAnd(ans, MkNot(MkXBit(j)));
            mask = mask << 1;
        }
    }
    else {
        std::cerr << "Cannot use MkX: maxLevel must be in the range [2..7]" << std::endl;
    }
    return ans;
}

// Bits of i -> vec{y} for addition relation
CFLOBDD MkY(unsigned int i) {
     CFLOBDD ans = MkTrue();
    if (2 <= CFLOBDD::maxLevel && CFLOBDD::maxLevel <= 7) {
        unsigned int bits = 1 << (CFLOBDD::maxLevel - 2);   // i.e., 2**(maxLevel-2)
        unsigned int mask = 1;
        for (unsigned int j = 0; j < bits; j++) {
            if (mask & i) ans = MkAnd(ans, MkYBit(j));
            else          ans = MkAnd(ans, MkNot(MkYBit(j)));
            mask = mask << 1;
        }
    }
    else {
        std::cerr << "Cannot use MkY: maxLevel must be in the range [2..7]" << std::endl;
    }
    return ans;
}

// Bits of i -> vec{z} for addition relation
CFLOBDD MkZ(unsigned int i) {
     CFLOBDD ans = MkTrue();
    if (2 <= CFLOBDD::maxLevel && CFLOBDD::maxLevel <= 7) {
        unsigned int bits = 1 << (CFLOBDD::maxLevel - 2);   // i.e., 2**(maxLevel-2)
        unsigned int mask = 1;
        for (unsigned int j = 0; j < bits; j++) {
            if (mask & i) ans = MkAnd(ans, MkZBit(j));
            else          ans = MkAnd(ans, MkNot(MkZBit(j)));
            mask = mask << 1;
        }
    }
    else {
        std::cerr << "Cannot use MkZ: maxLevel must be in the range [2..7]" << std::endl;
    }
    return ans;
}
}
// Create a representation of the multiplication relation for natural numbers
// with numBits bits: { (xi yi zi _)* | vec{x} * vec{y} = vec{z} }
// We create entries only for products without overflow, which includes
// products of the form x * 1 and 1 * y.
// The current method is brute force.
// As an optimization, we only create entries for x * y where x >= y.
CFLOBDD MkMultiplicationInterleavedBruteForce() {
    CFLOBDD ans = MkFalse();
    if (2 <= CFLOBDD::maxLevel && CFLOBDD::maxLevel <= 7) {
        unsigned int numBits = 1 << (CFLOBDD::maxLevel - 2);   // i.e., 2**(maxLevel-2)

        CFLOBDD ZeroZ = MkZ(0);

        // Add all products of the form 0 * y = 0
//        std::cout << "Add all products of the form 0 * y = 0" << std::endl;
//        CFLOBDD ZeroXZ = MkAnd(MkX(0), ZeroZ);
//        for (unsigned int j = 0; j < ((unsigned int)(1 << numBits)); j++) {
//            CFLOBDD Y = MkY(j);
//            ans = MkOr(ans, MkAnd(ZeroXZ,Y));
//        }

        // Add all products of the form x * 0 = 0 -- i.e., x is unconstrained
        std::cout << "Add all products of the form x * 0 = 0" << std::endl;
        ans = MkOr(ans, MkAnd(MkY(0), ZeroZ));
/*        CFLOBDD ZeroYZ = MkAnd(MkY(0), ZeroZ);
        for (unsigned int j = 0; j < ((unsigned int)(1 << numBits)); j++) {
            CFLOBDD X = MkX(j);
            ans = MkOr(ans, MkAnd(X, ZeroYZ));
        }
*/
        //        // Add all products of the form 1 * y = y, y in [1 .. 2^numBits - 1]
//        std::cout << "Add all products of the form 1 * y = y, y in [1 .. 2^numBits - 1]" << std::endl;
//        CFLOBDD OneX = MkX(1);
//        for (unsigned int j = 1; j < ((unsigned int)(1 << numBits)); j++) {
//            CFLOBDD Y = MkY(j);
//            CFLOBDD Z = MkZ(j);
//            ans = MkOr(ans, MkAnd(OneX,MkAnd(Y,Z)));
//        }

        // Add all products of the form x * 1 = x, x in [1 .. 2^numBits - 1]
        std::cout << "Add all products of the form x * 1 = x, x in [1 .. 2^numBits - 1]" << std::endl;
        CFLOBDD XeqZ = MkTrue();
        for (unsigned int b = 0; b < numBits; b++) {
            CFLOBDD Xb = MkXBit(b);
            CFLOBDD Zb = MkZBit(b);
            XeqZ = MkAnd(XeqZ, MkIff(Xb, Zb));
        }
        ans = MkOr(ans, MkAnd(MkY(1), XeqZ));

/*        // Add all products of the form x * 1 = x, x in [1 .. 2^numBits - 1]
        std::cout << "Add all products of the form x * 1 = x, x in [1 .. 2^numBits - 1]" << std::endl;
        CFLOBDD OneY = MkY(1);
        for (unsigned int i = 1; i < ((unsigned int)(1 << numBits)); i++) {
            CFLOBDD X = MkX(i);
            CFLOBDD Z = MkZ(i);
            ans = MkOr(ans, MkAnd(X,MkAnd(OneY,Z)));
        }
        ans = MkOr(ans, MkAnd(MkX(0), MkAnd(OneY, MkZ(0))));
*/
        // Add all products of all other forms that cannot overflow: i,j in [2 .. 2^numBits - 1]
        std::cout << "Add all products of all other forms that cannot overflow: i,j in [2 .. 2^numBits - 1]" << std::endl;
        for (unsigned long long i = 2ULL; i < ((unsigned long long)(1 << numBits))/2ULL; i++) {
            std::cout << "i = " << i << std::endl;
            CFLOBDD X = MkX(i);
//            for (unsigned long long j = 2ULL; j < ((unsigned long long)(1 << numBits))/2ULL; j++) {
            for (unsigned long long j = 2ULL; j <= i; j++) {
                unsigned long long k = i*j;
                if (k < ((unsigned long long)(1 << numBits))) { // no overflow
                    CFLOBDD Y = MkY(j);
                    CFLOBDD Z = MkZ(k);
                    ans = MkOr(ans, MkAnd(X,MkAnd(Y,Z)));
//                    std::cout << " No overflow on [" << i << ", " << j << ", " << k << "]" << std::endl;
                }
//                else
//                    std::cout << "Overflow on (" << i << ", " << j << ", " << k << ")" << std::endl;
            }
        }
    }
    else {
        std::cerr << "Cannot use MkMultiplicationInterleavedBruteForce: maxLevel must be in the range [2..7]" << std::endl;
    }
    return ans;
}

//
// HighOrderBitPosition
//
// Find the position of the high-order bit of z.
//
static int HighOrderBitPosition(unsigned int z) {
    int pos = -1;  
    while (z > 0) {
        z = z >> 1;
        pos++;
    }
    return pos;
}

namespace CFL_OBDD {
//
// Factor
//
// Return true if z has a non-trivial factorization w.r.t. relation R,
// in which case the factors are returned in f1 and f2 (as CFLOBDDs)
// and v1 and v2 (as unsigned ints).
//
// The search for factors of z is biased away from 1 and z.
//     Search from high-order bits to low-order bits
//     Let k denote the position of the high-order bit of z.
//     For there to be a non-trivial factoring of z, the position
//        of the high-order bit of each factor is strictly less than k;
//        i.e., the k-th bit of both factors must be 0.
//        Consequently, we try the 0,0 case first
//
bool FactorZ(CFLOBDD R, unsigned int z, CFLOBDD &f1, CFLOBDD &f2, unsigned int &v1, unsigned int &v2) {
    assert(CFLOBDD::maxLevel >= 2);
    unsigned int numBits = 1 << (CFLOBDD::maxLevel - 2);   // i.e., 2**(maxLevel-2)
    
    CFLOBDD False = MkFalse();
    CFLOBDD temp;
    CFLOBDD Z = MkZ(z);
    CFLOBDD curRel = MkAnd(R, Z);  // assert that the product is z
    CFLOBDD g1 = MkTrue();         // g1 initially unconstrained
    CFLOBDD g2 = g1;               // g2 initially unconstrained
    unsigned int w1 = 0;
    unsigned int w2 = 0;

    // Loop through the bit positions; find bit values for g1 and g2
    int hobp = HighOrderBitPosition(z);
    if (4*hobp >= (1 << CFLOBDD::maxLevel)) {
        f1 = MkX(0);
        f2 = MkY(0);
        v1 = 0;
        v2 = 0;
        return false;
    }

    // Zero out the high-order bits
    for (unsigned int j = hobp+1; j < numBits; j++) {
        CFLOBDD curXTrue = MkXBit(j);
        CFLOBDD curXFalse = MkNot(curXTrue);
        CFLOBDD curYTrue = MkYBit(j);
        CFLOBDD curYFalse = MkNot(curYTrue);
        g1 = MkAnd(g1, curXFalse);
        g2 = MkAnd(g2, curYFalse);
        w1 = w1 & ~((unsigned int)(1 << j));
        w2 = w2 & ~((unsigned int)(1 << j));
    }

    // Search for the values of the rest of the bits, from high-order to low-order
    for (unsigned int j = 0; j <= hobp; j++) {
        unsigned int i = hobp - j;
        CFLOBDD curXTrue = MkXBit(i);
        CFLOBDD curXFalse = MkNot(curXTrue);
        CFLOBDD curYTrue = MkYBit(i);
        CFLOBDD curYFalse = MkNot(curYTrue);
        
        // One of the following four possibilities must be true
        // 0,0 case performed first
        temp = MkAnd(curRel, MkAnd(curXFalse, curYFalse));
        if (temp != False) {
            curRel = temp;
            g1 = MkAnd(g1, curXFalse);
            g2 = MkAnd(g2, curYFalse);
            w1 = w1 & ~((unsigned int)(1 << i));
            w2 = w2 & ~((unsigned int)(1 << i));
            continue;
        }

        // 1,0
        temp = MkAnd(curRel, MkAnd(curXTrue, curYFalse));
        if (temp != False) {
            curRel = temp;
            g1 = MkAnd(g1, curXTrue);
            g2 = MkAnd(g2, curYFalse);
            w1 = w1 | ((unsigned int)(1 << i));
            w2 = w2 & ~((unsigned int)(1 << i));
            continue;
        }
        // 0,1
        temp = MkAnd(curRel, MkAnd(curXFalse, curYTrue));
        if (temp != False) {
            curRel = temp;
            g1 = MkAnd(g1, curXFalse);
            g2 = MkAnd(g2, curYTrue);
            w1 = w1 & ~((unsigned int)(1 << i));
            w2 = w2 | ((unsigned int)(1 << i));
            continue;
        }
        // 1,1
        temp = MkAnd(curRel, MkAnd(curXTrue, curYTrue));
        if (temp != False) {
            curRel = temp;
            g1 = MkAnd(g1, curXTrue);
            g2 = MkAnd(g2, curYTrue);
            w1 = w1 | ((unsigned int)(1 << i));
            w2 = w2 | ((unsigned int)(1 << i));
            continue;
        }
        assert(false);
    }
    f1 = g1;
    f2 = g2;
    v1 = w1;
    v2 = w2;
    CFLOBDD XOne = MkX(1);
    CFLOBDD YOne = MkY(1);
    return (g1 != XOne && g2 != YOne);  // Return true if neither factor is 1
}


// Create representation of parity function
CFLOBDD MkParity()
{
  return CFLOBDD(MkParityTop());
}



// Create the representation of a step function
CFLOBDD MkStepUpOneFourth()
{
  assert(CFLOBDD::maxLevel >= 1);
  return CFLOBDD(MkStepUpOneFourthTop());
}

// Create the representation of a step function
CFLOBDD MkStepDownOneFourth()
{
  assert(CFLOBDD::maxLevel >= 1);
  return CFLOBDD(MkStepDownOneFourthTop());
}

#ifdef ARBITRARY_STEP_FUNCTIONS
// Create the representation of a step function
CFLOBDD MkStepUp(unsigned int i)
{
  assert(CFLOBDD::maxLevel <= 5);
  return CFLOBDD(MkStepUpTop(i));
}

// Create the representation of a step function
CFLOBDD MkStepDown(unsigned int i)
{
  assert(CFLOBDD::maxLevel <= 5);
  return CFLOBDD(MkStepDownTop(i));
}

// Create the representation of an up-pulse function
CFLOBDD MkPulseUp(unsigned int i, unsigned int j)
{
  assert(i < j);
  return MkAnd(MkStepUp(i), MkStepDown(j));
}

// Create the representation of a down-pulse function
CFLOBDD MkPulseDown(unsigned int i, unsigned int j)
{
  assert(i < j);
  return MkOr(MkStepDown(i), MkStepUp(j));
}

#endif

CFLOBDD shiftAtoB(const CFLOBDD f, const unsigned int levelAtWhichToShift)
{
	assert(levelAtWhichToShift <= CFLOBDD::maxLevel);
	CFLOBDD ans(shiftAtoBTop(f.root, levelAtWhichToShift));
	return ans;
}

CFLOBDD shiftBtoA(const CFLOBDD f, const unsigned int levelAtWhichToShift)
{
	assert(levelAtWhichToShift <= CFLOBDD::maxLevel);
	return CFLOBDD(shiftBtoATop(f.root, levelAtWhichToShift));
}

CFLOBDD shiftAtoBAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, const CFLOBDD f)
{
	return CFLOBDD(shiftAtoBAtLevelOneTop(visitedNodes, totalVisitCount, redundantVisitCount, f.root));
}

CFLOBDD shiftBtoAAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, const CFLOBDD f)
{
	return CFLOBDD(shiftBtoAAtLevelOneTop(visitedNodes, totalVisitCount, redundantVisitCount, f.root));
}

CFLOBDD duplicateAinBAtLevelOne(const CFLOBDD f)
{
	return CFLOBDD(duplicateAinBAtLevelOneTop(f.root));
}

// Unary operations on CFLOBDDs --------------------------------------

// Implements \f.!f
CFLOBDD MkNot(CFLOBDD f)
{
  return CFLOBDD(MkNot(f.root));
}

// Binary operations on CFLOBDDs --------------------------------------

// \f.\g.(f && g)
CFLOBDD MkAnd(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkAnd(f.root, g.root));
}

// \f.\g.!(f && g)
CFLOBDD MkNand(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkNand(f.root, g.root));
}

// \f.\g.(f || g)
CFLOBDD MkOr(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkOr(f.root, g.root));
}

// \f.\g.!(f || g)
CFLOBDD MkNor(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkNor(f.root, g.root));
}

// \f.\g.(f == g)
CFLOBDD MkIff(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkIff(f.root, g.root));
}

// \f.\g.(f != g)
CFLOBDD MkExclusiveOr(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkExclusiveOr(f.root, g.root));
}

// \f.\g.(!f || g)
CFLOBDD MkImplies(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkImplies(f.root, g.root));
}

// \f.\g.(f && !g)
CFLOBDD MkMinus(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkMinus(f.root, g.root));
}

// \f.\g.(!g || f)
CFLOBDD MkQuotient(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkQuotient(f.root, g.root));
}

// \f.\g.(g && !f)
CFLOBDD MkNotQuotient(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkNotQuotient(f.root, g.root));
}

// \f.\g.f
CFLOBDD MkFirst(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkFirst(f.root, g.root));
}

// \f.\g.!f
CFLOBDD MkNotFirst(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkNotFirst(f.root, g.root));
}

// \f.\g.g
CFLOBDD MkSecond(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkSecond(f.root, g.root));
}

// \f.\g.!g
CFLOBDD MkNotSecond(CFLOBDD f, CFLOBDD g)
{
  return CFLOBDD(MkNotSecond(f.root, g.root));
}

/*
// \f.\g.(f + g)
CFLOBDD MkPlus(CFLOBDD f, CFLOBDD g)
{
	return CFLOBDD(MkPlus(f.root, g.root));
}

// \f.\g.(f * g)
CFLOBDD MkTimes(CFLOBDD f, CFLOBDD g)
{
	return CFLOBDD(MkTimes(f.root, g.root));
}
*/

// // N-ary operations on CFLOBDDs -----------------------------------

// // \f1. ... \fk.(f1 && ... && fk)
// CFLOBDD  MkAnd(int N, ...)
// {
// 	CFLOBDD  temp;

// 	assert(N >= 2);
// 	va_list ap;
// 	va_start(ap, N);
// 	temp = va_arg(ap, CFLOBDD);
// 	for (int k = 1; k < N; k++) {
// 		temp = MkAnd(temp, va_arg(ap, CFLOBDD));
// 	}
// 	va_end(ap);
// 	return temp;
// }

// // \f1. ... \fk.!(f1 && ... && fk)
// CFLOBDD  MkNand(int N, ...)
// {
// 	CFLOBDD  temp;

// 	assert(N >= 2);
// 	va_list ap;
// 	va_start(ap, N);
// 	temp = va_arg(ap, CFLOBDD);
// 	for (int k = 1; k < N; k++) {
// 		temp = MkAnd(temp, va_arg(ap, CFLOBDD));
// 	}
// 	va_end(ap);
// 	return MkNot(temp);
// }

// // \f1. ... \fk.(f1 || ... || fk)
// CFLOBDD  MkOr(int N, ...)
// {
// 	CFLOBDD  temp;

// 	assert(N >= 2);
// 	va_list ap;
// 	va_start(ap, N);
// 	temp = va_arg(ap, CFLOBDD);
// 	for (int k = 1; k < N; k++) {
// 		temp = MkOr(temp, va_arg(ap, CFLOBDD));
// 	}
// 	va_end(ap);
// 	return temp;
// }

// // \f1. ... \fk.!(f1 || ... || fk)
// CFLOBDD  MkNor(int N, ...)
// {
// 	CFLOBDD  temp;

// 	assert(N >= 2);
// 	va_list ap;
// 	va_start(ap, N);
// 	temp = va_arg(ap, CFLOBDD);
// 	for (int k = 1; k < N; k++) {
// 		temp = MkOr(temp, va_arg(ap, CFLOBDD));
// 	}
// 	va_end(ap);
// 	return MkNot(temp);
// }

// Ternary operations on CFLOBDDs --------------------------------------

// \a.\b.\c.(a && b) || (!a && c)
CFLOBDD MkIfThenElse(CFLOBDD f, CFLOBDD g, CFLOBDD h)
{
  return CFLOBDD(MkIfThenElse(f.root, g.root, h.root));
}

// \a.\b.\c.(b && !a) || (c && !a) || (b && c)
CFLOBDD MkNegMajority(CFLOBDD f, CFLOBDD g, CFLOBDD h)
{
  return CFLOBDD(MkNegMajority(f.root, g.root, h.root));
}

// \f. f | (x_i = val)
CFLOBDD MkRestrict(CFLOBDD f, unsigned int i, bool val)
{
  return CFLOBDD(MkRestrict(f.root, i, val));
}

// \f. exists x_i : f
CFLOBDD MkExists(CFLOBDD f, unsigned int i)
{
  return CFLOBDD(MkExists(f.root, i));
}

// \f. forall x_i : f
CFLOBDD MkForall(CFLOBDD f, unsigned int i)
{
  return CFLOBDD(MkForall(f.root, i));
}

// Other operations on CFLOBDDs --------------------------------------

bool DependsOn(CFLOBDD f, int i)
{
	// horrible performance!
	if (MkRestrict(f, i, false) == f) {
		assert(MkRestrict(f, i, true) == f); // or this algorithm is bad
		return false;
	}
	assert(MkRestrict(f, i, true) != f);
	return true;
}

bool IsPositiveCube(CFLOBDD f)
{
	return IsPositiveCubeInt(f, 0);
}

bool IsPositiveCubeInt(CFLOBDD f, int least)
{
	// base case:
	if (f == MkFalse()) return false;
	if (f == MkTrue()) return true;

	// recursive step:
	for (int xi = least; xi < CFLOBDD::maxLevel; xi++) {
		if (DependsOn(f, xi)) {
			if (MkRestrict(f, xi, false) != MkFalse())
				return false;
			if (!IsPositiveCubeInt(MkRestrict(f, xi, true), xi + 1))
				return false;
		}
	}
	return true;
}

bool SupportSetIs(CFLOBDD f, const apvector<int> &ss)
{
	assert(&ss != NULL);
	apvector<int> *truess = GetSupportSet(f);
	bool rc = (ss == *truess);
	delete truess;
	return rc;
}

apvector<int> *GetSupportSet(CFLOBDD f)
// Returns a new, sorted-in-increasing-order array
// containing the support set (the set of VarNums
// corresponding to indices of projection functions upon which bdd
// depends -- for (x0+x1*x3), it would be [0, 1, 3]).
//
// Delete it when you are done.
{
	// Decent algorithm: find the index associated with each fork node on EVERY path.  kind of ugly, but better than this quickie:

	apvector<int> *vec = new apvector<int>;
	for (unsigned int i = 0; i < CFLOBDD::maxLevel; i++) {
		if (DependsOn(f, i)) {
			vec->AddToEnd(i);
		}
	}
	assert(vec != NULL);
	return vec;

	/*
	init vector
	for xi in ( x0 to highest projection function xn (inclusive) ) {
	if (bdd.DependsOn(xi)) {
	add xi to vector
	}
	}
	return vec;
	*/
}

double ComputeProbability(CFLOBDD f, std::vector<double>& probs)
{
  return ComputeProbabilityTop(f.root, probs);
}

std::vector<double> ComputeProbabilityOfList(CFLOBDD f, std::vector<std::vector<double>>& probs)
{
  return ComputeProbabilityOfListTop(f.root, probs);
}

// std::vector<double> ComputeEntropyOfList(CFLOBDD f, std::vector<std::vector<double>>& probs)
// {
//   return ComputeEntropyOfListTop(f.root, probs);
// }

void PrintCFLOBDD(CFLOBDD f)
{
  std::cout << f << std::endl;
}	


} // namespace CFL_OBDD

//********************************************************************
// Operations to gather statistics
//********************************************************************

// Statistics on Connections ----------------------------------

template<>
void CFLOBDD::DumpConnections(std::ostream & out /* = std::cout */)
{
  Hashset<CFLOBDDNodeHandle> *visited = new Hashset<CFLOBDDNodeHandle>;
  root->DumpConnections(visited, out);
  delete visited;
}

template<>
Hashset<CFLOBDDNodeHandle> *CFLOBDD::visitedDuringGroupDumpConnections = NULL;

void GroupDumpConnectionsStart()
{
  CFLOBDD::visitedDuringGroupDumpConnections = new Hashset<CFLOBDDNodeHandle>;
}

void GroupDumpConnectionsEnd()
{
  delete CFLOBDD::visitedDuringGroupDumpConnections;
  CFLOBDD::visitedDuringGroupDumpConnections = NULL;
}

template<>
void CFLOBDD::GroupDumpConnections(std::ostream & out /* = std::cout */)
{
  root->DumpConnections(CFLOBDD::visitedDuringGroupDumpConnections, out);
}

// Statistics on size ----------------------------------------------------

/*
template<>
void CFLOBDD::CountNodesAndEdges(unsigned int &nodeCount, unsigned int &edgeCount)
{
  Hashset<CFLOBDDNodeHandle> *visitedNodes = new Hashset<CFLOBDDNodeHandle>;
  Hashset<CFLOBDDReturnMapBody> *visitedEdges = new Hashset<CFLOBDDReturnMapBody>;
  nodeCount = 0;
  edgeCount = 0;
  root->CountNodesAndEdges(visitedNodes, visitedEdges, nodeCount, edgeCount);
  delete visitedNodes;
  delete visitedEdges;
}
*/

template<>
Hashset<CFLOBDDNodeHandle> *CFLOBDD::visitedNodesDuringGroupCountNodesAndEdges = NULL;
template<>
Hashset<CFLOBDDReturnMapBody> *CFLOBDD::visitedEdgesDuringGroupCountNodesAndEdges = NULL;

void CFL_OBDD::GroupCountNodesAndEdgesStart(unsigned int &nodeCount, unsigned int &edgeCount)
{
  nodeCount = 0;
  edgeCount = 0;
  CFLOBDD::visitedNodesDuringGroupCountNodesAndEdges = new Hashset<CFLOBDDNodeHandle>;
  CFLOBDD::visitedEdgesDuringGroupCountNodesAndEdges = new Hashset<CFLOBDDReturnMapBody>;
}

void CFL_OBDD::GroupCountNodesAndEdgesEnd()
{
  delete CFLOBDD::visitedNodesDuringGroupCountNodesAndEdges;
  delete CFLOBDD::visitedEdgesDuringGroupCountNodesAndEdges;
  CFLOBDD::visitedNodesDuringGroupCountNodesAndEdges = NULL;
  CFLOBDD::visitedEdgesDuringGroupCountNodesAndEdges = NULL;
}

template<>
void CFLOBDD::GroupCountNodesAndEdges(unsigned int &nodeCount, unsigned int &edgeCount)
{
    unsigned int returnEdgesCount, returnEdgesObjCount;
  root->CountNodesAndEdges(CFLOBDD::visitedNodesDuringGroupCountNodesAndEdges,
                           CFLOBDD::visitedEdgesDuringGroupCountNodesAndEdges,
                           nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount
                          );
}
