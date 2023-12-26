
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "assignment.h"
#include "bool_op.h"
#include "weighted_cflobdd_t.h"
#include "wmatrix1234_complex_fb_mul.h"

using namespace CFL_OBDD;

//********************************************************************
// CFLOBDD_T specializations
//********************************************************************

namespace CFL_OBDD	{
	template class WEIGHTED_CFLOBDD_T<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;    // aka CFLOBDD
}

#include "weighted_cflobdd_t.h"
#include "weighted_cflobdd_top_node_complex_fb_mul.h"

// Constructors/Destructor -------------------------------------------

// Default constructor
template<>
WEIGHTED_CFLOBDD_T<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WEIGHTED_CFLOBDD_T()
{
	root = NULL;//MkTrueFloatMulTop();
}


namespace CFL_OBDD {

	// CFLOBDD-creation operations --------------------------------------

// 	// Create representation of \x.true
// 	CFLOBDD_FLOAT_BOOST MkTrue()
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkTrueTop());
// 	}

// 	// Create representation of \x.false
// 	CFLOBDD_FLOAT_BOOST MkFalse()
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkFalseTop());
// 	}

// 	// Create representation of \x.x_i
// 	CFLOBDD_FLOAT_BOOST MkProjection(unsigned int i)
// 	{
// 		assert(i < (1 << CFLOBDD_MAX_LEVEL));   // i.e., i < 2**maxLevel
// 		return CFLOBDD_FLOAT_BOOST(MkDistinction(i));
// 	}

// 	CFLOBDD_FLOAT_BOOST MkCompose(CFLOBDD_FLOAT_BOOST f, int i, CFLOBDD_FLOAT_BOOST g)
// 	{
// #ifdef NDEBUG
// 		return  MkComposeTop(f.root, i, g.root);
// #else
// 		CFLOBDD_FLOAT_BOOST comp = MkComposeTop(f.root, i, g.root);
// 		assert(comp.root->count == 1);
// 		return comp;
// #endif
// 	}

// 	// MkXBit, MkYBit, MkZBit -----------------------------------------------
// 	//
// 	// MISSING: introduce a level of indirection so that the bits can be
// 	// in a different order than interleaved.
// 	CFLOBDD_FLOAT_BOOST MkXBit(unsigned int j) {
// 		return MkProjection(j * 4);
// 	}

// 	CFLOBDD_FLOAT_BOOST MkYBit(unsigned int j) {
// 		return MkProjection(j * 4 + 1);
// 	}

// #define EVEN(x) (((x)/2) * 2 == (x))

// 	CFLOBDD_FLOAT_BOOST MkZBit(unsigned int j) {
// 		//    return MkProjection((j/2) * 4 + 2 + (EVEN(j) ? 0 : 1));
// 		return MkProjection(j * 4 + 2);
// 	}

// 	// Create a representation of the addition relation for natural numbers
// 	// { (xi yi zi _)* | vec{x} * vec{y} = vec{z} }.
// 	// We create entries only for sums where the highest bit of vec{x} and vec{y}
// 	// are 0 -- TEMPORARILY REMOVED
// 	CFLOBDD_FLOAT_BOOST MkAdditionInterleavedRecursive()
// 	{
// 		CFLOBDD_FLOAT_BOOST quasiAns = CFLOBDD_FLOAT_BOOST(MkAdditionInterleavedRecursiveTop());

// 		return quasiAns;

// 		// Restrict the highest bits of x and y to 0
// 		// CFLOBDD_FLOAT_BOOST xRestriction = MkNot(MkXBit(1 << ((CFLOBDD_MAX_LEVEL - 2) - 1)));
// 		// CFLOBDD_FLOAT_BOOST yRestriction = MkNot(MkYBit(1 << ((CFLOBDD_MAX_LEVEL - 2) - 1)));
// 		// return MkAnd(quasiAns, MkAnd(xRestriction, yRestriction));
// 	}

// 	// Create a representation of the addition relation for natural numbers
// 	// { (xi yi zi _)* | vec{x} * vec{y} = vec{z} }.
// 	// We create entries only for sums where the highest bit of vec{x} and vec{y}
// 	// are 0 -- TEMPORARILY REMOVED
// 	CFLOBDD_FLOAT_BOOST MkAdditionInterleaved()
// 	{
// 		CFLOBDD_FLOAT_BOOST quasiAns = CFLOBDD_FLOAT_BOOST(MkAdditionInterleavedTop());

// 		return quasiAns;

// 		// Restrict the highest bits of x and y to 0
// 		// CFLOBDD_FLOAT_BOOST xRestriction = MkNot(MkXBit(1 << ((CFLOBDD_MAX_LEVEL - 2) - 1)));
// 		// CFLOBDD_FLOAT_BOOST yRestriction = MkNot(MkYBit(1 << ((CFLOBDD_MAX_LEVEL - 2) - 1)));
// 		// return MkAnd(quasiAns, MkAnd(xRestriction, yRestriction));
// 	}

// 	// Create a representation of the addition relation for natural numbers
// 	// { (xi yi zi _)* | vec{x} * vec{y} = vec{z} } by brute force.
// 	// We create entries only for sums where the highest bit of vec{x} and vec{y}
// 	// are 0 -- TEMPORARILY REMOVED
// 	CFLOBDD_FLOAT_BOOST MkAdditionInterleavedBruteForce() {
// 		CFLOBDD_FLOAT_BOOST ans = MkFalse();
// 		if (2 <= CFLOBDD_MAX_LEVEL && CFLOBDD_MAX_LEVEL <= 7) {
// 			unsigned int numBits = 1 << (CFLOBDD_MAX_LEVEL - 2);   // i.e., 2**(maxLevel-2)
// 			for (unsigned int i = 0; i < ((unsigned int)(1 << (numBits))); i++) {  // TWR: (numBits - 1)
// 				CFLOBDD_FLOAT_BOOST X = MkX(i);
// 				for (unsigned int j = 0; j < ((unsigned int)(1 << (numBits))); j++) {  // TWR: (numBits - 1)
// 					CFLOBDD_FLOAT_BOOST Y = MkY(j);
// 					CFLOBDD_FLOAT_BOOST Z = MkZ(i + j);
// 					ans = MkOr(ans, MkAnd(X, MkAnd(Y, Z)));
// 					// std::cout << i << " + " << j << " = " << i + j << std::endl;
// 				}
// 			}
// 		}
// 		else {
// 			std::cerr << "Cannot use MkAdditionInterleavedBruteForce: maxLevel must be in the range [2..7]" << std::endl;
// 		}
// 		return ans;
// 	}

// 	// Bits of i -> vec{x} for addition relation
// 	CFLOBDD_FLOAT_BOOST MkX(unsigned int i) {
// 		CFLOBDD_FLOAT_BOOST ans = MkTrue();
// 		if (2 <= CFLOBDD_MAX_LEVEL && CFLOBDD_MAX_LEVEL <= 7) {
// 			unsigned int bits = 1 << (CFLOBDD_MAX_LEVEL - 2);   // i.e., 2**(maxLevel-2)
// 			unsigned int mask = 1;
// 			for (unsigned int j = 0; j < bits; j++) {
// 				if (mask & i) ans = MkAnd(ans, MkXBit(j));
// 				else          ans = MkAnd(ans, MkNot(MkXBit(j)));
// 				mask = mask << 1;
// 			}
// 		}
// 		else {
// 			std::cerr << "Cannot use MkX: maxLevel must be in the range [2..7]" << std::endl;
// 		}
// 		return ans;
// 	}

// 	// Bits of i -> vec{y} for addition relation
// 	CFLOBDD_FLOAT_BOOST MkY(unsigned int i) {
// 		CFLOBDD_FLOAT_BOOST ans = MkTrue();
// 		if (2 <= CFLOBDD_MAX_LEVEL && CFLOBDD_MAX_LEVEL <= 7) {
// 			unsigned int bits = 1 << (CFLOBDD_MAX_LEVEL - 2);   // i.e., 2**(maxLevel-2)
// 			unsigned int mask = 1;
// 			for (unsigned int j = 0; j < bits; j++) {
// 				if (mask & i) ans = MkAnd(ans, MkYBit(j));
// 				else          ans = MkAnd(ans, MkNot(MkYBit(j)));
// 				mask = mask << 1;
// 			}
// 		}
// 		else {
// 			std::cerr << "Cannot use MkY: maxLevel must be in the range [2..7]" << std::endl;
// 		}
// 		return ans;
// 	}

// 	// Bits of i -> vec{z} for addition relation
// 	CFLOBDD_FLOAT_BOOST MkZ(unsigned int i) {
// 		CFLOBDD_FLOAT_BOOST ans = MkTrue();
// 		if (2 <= CFLOBDD_MAX_LEVEL && CFLOBDD_MAX_LEVEL <= 7) {
// 			unsigned int bits = 1 << (CFLOBDD_MAX_LEVEL - 2);   // i.e., 2**(maxLevel-2)
// 			unsigned int mask = 1;
// 			for (unsigned int j = 0; j < bits; j++) {
// 				if (mask & i) ans = MkAnd(ans, MkZBit(j));
// 				else          ans = MkAnd(ans, MkNot(MkZBit(j)));
// 				mask = mask << 1;
// 			}
// 		}
// 		else {
// 			std::cerr << "Cannot use MkZ: maxLevel must be in the range [2..7]" << std::endl;
// 		}
// 		return ans;
// 	}
// }
// // Create a representation of the multiplication relation for natural numbers
// // with numBits bits: { (xi yi zi _)* | vec{x} * vec{y} = vec{z} }
// // We create entries only for products without overflow, which includes
// // products of the form x * 1 and 1 * y.
// // The current method is brute force.
// // As an optimization, we only create entries for x * y where x >= y.
// CFLOBDD_FLOAT_BOOST MkMultiplicationInterleavedBruteForce() {
// 	CFLOBDD_FLOAT_BOOST ans = MkFalse();
// 	if (2 <= CFLOBDD_MAX_LEVEL && CFLOBDD_MAX_LEVEL <= 7) {
// 		unsigned int numBits = 1 << (CFLOBDD_MAX_LEVEL - 2);   // i.e., 2**(maxLevel-2)

// 		CFLOBDD_FLOAT_BOOST ZeroZ = MkZ(0);

// 		// Add all products of the form 0 * y = 0
// 		//        std::cout << "Add all products of the form 0 * y = 0" << std::endl;
// 		//        CFLOBDD_FLOAT_BOOST ZeroXZ = MkAnd(MkX(0), ZeroZ);
// 		//        for (unsigned int j = 0; j < ((unsigned int)(1 << numBits)); j++) {
// 		//            CFLOBDD_FLOAT_BOOST Y = MkY(j);
// 		//            ans = MkOr(ans, MkAnd(ZeroXZ,Y));
// 		//        }

// 		// Add all products of the form x * 0 = 0 -- i.e., x is unconstrained
// 		std::cout << "Add all products of the form x * 0 = 0" << std::endl;
// 		ans = MkOr(ans, MkAnd(MkY(0), ZeroZ));
// 		/*        CFLOBDD_FLOAT_BOOST ZeroYZ = MkAnd(MkY(0), ZeroZ);
// 		for (unsigned int j = 0; j < ((unsigned int)(1 << numBits)); j++) {
// 		CFLOBDD_FLOAT_BOOST X = MkX(j);
// 		ans = MkOr(ans, MkAnd(X, ZeroYZ));
// 		}
// 		*/
// 		//        // Add all products of the form 1 * y = y, y in [1 .. 2^numBits - 1]
// 		//        std::cout << "Add all products of the form 1 * y = y, y in [1 .. 2^numBits - 1]" << std::endl;
// 		//        CFLOBDD_FLOAT_BOOST OneX = MkX(1);
// 		//        for (unsigned int j = 1; j < ((unsigned int)(1 << numBits)); j++) {
// 		//            CFLOBDD_FLOAT_BOOST Y = MkY(j);
// 		//            CFLOBDD_FLOAT_BOOST Z = MkZ(j);
// 		//            ans = MkOr(ans, MkAnd(OneX,MkAnd(Y,Z)));
// 		//        }

// 		// Add all products of the form x * 1 = x, x in [1 .. 2^numBits - 1]
// 		std::cout << "Add all products of the form x * 1 = x, x in [1 .. 2^numBits - 1]" << std::endl;
// 		CFLOBDD_FLOAT_BOOST XeqZ = MkTrue();
// 		for (unsigned int b = 0; b < numBits; b++) {
// 			CFLOBDD_FLOAT_BOOST Xb = MkXBit(b);
// 			CFLOBDD_FLOAT_BOOST Zb = MkZBit(b);
// 			XeqZ = MkAnd(XeqZ, MkIff(Xb, Zb));
// 		}
// 		ans = MkOr(ans, MkAnd(MkY(1), XeqZ));

// 		/*        // Add all products of the form x * 1 = x, x in [1 .. 2^numBits - 1]
// 		std::cout << "Add all products of the form x * 1 = x, x in [1 .. 2^numBits - 1]" << std::endl;
// 		CFLOBDD_FLOAT_BOOST OneY = MkY(1);
// 		for (unsigned int i = 1; i < ((unsigned int)(1 << numBits)); i++) {
// 		CFLOBDD_FLOAT_BOOST X = MkX(i);
// 		CFLOBDD_FLOAT_BOOST Z = MkZ(i);
// 		ans = MkOr(ans, MkAnd(X,MkAnd(OneY,Z)));
// 		}
// 		ans = MkOr(ans, MkAnd(MkX(0), MkAnd(OneY, MkZ(0))));
// 		*/
// 		// Add all products of all other forms that cannot overflow: i,j in [2 .. 2^numBits - 1]
// 		std::cout << "Add all products of all other forms that cannot overflow: i,j in [2 .. 2^numBits - 1]" << std::endl;
// 		for (unsigned long long i = 2ULL; i < ((unsigned long long)(1 << numBits)) / 2ULL; i++) {
// 			std::cout << "i = " << i << std::endl;
// 			CFLOBDD_FLOAT_BOOST X = MkX(i);
// 			//            for (unsigned long long j = 2ULL; j < ((unsigned long long)(1 << numBits))/2ULL; j++) {
// 			for (unsigned long long j = 2ULL; j <= i; j++) {
// 				unsigned long long k = i*j;
// 				if (k < ((unsigned long long)(1 << numBits))) { // no overflow
// 					CFLOBDD_FLOAT_BOOST Y = MkY(j);
// 					CFLOBDD_FLOAT_BOOST Z = MkZ(k);
// 					ans = MkOr(ans, MkAnd(X, MkAnd(Y, Z)));
// 					//                    std::cout << " No overflow on [" << i << ", " << j << ", " << k << "]" << std::endl;
// 				}
// 				//                else
// 				//                    std::cout << "Overflow on (" << i << ", " << j << ", " << k << ")" << std::endl;
// 			}
// 		}
// 	}
// 	else {
// 		std::cerr << "Cannot use MkMultiplicationInterleavedBruteForce: maxLevel must be in the range [2..7]" << std::endl;
// 	}
// 	return ans;
// }

// //
// // HighOrderBitPosition
// //
// // Find the position of the high-order bit of z.
// //
// static int HighOrderBitPosition(unsigned int z) {
// 	int pos = -1;
// 	while (z > 0) {
// 		z = z >> 1;
// 		pos++;
// 	}
// 	return pos;
// }

// namespace CFL_OBDD {
// 	//
// 	// Factor
// 	//
// 	// Return true if z has a non-trivial factorization w.r.t. relation R,
// 	// in which case the factors are returned in f1 and f2 (as CFLOBDDs)
// 	// and v1 and v2 (as unsigned ints).
// 	//
// 	// The search for factors of z is biased away from 1 and z.
// 	//     Search from high-order bits to low-order bits
// 	//     Let k denote the position of the high-order bit of z.
// 	//     For there to be a non-trivial factoring of z, the position
// 	//        of the high-order bit of each factor is strictly less than k;
// 	//        i.e., the k-th bit of both factors must be 0.
// 	//        Consequently, we try the 0,0 case first
// 	//
// 	bool FactorZ(CFLOBDD_FLOAT_BOOST R, unsigned int z, CFLOBDD_FLOAT_BOOST &f1, CFLOBDD_FLOAT_BOOST &f2, unsigned int &v1, unsigned int &v2) {
// 		assert(CFLOBDD_MAX_LEVEL >= 2);
// 		unsigned int numBits = 1 << (CFLOBDD_MAX_LEVEL - 2);   // i.e., 2**(maxLevel-2)

// 		CFLOBDD_FLOAT_BOOST False = MkFalse();
// 		CFLOBDD_FLOAT_BOOST temp;
// 		CFLOBDD_FLOAT_BOOST Z = MkZ(z);
// 		CFLOBDD_FLOAT_BOOST curRel = MkAnd(R, Z);  // assert that the product is z
// 		CFLOBDD_FLOAT_BOOST g1 = MkTrue();         // g1 initially unconstrained
// 		CFLOBDD_FLOAT_BOOST g2 = g1;               // g2 initially unconstrained
// 		unsigned int w1 = 0;
// 		unsigned int w2 = 0;

// 		// Loop through the bit positions; find bit values for g1 and g2
// 		int hobp = HighOrderBitPosition(z);
// 		if (4 * hobp >= (1 << CFLOBDD_MAX_LEVEL)) {
// 			f1 = MkX(0);
// 			f2 = MkY(0);
// 			v1 = 0;
// 			v2 = 0;
// 			return false;
// 		}

// 		// Zero out the high-order bits
// 		for (unsigned int j = hobp + 1; j < numBits; j++) {
// 			CFLOBDD_FLOAT_BOOST curXTrue = MkXBit(j);
// 			CFLOBDD_FLOAT_BOOST curXFalse = MkNot(curXTrue);
// 			CFLOBDD_FLOAT_BOOST curYTrue = MkYBit(j);
// 			CFLOBDD_FLOAT_BOOST curYFalse = MkNot(curYTrue);
// 			g1 = MkAnd(g1, curXFalse);
// 			g2 = MkAnd(g2, curYFalse);
// 			w1 = w1 & ~((unsigned int)(1 << j));
// 			w2 = w2 & ~((unsigned int)(1 << j));
// 		}

// 		// Search for the values of the rest of the bits, from high-order to low-order
// 		for (unsigned int j = 0; j <= hobp; j++) {
// 			unsigned int i = hobp - j;
// 			CFLOBDD_FLOAT_BOOST curXTrue = MkXBit(i);
// 			CFLOBDD_FLOAT_BOOST curXFalse = MkNot(curXTrue);
// 			CFLOBDD_FLOAT_BOOST curYTrue = MkYBit(i);
// 			CFLOBDD_FLOAT_BOOST curYFalse = MkNot(curYTrue);

// 			// One of the following four possibilities must be true
// 			// 0,0 case performed first
// 			temp = MkAnd(curRel, MkAnd(curXFalse, curYFalse));
// 			if (temp != False) {
// 				curRel = temp;
// 				g1 = MkAnd(g1, curXFalse);
// 				g2 = MkAnd(g2, curYFalse);
// 				w1 = w1 & ~((unsigned int)(1 << i));
// 				w2 = w2 & ~((unsigned int)(1 << i));
// 				continue;
// 			}

// 			// 1,0
// 			temp = MkAnd(curRel, MkAnd(curXTrue, curYFalse));
// 			if (temp != False) {
// 				curRel = temp;
// 				g1 = MkAnd(g1, curXTrue);
// 				g2 = MkAnd(g2, curYFalse);
// 				w1 = w1 | ((unsigned int)(1 << i));
// 				w2 = w2 & ~((unsigned int)(1 << i));
// 				continue;
// 			}
// 			// 0,1
// 			temp = MkAnd(curRel, MkAnd(curXFalse, curYTrue));
// 			if (temp != False) {
// 				curRel = temp;
// 				g1 = MkAnd(g1, curXFalse);
// 				g2 = MkAnd(g2, curYTrue);
// 				w1 = w1 & ~((unsigned int)(1 << i));
// 				w2 = w2 | ((unsigned int)(1 << i));
// 				continue;
// 			}
// 			// 1,1
// 			temp = MkAnd(curRel, MkAnd(curXTrue, curYTrue));
// 			if (temp != False) {
// 				curRel = temp;
// 				g1 = MkAnd(g1, curXTrue);
// 				g2 = MkAnd(g2, curYTrue);
// 				w1 = w1 | ((unsigned int)(1 << i));
// 				w2 = w2 | ((unsigned int)(1 << i));
// 				continue;
// 			}
// 			assert(false);
// 		}
// 		f1 = g1;
// 		f2 = g2;
// 		v1 = w1;
// 		v2 = w2;
// 		CFLOBDD_FLOAT_BOOST XOne = MkX(1);
// 		CFLOBDD_FLOAT_BOOST YOne = MkY(1);
// 		return (g1 != XOne && g2 != YOne);  // Return true if neither factor is 1
// 	}


// 	// Create representation of parity function
// 	CFLOBDD_FLOAT_BOOST MkParity()
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkParityTop());
// 	}



// 	// Create the representation of a step function
// 	CFLOBDD_FLOAT_BOOST MkStepUpOneFourth()
// 	{
// 		assert(CFLOBDD_MAX_LEVEL >= 1);
// 		return CFLOBDD_FLOAT_BOOST(MkStepUpOneFourthTop());
// 	}

// 	// Create the representation of a step function
// 	CFLOBDD_FLOAT_BOOST MkStepDownOneFourth()
// 	{
// 		assert(CFLOBDD_MAX_LEVEL >= 1);
// 		return CFLOBDD_FLOAT_BOOST(MkStepDownOneFourthTop());
// 	}

// #ifdef ARBITRARY_STEP_FUNCTIONS
// 	// Create the representation of a step function
// 	CFLOBDD_FLOAT_BOOST MkStepUp(unsigned int i)
// 	{
// 		assert(CFLOBDD_MAX_LEVEL <= 5);
// 		return CFLOBDD_FLOAT_BOOST(MkStepUpTop(i));
// 	}

// 	// Create the representation of a step function
// 	CFLOBDD_FLOAT_BOOST MkStepDown(unsigned int i)
// 	{
// 		assert(CFLOBDD_MAX_LEVEL <= 5);
// 		return CFLOBDD_FLOAT_BOOST(MkStepDownTop(i));
// 	}

// 	// Create the representation of an up-pulse function
// 	CFLOBDD_FLOAT_BOOST MkPulseUp(unsigned int i, unsigned int j)
// 	{
// 		assert(i < j);
// 		return MkAnd(MkStepUp(i), MkStepDown(j));
// 	}

// 	// Create the representation of a down-pulse function
// 	CFLOBDD_FLOAT_BOOST MkPulseDown(unsigned int i, unsigned int j)
// 	{
// 		assert(i < j);
// 		return MkOr(MkStepDown(i), MkStepUp(j));
// 	}

// #endif

// 	CFLOBDD_FLOAT_BOOST shiftAtoB(const CFLOBDD_FLOAT_BOOST f, const unsigned int levelAtWhichToShift)
// 	{
// 		assert(levelAtWhichToShift <= CFLOBDD_MAX_LEVEL);
// 		CFLOBDD_FLOAT_BOOST ans(shiftAtoBTop(f.root, levelAtWhichToShift));
// 		return ans;
// 	}

// 	CFLOBDD_FLOAT_BOOST shiftBtoA(const CFLOBDD_FLOAT_BOOST f, const unsigned int levelAtWhichToShift)
// 	{
// 		assert(levelAtWhichToShift <= CFLOBDD_MAX_LEVEL);
// 		return CFLOBDD_FLOAT_BOOST(shiftBtoATop(f.root, levelAtWhichToShift));
// 	}

// 	CFLOBDD_FLOAT_BOOST shiftAtoBAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, const CFLOBDD_FLOAT_BOOST f)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(shiftAtoBAtLevelOneTop(visitedNodes, totalVisitCount, redundantVisitCount, f.root));
// 	}

// 	CFLOBDD_FLOAT_BOOST shiftBtoAAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, const CFLOBDD_FLOAT_BOOST f)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(shiftBtoAAtLevelOneTop(visitedNodes, totalVisitCount, redundantVisitCount, f.root));
// 	}

// 	CFLOBDD_FLOAT_BOOST duplicateAinBAtLevelOne(const CFLOBDD_FLOAT_BOOST f)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(duplicateAinBAtLevelOneTop(f.root));
// 	}

// 	// Unary operations on CFLOBDDs --------------------------------------

// 	// Implements \f.!f
// 	CFLOBDD_FLOAT_BOOST MkNot(CFLOBDD_FLOAT_BOOST f)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkNot(f.root));
// 	}

// 	// Binary operations on CFLOBDDs --------------------------------------

// 	// \f.\g.(f && g)
// 	CFLOBDD_FLOAT_BOOST MkAnd(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkAnd(f.root, g.root));
// 	}

// 	// \f.\g.!(f && g)
// 	CFLOBDD_FLOAT_BOOST MkNand(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkNand(f.root, g.root));
// 	}

// 	// \f.\g.(f || g)
// 	CFLOBDD_FLOAT_BOOST MkOr(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkOr(f.root, g.root));
// 	}

// 	// \f.\g.!(f || g)
// 	CFLOBDD_FLOAT_BOOST MkNor(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkNor(f.root, g.root));
// 	}

// 	// \f.\g.(f == g)
// 	CFLOBDD_FLOAT_BOOST MkIff(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkIff(f.root, g.root));
// 	}

// 	// \f.\g.(f != g)
// 	CFLOBDD_FLOAT_BOOST MkExclusiveOr(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkExclusiveOr(f.root, g.root));
// 	}

// 	// \f.\g.(!f || g)
// 	CFLOBDD_FLOAT_BOOST MkImplies(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkImplies(f.root, g.root));
// 	}

// 	// \f.\g.(f && !g)
// 	CFLOBDD_FLOAT_BOOST MkMinus(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkMinus(f.root, g.root));
// 	}

// 	// \f.\g.(!g || f)
// 	CFLOBDD_FLOAT_BOOST MkQuotient(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkQuotient(f.root, g.root));
// 	}

// 	// \f.\g.(g && !f)
// 	CFLOBDD_FLOAT_BOOST MkNotQuotient(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkNotQuotient(f.root, g.root));
// 	}

// 	// \f.\g.f
// 	CFLOBDD_FLOAT_BOOST MkFirst(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkFirst(f.root, g.root));
// 	}

// 	// \f.\g.!f
// 	CFLOBDD_FLOAT_BOOST MkNotFirst(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkNotFirst(f.root, g.root));
// 	}

// 	// \f.\g.g
// 	CFLOBDD_FLOAT_BOOST MkSecond(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkSecond(f.root, g.root));
// 	}

// 	// \f.\g.!g
// 	CFLOBDD_FLOAT_BOOST MkNotSecond(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkNotSecond(f.root, g.root));
// 	}

// 	/*
// 	// \f.\g.(f + g)
// 	CFLOBDD_FLOAT_BOOST MkPlus(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 	return CFLOBDD_FLOAT_BOOST(MkPlus(f.root, g.root));
// 	}

// 	// \f.\g.(f * g)
// 	CFLOBDD_FLOAT_BOOST MkTimes(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g)
// 	{
// 	return CFLOBDD_FLOAT_BOOST(MkTimes(f.root, g.root));
// 	}
// 	*/

// 	// N-ary operations on CFLOBDD_FLOAT_BOOSTs -----------------------------------

// 	// // \f1. ... \fk.(f1 && ... && fk)
// 	// CFLOBDD_FLOAT_BOOST  MkAnd(int N, ...)
// 	// {
// 	// 	CFLOBDD_FLOAT_BOOST  temp;

// 	// 	assert(N >= 2);
// 	// 	va_list ap;
// 	// 	va_start(ap, N);
// 	// 	temp = va_arg(ap, CFLOBDD_FLOAT_BOOST);
// 	// 	for (int k = 1; k < N; k++) {
// 	// 		temp = MkAnd(temp, va_arg(ap, CFLOBDD_FLOAT_BOOST));
// 	// 	}
// 	// 	va_end(ap);
// 	// 	return temp;
// 	// }

// 	// // \f1. ... \fk.!(f1 && ... && fk)
// 	// CFLOBDD_FLOAT_BOOST  MkNand(int N, ...)
// 	// {
// 	// 	CFLOBDD_FLOAT_BOOST  temp;

// 	// 	assert(N >= 2);
// 	// 	va_list ap;
// 	// 	va_start(ap, N);
// 	// 	temp = va_arg(ap, CFLOBDD_FLOAT_BOOST);
// 	// 	for (int k = 1; k < N; k++) {
// 	// 		temp = MkAnd(temp, va_arg(ap, CFLOBDD_FLOAT_BOOST));
// 	// 	}
// 	// 	va_end(ap);
// 	// 	return MkNot(temp);
// 	// }

// 	// // \f1. ... \fk.(f1 || ... || fk)
// 	// CFLOBDD_FLOAT_BOOST  MkOr(int N, ...)
// 	// {
// 	// 	CFLOBDD_FLOAT_BOOST  temp;

// 	// 	assert(N >= 2);
// 	// 	va_list ap;
// 	// 	va_start(ap, N);
// 	// 	temp = va_arg(ap, CFLOBDD_FLOAT_BOOST);
// 	// 	for (int k = 1; k < N; k++) {
// 	// 		temp = MkOr(temp, va_arg(ap, CFLOBDD_FLOAT_BOOST));
// 	// 	}
// 	// 	va_end(ap);
// 	// 	return temp;
// 	// }

// 	// // \f1. ... \fk.!(f1 || ... || fk)
// 	// CFLOBDD_FLOAT_BOOST  MkNor(int N, ...)
// 	// {
// 	// 	CFLOBDD_FLOAT_BOOST  temp;

// 	// 	assert(N >= 2);
// 	// 	va_list ap;
// 	// 	va_start(ap, N);
// 	// 	temp = va_arg(ap, CFLOBDD_FLOAT_BOOST);
// 	// 	for (int k = 1; k < N; k++) {
// 	// 		temp = MkOr(temp, va_arg(ap, CFLOBDD_FLOAT_BOOST));
// 	// 	}
// 	// 	va_end(ap);
// 	// 	return MkNot(temp);
// 	// }

// 	// Ternary operations on CFLOBDDs --------------------------------------

// 	// \a.\b.\c.(a && b) || (!a && c)
// 	CFLOBDD_FLOAT_BOOST MkIfThenElse(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g, CFLOBDD_FLOAT_BOOST h)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkIfThenElse(f.root, g.root, h.root));
// 	}

// 	// \a.\b.\c.(b && !a) || (c && !a) || (b && c)
// 	CFLOBDD_FLOAT_BOOST MkNegMajority(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g, CFLOBDD_FLOAT_BOOST h)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkNegMajority(f.root, g.root, h.root));
// 	}

// 	// \f. f | (x_i = val)
// 	CFLOBDD_FLOAT_BOOST MkRestrict(CFLOBDD_FLOAT_BOOST f, unsigned int i, bool val)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkRestrict(f.root, i, val));
// 	}

// 	// \f. exists x_i : f
// 	CFLOBDD_FLOAT_BOOST MkExists(CFLOBDD_FLOAT_BOOST f, unsigned int i)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkExists(f.root, i));
// 	}

// 	// \f. forall x_i : f
// 	CFLOBDD_FLOAT_BOOST MkForall(CFLOBDD_FLOAT_BOOST f, unsigned int i)
// 	{
// 		return CFLOBDD_FLOAT_BOOST(MkForall(f.root, i));
// 	}

// 	// Other operations on CFLOBDD_FLOAT_BOOSTs --------------------------------------

// 	bool DependsOn(CFLOBDD_FLOAT_BOOST f, int i)
// 	{
// 		// horrible performance!
// 		if (MkRestrict(f, i, false) == f) {
// 			assert(MkRestrict(f, i, true) == f); // or this algorithm is bad
// 			return false;
// 		}
// 		assert(MkRestrict(f, i, true) != f);
// 		return true;
// 	}

// 	bool IsPositiveCube(CFLOBDD_FLOAT_BOOST f)
// 	{
// 		return IsPositiveCubeInt(f, 0);
// 	}

// 	bool IsPositiveCubeInt(CFLOBDD_FLOAT_BOOST f, int least)
// 	{
// 		// base case:
// 		if (f == MkFalse()) return false;
// 		if (f == MkTrue()) return true;

// 		// recursive step:
// 		for (int xi = least; xi < CFLOBDD_MAX_LEVEL; xi++) {
// 			if (DependsOn(f, xi)) {
// 				if (MkRestrict(f, xi, false) != MkFalse())
// 					return false;
// 				if (!IsPositiveCubeInt(MkRestrict(f, xi, true), xi + 1))
// 					return false;
// 			}
// 		}
// 		return true;
// 	}

// 	bool SupportSetIs(CFLOBDD_FLOAT_BOOST f, const apvector<int> &ss)
// 	{
// 		assert(&ss != NULL);
// 		apvector<int> *truess = GetSupportSet(f);
// 		bool rc = (ss == *truess);
// 		delete truess;
// 		return rc;
// 	}

// 	apvector<int> *GetSupportSet(CFLOBDD_FLOAT_BOOST f)
// 		// Returns a new, sorted-in-increasing-order array
// 		// containing the support set (the set of VarNums
// 		// corresponding to indices of projection functions upon which bdd
// 		// depends -- for (x0+x1*x3), it would be [0, 1, 3]).
// 		//
// 		// Delete it when you are done.
// 	{
// 		// Decent algorithm: find the index associated with each fork node on EVERY path.  kind of ugly, but better than this quickie:

// 		apvector<int> *vec = new apvector<int>;
// 		for (unsigned int i = 0; i < CFLOBDD_MAX_LEVEL; i++) {
// 			if (DependsOn(f, i)) {
// 				vec->AddToEnd(i);
// 			}
// 		}
// 		assert(vec != NULL);
// 		return vec;

// 		/*
// 		init vector
// 		for xi in ( x0 to highest projection function xn (inclusive) ) {
// 		if (bdd.DependsOn(xi)) {
// 		add xi to vector
// 		}
// 		}
// 		return vec;
// 		*/
// 	}



} // namespace CFL_OBDD

//********************************************************************
// Operations to gather statistics
//********************************************************************

// Statistics on Connections ----------------------------------

// template<>
// void CFLOBDD_FLOAT_BOOST::DumpConnections(std::ostream & out /* = std::cout */)
// {
// 	Hashset<CFLOBDDNodeHandle> *visited = new Hashset<CFLOBDDNodeHandle>;
// 	root->DumpConnections(visited, out);
// 	delete visited;
// }

// template<>
// Hashset<CFLOBDDNodeHandle> *CFLOBDD_FLOAT_BOOST::visitedDuringGroupDumpConnections = NULL;

// template<>
// void CFLOBDD_FLOAT_BOOST::GroupDumpConnections(std::ostream & out /* = std::cout */)
// {
// 	root->DumpConnections(CFLOBDD_FLOAT_BOOST::visitedDuringGroupDumpConnections, out);
// }

// // Statistics on size ----------------------------------------------------

// /*
// template<>
// void CFLOBDD_FLOAT_BOOST::CountNodesAndEdges(unsigned int &nodeCount, unsigned int &edgeCount)
// {
// Hashset<CFLOBDDNodeHandle> *visitedNodes = new Hashset<CFLOBDDNodeHandle>;
// Hashset<CFLOBDDReturnMapBody> *visitedEdges = new Hashset<CFLOBDDReturnMapBody>;
// nodeCount = 0;
// edgeCount = 0;
// root->CountNodesAndEdges(visitedNodes, visitedEdges, nodeCount, edgeCount);
// delete visitedNodes;
// delete visitedEdges;
// }
// */

// template<>
// Hashset<CFLOBDDNodeHandle> *CFLOBDD_FLOAT_BOOST::visitedNodesDuringGroupCountNodesAndEdges = NULL;
// template<>
// Hashset<CFLOBDDReturnMapBody> *CFLOBDD_FLOAT_BOOST::visitedEdgesDuringGroupCountNodesAndEdges = NULL;

// template<>
// void CFLOBDD_FLOAT_BOOST::GroupCountNodesAndEdges(unsigned int &nodeCount, unsigned int &edgeCount)
// {
// 	unsigned int returnEdgesCount, returnEdgesObjCount;
// 	root->CountNodesAndEdges(CFLOBDD_FLOAT_BOOST::visitedNodesDuringGroupCountNodesAndEdges,
// 		CFLOBDD_FLOAT_BOOST::visitedEdgesDuringGroupCountNodesAndEdges,
// 		nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount
// 		);
// }
