#ifndef MULTIPLICATION_CRT_GUARD
#define MULTIPLICATION_CRT_GUARD

//
//    Copyright (c) 2024 Thomas W. Reps
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

// Multiplication via Chinese Remainder Theorem for CFLOBDDs
// Based on "Multiplication_via_CRT.pdf"
//
// This implements a space-efficient representation of the (w x w)-bit
// multiplication relation using the Chinese Remainder Theorem. Instead of
// a single CFLOBDD that would be exponentially large, we use a collection
// of CFLOBDDs, each representing multiplication modulo a prime.
//
// Key differences from interleaved multiplication:
// - Variables are NOT interleaved: all x bits, then all y bits (concatenated)
// - Uses modular arithmetic with coprime moduli
// - Bit order is high-order to low-order

#include "cflobdd_int.h"

namespace CFL_OBDD {

// Constants from the document
const unsigned int logLogOfMaxModulus = 3;  // Each modulus must be <= 256
const unsigned int maxModulus = 1 << (1 << logLogOfMaxModulus);  // = 256

// To indicate the role of a Grouping in a recursive call (Figure 7, line 13)
enum Position { TopLevel, A, B };

// First 26 odd primes (Figure 12)
const unsigned int numberOfMultRelations = 26;
extern const unsigned int Moduli[numberOfMultRelations];

//
// ProtoCFLOBDDNumsModK
//
// Utility procedure for construction of a Grouping that represents numbers mod k.
// (Figure 9 in the document)
//
// Parameters:
//   lev - the level of the grouping to construct
//   k   - the modulus
//
// Returns:
//   A CFLOBDDNodeHandle representing the grouping
//
extern CFLOBDDNodeHandle ProtoCFLOBDDNumsModK(unsigned int lev, unsigned int k);

//
// NumsModK
//
// Construction of a CFLOBDD that represents numbers mod k (either in the
// topmost Grouping's A-connection or B-connection, as specified by the
// Position argument). (Figure 8 in the document)
//
// After following (in the A- or B-connection, as specified) the assignment a
// for some number q, where the earliest bit in a corresponds to the
// most-significant bit of q, the exit-vertex position is (q mod k).
//
// Parameters:
//   k - the modulus (must be in range [2, maxModulus])
//   p - position (A or B, not TopLevel)
//
// Returns:
//   A CFLOBDD with return values 0..k-1 representing the value mod k
//
extern CFLOBDD NumsModK(unsigned int k, Position p);

//
// MultModK
//
// Construction of a CFLOBDD that represents the multiplication relation mod k.
// The bits of the two arguments are concatenated (not interleaved).
// (Figure 10 in the document)
//
// Parameters:
//   k - the modulus
//
// Returns:
//   A CFLOBDD representing { (x, y) | (x * y) mod k }
//
extern CFLOBDD MultModK(unsigned int k);

//
// MultRelation class
//
// Represents the multiplication relation via the Chinese Remainder Theorem.
// (Figures 11-12 in the document)
//
// With the first 26 odd primes, we can represent 133-bit numbers, which is
// sufficient for handling multiplication mod 2^64. Requires CFLOBDDMaxLevel >= 9.
//
class MultRelation {
public:
    static const unsigned int numBits = 1 << (CFLOBDDMaxLevel - 1);
    unsigned int ModuliArray[numberOfMultRelations];
    CFLOBDD ModularMultRelations[numberOfMultRelations];

    MultRelation();  // Constructor
    bool operator==(const MultRelation& other) const;
    static bool VerifyShiftAndAddMultiplication();
};

//
// Factor
//
// Operation to identify all factorizations of a number v.
// (Figure 13 in the document)
//
// If v is a prime, the CFLOBDD returned has exactly two assignments that lead
// to the terminal 1: [x -> 1, y -> v] and [x -> v, y -> 1].
//
// If v is composite, in the CFLOBDD returned, the assignments that lead to
// the terminal 1 are the ones with the property [x -> a, y -> b] such that
// a * b = v mod n, where n is the product of the moduli in MultRelation R.
//
// Parameters:
//   v - the number to factor
//
// Returns:
//   A CFLOBDD representing all factorizations of v
//
extern CFLOBDD FactorViaCRT(unsigned int v);


// -----------------------------------------------------------------------------
// ShiftAndAddMultiplicationModK
//
// Operation to build the multiplication relation mod k by a shift-and-add construction
// -----------------------------------------------------------------------------
extern CFLOBDD ShiftAndAddMultiplicationModK(unsigned int k);

// -----------------------------------------------------------------------------
// VerifyShiftAndAddMultiplicationModK
//
// Operation to check whether a shift-and-add multiplier produces the correct result
//
// Parameters:
//   k - the modulus to use
//
// Returns:
//   A bool indicating whether the result == the specification (i.e., 1 means "correct")
// -----------------------------------------------------------------------------
extern bool VerifyShiftAndAddMultiplicationModK(unsigned int k);

} // namespace CFL_OBDD

#endif
