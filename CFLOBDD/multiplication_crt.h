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
#include <cstdint>
#include <boost/multiprecision/cpp_int.hpp>

namespace CFL_OBDD {

// Constants from the document
const unsigned int logLogOfMaxModulus = 4;  // Each modulus must be <= 65536
const unsigned int maxModulus = 1 << (1 << logLogOfMaxModulus);  // = 65536

// To indicate the role of a Grouping in a recursive call (Figure 7, line 13)
enum Position { TopLevel, A, B, AA, AB, BA, BB };

// =============================================================================
// Bit-width configuration
// =============================================================================

#ifndef NUM_BITS
#define NUM_BITS 2048   // BITWIDTH
#endif

// =============================================================================
// Constants for CRT-based multiplication (Figure 12 in Multiplication_via_CRT.pdf)
// =============================================================================

// Number of primes needed per bit-width (product must exceed 2^(2*bitwidth))
#if NUM_BITS == 32
  const unsigned int numberOfMultRelations = 16;
#elif NUM_BITS == 64
  const unsigned int numberOfMultRelations = 26;
#elif NUM_BITS == 128
  const unsigned int numberOfMultRelations = 44;
#elif NUM_BITS == 256
  const unsigned int numberOfMultRelations = 75;
#elif NUM_BITS == 512
  const unsigned int numberOfMultRelations = 131;
#elif NUM_BITS == 1024
  const unsigned int numberOfMultRelations = 233;
#elif NUM_BITS == 2048
  const unsigned int numberOfMultRelations = 418;
#elif NUM_BITS == 4096
  const unsigned int numberOfMultRelations = 758;
#else
  #error "Unsupported NUM_BITS value. Use 32, 64, 128, 256, 512, 1024, 2048, or 4096."
#endif

// Number of levels needed per bit-width (2**(level))
#if NUM_BITS == 32
  const unsigned int virtualMaxLevel = 6;
#elif NUM_BITS == 64
  const unsigned int virtualMaxLevel = 7;
#elif NUM_BITS == 128
  const unsigned int virtualMaxLevel = 8;
#elif NUM_BITS == 256
  const unsigned int virtualMaxLevel = 9;
#elif NUM_BITS == 512
  const unsigned int virtualMaxLevel = 10;
#elif NUM_BITS == 1024
  const unsigned int virtualMaxLevel = 11;
#elif NUM_BITS == 2048
  const unsigned int virtualMaxLevel = 12;
#elif NUM_BITS == 4096
  const unsigned int virtualMaxLevel = 13;
#else
  #error "Unsupported NUM_BITS value. Use 32, 64, 128, 256, 512, 1024, 2048, or 4096."
#endif

// Compile-time check: virtualMaxLevel must fit within CFLOBDDMaxLevel
static_assert(virtualMaxLevel <= CFLOBDDMaxLevel,
    "NUM_BITS requires virtualMaxLevel > CFLOBDD_MAX_LEVEL. "
    "Increase CFLOBDD_MAX_LEVEL in cflobdd_node.h.");

// Topmost-embedding constants (work for any CFLOBDDMaxLevel >= virtualMaxLevel)
const unsigned int bottomLevel = CFLOBDDMaxLevel - virtualMaxLevel;
const unsigned int stride = 1u << bottomLevel;
const unsigned int totalVars = 1u << CFLOBDDMaxLevel;
const unsigned int halfVars = 1u << (CFLOBDDMaxLevel - 1);

/// Total number of odd primes available
const unsigned int numberOfOddPrimes = 999;

/// Array of first 999 odd primes
extern const unsigned int AllModuli[numberOfOddPrimes];

/// Pointer to the moduli array (for backward compatibility)
extern const unsigned int* Moduli;


// =============================================================================
// INPUT_TYPE and OUTPUT_TYPE for varying bit-widths
// =============================================================================

namespace mp = boost::multiprecision;

#if NUM_BITS == 32
  typedef uint32_t INPUT_TYPE;
  typedef uint64_t OUTPUT_TYPE;
#elif NUM_BITS == 64
  typedef uint64_t INPUT_TYPE;
  typedef mp::uint128_t OUTPUT_TYPE;
#elif NUM_BITS == 128
  typedef mp::uint128_t INPUT_TYPE;
  typedef mp::uint256_t OUTPUT_TYPE;
#elif NUM_BITS == 256
  typedef mp::uint256_t INPUT_TYPE;
  typedef mp::uint512_t OUTPUT_TYPE;
#elif NUM_BITS == 512
  typedef mp::uint512_t INPUT_TYPE;
  typedef mp::uint1024_t OUTPUT_TYPE;
#elif NUM_BITS == 1024
  typedef mp::uint1024_t INPUT_TYPE;
  typedef mp::number<mp::cpp_int_backend<2048, 2048,
      mp::unsigned_magnitude, mp::unchecked, void>> OUTPUT_TYPE;
#elif NUM_BITS == 2048
  typedef mp::number<mp::cpp_int_backend<2048, 2048,
      mp::unsigned_magnitude, mp::unchecked, void>> INPUT_TYPE;
  typedef mp::number<mp::cpp_int_backend<4096, 4096,
      mp::unsigned_magnitude, mp::unchecked, void>> OUTPUT_TYPE;
#elif NUM_BITS == 4096
  typedef mp::number<mp::cpp_int_backend<4096, 4096,
      mp::unsigned_magnitude, mp::unchecked, void>> INPUT_TYPE;
  typedef mp::number<mp::cpp_int_backend<8192, 8192,
      mp::unsigned_magnitude, mp::unchecked, void>> OUTPUT_TYPE;
#endif


//
// PrintSize
//
// Utility procedure to print the number of groupings and edges
//
extern void PrintSize(CFLOBDD C);

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
// Uses the first numberOfMultRelations odd primes, whose product exceeds
// 2^(2*NUM_BITS), which is sufficient for handling NUM_BITS-bit multiplication.
//
class MultRelation {
public:
    static const unsigned int numBits = NUM_BITS;
    unsigned int ModuliArray[numberOfMultRelations];
    CFLOBDD ModularMultRelations[numberOfMultRelations];

    MultRelation();  // Constructor
    bool operator==(const MultRelation& other) const;
    static bool VerifyShiftAndAddMultiplication();

    // Multiply two m-bit values using CRT and Garner's algorithm
    // Return the 2*m-bit product
    OUTPUT_TYPE Multiply(INPUT_TYPE a, INPUT_TYPE b) const;

    // Create an assignment mapping a_val to the A-position variables
    // and b_val to the B-position variables, using stride-based embedding
    static SH_OBDD::Assignment MakeAssignment(INPUT_TYPE a_val, INPUT_TYPE b_val);
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

// -----------------------------------------------------------------------------
// VerifyShiftAndAddMultiplicationModuliwise
//
// Operation to check whether a shift-and-add multiplier produces the correct result
//
// Iterate through moduli, and check whether a shift-and-add multiplier produces
// the correct result for all moduli
// -----------------------------------------------------------------------------
extern bool VerifyShiftAndAddMultiplicationModuliwise();

// -----------------------------------------------------------------------------
// SubtractiveKaratsubaOneLevel
//
// Create the CFLOBDD for the multiplication relation mod k by a one-level
// subtractive Karatsuba multiplier produces
// -----------------------------------------------------------------------------
extern CFLOBDD SubtractiveKaratsubaOneLevel(unsigned int k);

// -----------------------------------------------------------------------------
// VerifySubtractiveKaratsubaOneLevel
//
// Check whether one level of a subtractive Karatsuba multiplier produces 
// the correct result for input modulus k
// -----------------------------------------------------------------------------
extern bool VerifySubtractiveKaratsubaOneLevel(unsigned int k);

// -----------------------------------------------------------------------------
// VerifySubtractiveKaratsubaOneLevelModuliwise
//
// Iterate through moduli, and check whether a one-level Karatsuba multiplier produces
// the correct result for all moduli
// -----------------------------------------------------------------------------
extern bool VerifySubtractiveKaratsubaOneLevelModuliwise();

// Multiply two m-bit values using shift-and-add multiplication
// Return the 2*m-bit product
extern OUTPUT_TYPE ShiftAndAddMultiplication(INPUT_TYPE a, INPUT_TYPE b);

} // namespace CFL_OBDD

#endif
