/**
  @file addMultiplication.hh

  @brief ADD-based multiplication modulo k using the CRT approach.

  @details This file provides functions to build ADDs (Algebraic Decision
  Diagrams) that represent:
  - NumsModK: λx ∈ B^n . val(x) mod k
  - MultModK: λ(x,y) ∈ B^n × B^n . (val(x) × val(y)) mod k

  These are the ADD equivalents of the CFLOBDD functions described in
  "Multiplication_via_CRT.pdf" (Figures 8-10).

  The key insight from Figure 2 of the document is that the MTBDD/ADD for
  NumsModK acts like a finite-state automaton:
  - States represent value mod k of the prefix read so far
  - Going left (0) multiplies by 2 mod k
  - Going right (1) multiplies by 2 and adds 1, mod k
  - There are at most k nodes per ply
  - Total size is O(k × n) where n is the number of bits

  For MultModK, we:
  1. Build an ADD for NumsModK over the first n bits (for input x)
  2. Build an ADD for NumsModK over the second n bits (for input y)
  3. Combine them with multiplication mod k at the leaves

  The variables are non-interleaved: all x bits first (variables 0..n-1),
  then all y bits (variables n..2n-1).

  @author Based on "Multiplication_via_CRT.pdf" by Thomas Reps

*/

#ifndef ADD_MULTIPLICATION_HH_
#define ADD_MULTIPLICATION_HH_

#include "../cudd-3.0.0/cplusplus/cuddObj.hh"
#include <boost/multiprecision/cpp_int.hpp>

// =============================================================================
// Bit-width configuration
// =============================================================================

#ifndef NUM_BITS
#define NUM_BITS 64
#endif

// =============================================================================
// Constants for CRT-based multiplication (Figure 12 in Multiplication_via_CRT.pdf)
// =============================================================================

/// Number of primes needed per bit-width (product must exceed 2^(2*bitwidth))
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

/// Total number of odd primes available
const unsigned int numberOfOddPrimes = 999;

/// Array of first 999 odd primes
extern const unsigned int AllModuli[numberOfOddPrimes];

/// Pointer to the moduli array (for backward compatibility)
extern const unsigned int* Moduli;

// =============================================================================
// Position enumeration for Karatsuba splitting
// =============================================================================

/// Position enumeration for variable placement in ADDs
/// TopLevel: not used for NumsModK
/// A, B: standard positions for two input operands (each n bits)
/// AA, AB, BA, BB: half-width positions for Karatsuba decomposition
enum Position { TopLevel, A, B, AA, AB, BA, BB };

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


// =============================================================================
// MultRelation class
// =============================================================================

/**
  @brief Class to hold ADDs for modular multiplication relations.

  @details This class stores ADDs for multiplication mod k for each of the
  first 26 odd primes. The product of these primes exceeds 2^133, which is
  sufficient to uniquely determine the product of two 64-bit numbers via
  the Chinese Remainder Theorem.

  This is the ADD equivalent of MultRelation in the CFLOBDD implementation
  (Figure 12 in the document).

*/
class MultRelation {
public:
    /// Number of bits per operand (configured by NUM_BITS)
    static const int numBits = NUM_BITS;

    /// Array of moduli (first 26 odd primes)
    unsigned int ModuliArray[numberOfMultRelations];

    /// ADDs for multiplication mod each prime
    ADD ModularMultRelations[numberOfMultRelations];

    /// CUDD manager (shared by all ADDs)
    Cudd* mgr;

    /**
      @brief Constructor that builds ADDs for all 26 primes.

      @param manager  CUDD manager to use for building ADDs
    */
    MultRelation(Cudd& manager);

    /**
      @brief Default constructor (creates its own manager).
    */
    MultRelation();

    /**
      @brief Destructor.
    */
    ~MultRelation();

    /**
      @brief Compares two MultRelation objects for equality.
    */
    bool operator==(const MultRelation& other) const;

    /**
      @brief Multiplies two values using CRT and Garner's algorithm.

      @details Evaluates each ADD to get remainders, then reconstructs the
      double-width result using Garner's algorithm.

      @param a  First operand (INPUT_TYPE)
      @param b  Second operand (INPUT_TYPE)

      @return The product (OUTPUT_TYPE)
    */
    OUTPUT_TYPE Multiply(INPUT_TYPE a, INPUT_TYPE b) const;

    /**
      @brief Prints statistics about all ADDs in the relation.
    */
    void PrintStats() const;
};

/**
  @brief Builds an ADD that computes val(x) mod k for an n-bit input.

  @details This function creates an ADD over variables top..top+n-1 that
  computes the value of an n-bit unsigned binary number modulo k.

  The bit order is high-order to low-order (MSB at variable 'top').

  This is equivalent to NumsModK in the CFLOBDD implementation (Figure 8),
  and is a wrapper around Cudd::addResidue with CUDD_RESIDUE_MSB option.

  @param mgr      CUDD manager
  @param n        Number of bits
  @param k        Modulus (must be >= 2)
  @param top      Index of the top (most significant) variable

  @return ADD representing val(x) mod k

  @see ADD_MultModK

*/
ADD ADD_NumsModK(Cudd& mgr, int n, int k, int top);


/**
  @brief Builds an ADD that computes (val(x) × val(y)) mod k.

  @details This function creates an ADD over 2n variables that computes
  the product of two n-bit unsigned binary numbers modulo k.

  Variables 0..n-1 represent x (the first operand).
  Variables n..2n-1 represent y (the second operand).

  The bit order is high-order to low-order within each operand.

  This is equivalent to MultModK in the CFLOBDD implementation (Figure 10).

  @param mgr      CUDD manager
  @param n        Number of bits for each operand
  @param k        Modulus (must be >= 2)

  @return ADD representing (val(x) × val(y)) mod k

  @see ADD_NumsModK

*/
ADD ADD_MultModK(Cudd& mgr, int n, int k);


/**
  @brief Evaluates an ADD for a given 2n-bit input (x, y).

  @details Given two n-bit unsigned numbers x and y, evaluates the ADD
  and returns the leaf value.

  The bit order assumed is high-order to low-order (MSB at lower variable
  indices).

  @param mgr      CUDD manager
  @param f        ADD to evaluate
  @param n        Number of bits per operand
  @param x        First n-bit operand value (INPUT_TYPE)
  @param y        Second n-bit operand value (INPUT_TYPE)

  @return The leaf value (double) for the given input.

*/
CUDD_VALUE_TYPE ADD_Evaluate(Cudd& mgr, const ADD& f, int n,
                              INPUT_TYPE x, INPUT_TYPE y);


/**
  @brief Verifies that ADD_MultModK computes the correct result.

  @details Tests the ADD_MultModK function by comparing its output against
  direct computation of (x * y) mod k for a range of inputs.

  @param mgr      CUDD manager
  @param n        Number of bits for each operand
  @param k        Modulus

  @return true if all tests pass, false otherwise.

*/
bool verifyMultModK(Cudd& mgr, int n, int k);


/**
  @brief Prints statistics about an ADD.

  @details Prints the node count and leaf count of the ADD.

  @param f        ADD to analyze
  @param name     Name to display

*/
void printADDStats(const ADD& f, const char* name);


// =============================================================================
// Shift-and-Add Multiplication
// =============================================================================

/**
  @brief Builds multiplication mod k using the shift-and-add algorithm.

  @details Algorithm:
  result = 0
  for i = 0 to n-1 (high-order to low-order):
      p = projection of bit i of first input
      addend = p * b mod k
      result = result * 2 mod k
      result = result + addend mod k

  This mirrors ShiftAndAddMultiplicationModK in multiplication_crt.cpp.

  @param mgr  CUDD manager
  @param n    Number of bits per operand
  @param k    Modulus

  @return ADD representing shift-and-add multiplication mod k
*/
ADD ADD_ShiftAndAddMultiplicationModK(Cudd& mgr, int n, int k);

/**
  @brief Verifies shift-and-add multiplication for a single modulus.

  @param mgr  CUDD manager
  @param n    Number of bits per operand
  @param k    Modulus

  @return true if verification passes
*/
bool ADD_VerifyShiftAndAddMultiplicationModK(Cudd& mgr, int n, int k);

/**
  @brief Builds specification ADDs for all moduli with aggregate and final-modulus timing.

  @param mgr  CUDD manager
  @param n    Number of bits per operand
*/
void ADD_BuildMultiplicationSpecsModuliwise(Cudd& mgr, int n);

/**
  @brief Verifies shift-and-add multiplication for all moduli.

  @param mgr  CUDD manager
  @param n    Number of bits per operand

  @return true if all verifications pass
*/
bool ADD_VerifyShiftAndAddMultiplicationModuliwise(Cudd& mgr, int n);


// =============================================================================
// Karatsuba Multiplication
// =============================================================================

/**
  @brief Builds one level of subtractive Karatsuba multiplication mod k.

  @details Algorithm:
  Split X = X1*2^m + X0, Y = Y1*2^m + Y0 where m = n/2
  Z0 = X0 * Y0 mod k
  Z2 = X1 * Y1 mod k
  Z1 = (X1 - X0) * (Y0 - Y1) mod k
  result = Z2*2^(2m) + (Z1 + Z2 + Z0)*2^m + Z0 (all mod k)

  This mirrors SubtractiveKaratsubaOneLevel in multiplication_crt.cpp.

  @param mgr  CUDD manager
  @param n    Number of bits per operand (must be even)
  @param k    Modulus

  @return ADD representing Karatsuba multiplication mod k
*/
ADD ADD_SubtractiveKaratsubaOneLevel(Cudd& mgr, int n, int k);

/**
  @brief Verifies Karatsuba multiplication for a single modulus.

  @param mgr  CUDD manager
  @param n    Number of bits per operand
  @param k    Modulus

  @return true if verification passes
*/
bool ADD_VerifySubtractiveKaratsubaOneLevel(Cudd& mgr, int n, int k);

/**
  @brief Verifies Karatsuba multiplication for all moduli.

  @param mgr  CUDD manager
  @param n    Number of bits per operand

  @return true if all verifications pass
*/
bool ADD_VerifySubtractiveKaratsubaOneLevelModuliwise(Cudd& mgr, int n);


// =============================================================================
// Helper functions for modular arithmetic
// =============================================================================

/**
  @brief Computes val(x) mod k for a positioned input.

  @details Maps Position enum to variable ranges:
  - A: variables 0..n-1 (first input)
  - B: variables n..2n-1 (second input)
  - AA: variables 0..n/2-1 (high half of first input)
  - AB: variables n/2..n-1 (low half of first input)
  - BA: variables n..n+n/2-1 (high half of second input)
  - BB: variables n+n/2..2n-1 (low half of second input)

  @param mgr  CUDD manager
  @param n    Total number of bits per operand
  @param k    Modulus
  @param p    Position specifier

  @return ADD representing val(x) mod k for the specified position
*/
ADD ADD_NumsModKPositioned(Cudd& mgr, int n, int k, Position p);

/**
  @brief Adds two ADDs element-wise modulo k.

  @param mgr  CUDD manager
  @param a    First ADD operand
  @param b    Second ADD operand
  @param k    Modulus

  @return ADD representing (a + b) mod k
*/
ADD ADD_AddModK(Cudd& mgr, const ADD& a, const ADD& b, int k);

/**
  @brief Subtracts two ADDs element-wise modulo k.

  @param mgr  CUDD manager
  @param a    First ADD operand
  @param b    Second ADD operand
  @param k    Modulus

  @return ADD representing (a - b + k) mod k
*/
ADD ADD_SubtractModK(Cudd& mgr, const ADD& a, const ADD& b, int k);

/**
  @brief Multiplies two ADDs element-wise modulo k.

  @param mgr  CUDD manager
  @param a    First ADD operand
  @param b    Second ADD operand
  @param k    Modulus

  @return ADD representing (a * b) mod k
*/
ADD ADD_MultiplyModK(Cudd& mgr, const ADD& a, const ADD& b, int k);


#endif /* ADD_MULTIPLICATION_HH_ */
