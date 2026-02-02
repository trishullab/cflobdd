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

// =============================================================================
// Constants for CRT-based multiplication (Figure 12 in Multiplication_via_CRT.pdf)
// =============================================================================

/// Number of primes used for CRT-based multiplication verification
const unsigned int numberOfMultRelations = 26;

/// First 26 odd primes (product > 2^133, sufficient for 64-bit multiplication)
extern const unsigned int Moduli[numberOfMultRelations];


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
    /// Number of bits per operand (64 for standard multiplication)
    static const int numBits = 64;

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
      @brief Multiplies two 64-bit values using CRT and Garner's algorithm.

      @details Evaluates each ADD to get remainders, then reconstructs the
      128-bit result using Garner's algorithm.

      @param a  First 64-bit operand
      @param b  Second 64-bit operand

      @return The 128-bit product (as two 64-bit values: high and low)
    */
    void Multiply(unsigned long long a, unsigned long long b,
                  unsigned long long& resultHigh, unsigned long long& resultLow) const;

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
  @param x        First n-bit operand value
  @param y        Second n-bit operand value

  @return The leaf value (double) for the given input.

*/
CUDD_VALUE_TYPE ADD_Evaluate(Cudd& mgr, const ADD& f, int n,
                              unsigned long long x, unsigned long long y);


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


#endif /* ADD_MULTIPLICATION_HH_ */
