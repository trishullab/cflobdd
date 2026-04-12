/**
  @file sequiturToAdd.h

  @brief Convert a SEQUITUR grammar (Sequitur<int>) to a CUDD ADD.

  Uses the non-dense selector-variable approach:
  - Each rule P -> Q1 Q2 ... Qk allocates ceil(log2(k)) selector variables
  - The selector encodes which child Qi; indices >= k map to constant -1
  - Each nonterminal's ADD is built once and memoized
  - Repeated children share the same sub-ADD (canonicalized by CUDD)

  The resulting ADD maps Boolean variable assignments to int terminal values,
  representing the string yielded by the grammar as a function of position.
*/

#ifndef SEQUITUR_TO_ADD_H_
#define SEQUITUR_TO_ADD_H_

#include "../cudd-3.0.0/cplusplus/cuddObj.hh"
#include "../sequitur/sequitur.hpp"

using namespace jw;

// Result of converting a grammar node to an ADD
struct GrammarADDResult {
    ADD add;           // The ADD for this node's yield
    int maxVar;        // Maximum variable index used (-1 for constants)
    unsigned long yieldLength;  // Length of this node's yield
};

/// Convert a SEQUITUR grammar to a CUDD ADD.
/// Walks the grammar DAG bottom-up (memoized per rule ID), building an ADD
/// for each nonterminal using the non-dense selector-variable approach.
GrammarADDResult sequiturToADD(Cudd &mgr, const Sequitur<int> &seq);

/// Print statistics about the grammar and resulting ADD.
void printSequiturADDStats(const Sequitur<int> &seq,
                           const GrammarADDResult &result);

#endif // SEQUITUR_TO_ADD_H_
