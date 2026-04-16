/**
  @file traceToAdd.cc

  @brief Implementation of direct trace-to-ADD construction.
*/

#include "traceToAdd.h"

#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>

// Maximum number of Boolean variables (supports traces up to 2^63 tokens)
static const unsigned int MAX_TRACE_VARS = 63;

// =========================================================================
// CharTokenizer: each byte is a token with value 0..255
// =========================================================================

class CharTokenizer : public Tokenizer {
public:
    bool open(const std::string &filename) override {
        input.open(filename, std::ios::binary);
        return input.is_open();
    }

    bool next(unsigned int *value) override {
        char ch;
        if (input.get(ch)) {
            *value = static_cast<unsigned int>(static_cast<unsigned char>(ch));
            return true;
        }
        return false;
    }

    unsigned long count() override {
        return 0;  // unknown length for streaming
    }

    void close() override {
        input.close();
    }

private:
    std::ifstream input;
};

// =========================================================================
// traceToADD: streaming binary-counter merge algorithm
//
// The algorithm is analogous to incrementing a binary counter.
// We maintain a vector V of partial ADDs at increasing depths:
//
//   V[0] holds a depth-0 ADD (a single constant), or is empty
//   V[1] holds a depth-1 ADD (one variable, two leaves), or is empty
//   ...
//   V[d] holds the deepest ADD built so far
//
// The companion vector "occupied" tracks which slots of V hold valid
// ADDs.  (V entries are never removed; "occupied" controls validity.)
//
// insertToken works like adding 1 to a binary counter:
//   - The new constant ADD is the incoming "1"
//   - If V[0] is empty, store there (like 0+1=1, no carry)
//   - If V[0] is occupied, combine V[0] (left child) with the new ADD
//     (right child) using an ITE on the appropriate Boolean variable,
//     producing a "carry" ADD of depth 1.  Clear V[0].
//   - Propagate the carry upward: if V[1] is occupied, combine and
//     carry again; if V[1] is empty, store and stop.
//   - If the carry propagates past all existing slots, grow V by one
//     (allocating a new Boolean variable).
//
// Variables are allocated (bottom-up) by counting down from MAX_TRACE_VARS: depth-i
// combines use CUDD variable (MAX_TRACE_VARS - i).  This puts the
// innermost (LSB) variables at the highest CUDD indices and the
// outermost (MSB) at the lowest, matching CUDD's convention that
// lower-indexed variables are closer to the root.
//
// After the build, variables are shifted to use indices 0..h-1.
// If h is not a power of 2, additional -1 padding rounds it up
// so the result is compatible with ADD_to_CFLOBDD.
// =========================================================================

TraceADDResult traceToADD(Cudd &mgr, Tokenizer &tok)
{
    // V[i]: partial ADD of depth i (valid only when occupied[i] is true).
    // Entries are never removed from V; occupied[i] controls validity.
    std::vector<ADD> V;
    std::vector<bool> occupied;
    ADD dontCare = mgr.constant(-1);

    // numPlies: number of Boolean variables allocated so far.
    unsigned int numPlies = 0;

    auto insertToken = [&](ADD token) {
        // "carry" is the ADD being propagated upward, analogous to the
        // carry bit in binary addition.  It starts as the depth-0
        // constant for the new token.
        ADD carry = token;
        unsigned int ply = 0;
        while (ply < V.size() && occupied[ply]) {
            // Combine V[ply] (left child, earlier in trace) with carry
            // (right child, later in trace).  The Boolean variable for
            // ply i is (MAX_TRACE_VARS - i): high indices at the bottom,
            // low indices at the top.
            unsigned int varIdx = MAX_TRACE_VARS - ply;
            carry = mgr.addVar(varIdx).Ite(carry, V[ply]);
            occupied[ply] = false;
            ply++;
        }
        if (ply == V.size()) {
            // Carry propagated past all existing slots: allocate a new
            // ply to hold the result.  A new Boolean variable was used
            // during the combine at ply (V.size() - 1), so set
            // numPlies to (V.size() - 1).
            assert(V.size() <= MAX_TRACE_VARS && "Exceeded maximum trace variables");
            V.push_back(carry);
            occupied.push_back(true);
            numPlies = V.size() - 1;
        } else {
            // Found an empty slot: store the carry.
            V[ply] = carry;
            occupied[ply] = true;
        }
    };

    // Read tokens and build ADD incrementally
    unsigned long long tokenCount = 0;
    unsigned int val;
    while (tok.next(&val)) {
        insertToken(mgr.constant(static_cast<double>(val)));
        tokenCount++;
    }

    if (tokenCount == 0) {
        return {mgr.constant(-1), 0, 0};
    }

    // Pad with -1 until everything collapses into a single ADD.
    // This rounds the token count up to the next power of 2.
    while (true) {
        unsigned int numOccupied = 0;
        unsigned int topOccupied = 0;
        for (unsigned int i = 0; i < V.size(); i++) {
            if (occupied[i]) {
                topOccupied = i;
                numOccupied++;
            }
        }
        if (numOccupied == 1 && topOccupied == V.size() - 1)
            break;
        insertToken(dontCare);
    }

    unsigned int h = numPlies;  // number of Boolean variables used
    // The complete ADD is in the topmost occupied slot of V.
    unsigned int topSlot = V.size() - 1;
    ADD topADD = V[topSlot];

    // If h is not a power of 2, pad to the next power of 2 by adding
    // new variables on top.  Each new variable selects the existing ADD
    // on the left (var=0) and -1 on the right (var=1), doubling the
    // decision tree height by one ply each time.
    unsigned int targetH = 1;
    while (targetH < h) targetH <<= 1;
    for (unsigned int i = h; i < targetH; i++) {
        unsigned int varIdx = MAX_TRACE_VARS - i;
        topADD = mgr.addVar(varIdx).Ite(dontCare, topADD);
    }
    h = targetH;

    // Shift variables from [MAX_TRACE_VARS-h+1 .. MAX_TRACE_VARS] to [0 .. h-1].
    // Ply i used variable (MAX_TRACE_VARS - i); we want it at (h - 1 - i).
    // This preserves the ordering: the outermost variable (highest ply)
    // maps to index 0 (CUDD root), the innermost to index h-1 (CUDD bottom).
    std::vector<ADD> from(h), to(h);
    for (unsigned int i = 0; i < h; i++) {
        from[i] = mgr.addVar(MAX_TRACE_VARS - i);
        to[i] = mgr.addVar(h - 1 - i);
    }
    ADD result = topADD.SwapVariables(from, to);

    TraceADDResult ret;
    ret.add = result;
    ret.numVars = h;
    ret.traceLength = static_cast<unsigned long>(tokenCount);
    return ret;
}

// =========================================================================
// traceFileToADD: convenience wrapper
// =========================================================================

TraceADDResult traceFileToADD(Cudd &mgr, const std::string &filename)
{
    CharTokenizer tok;
    if (!tok.open(filename)) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return {mgr.constant(-1), 0, 0};
    }
    TraceADDResult result = traceToADD(mgr, tok);
    tok.close();
    return result;
}
