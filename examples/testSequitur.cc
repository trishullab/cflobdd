/**
  @file testSequitur.cc

  @brief Test driver for SEQUITUR grammar -> ADD -> CFLOBDD conversion.

  Reads a file, builds a SEQUITUR grammar over int-valued symbols,
  converts the grammar to a CUDD ADD, and optionally converts to CFLOBDD.

  Usage: testSequitur <filename>

  The file is read character-by-character; each character's ASCII value
  is used as the int terminal value.
*/

#include "sequiturToAdd.hh"
#include "addToCflobdd.hh"
#include "../CFLOBDD/cflobdd_int.h"
#include "../CFLOBDD/cflobdd_top_node_int.h"

#include <fstream>
#include <iostream>
#include <string>
#include <chrono>

using namespace CFL_OBDD;
using namespace std::chrono;

// =========================================================================
// CFLOBDD subsystem initialization
// =========================================================================
static bool cflobddInitialized = false;

static void InitCFLOBDD() {
    if (cflobddInitialized) return;
    CFLOBDDNodeHandle::InitNoDistinctionTable();
    CFLOBDDNodeHandle::InitReduceCache();
    InitPairProductCache();
    InitTripleProductCache();
    cflobddInitialized = true;
}

// =========================================================================
// Enumerate all paths through the ADD (decision tree traversal), collect
// leaf values, squeeze out -1 entries.  The result should equal the
// original string.
// =========================================================================
static void enumerateADD(DdNode *node, int numVars, int depth,
                         std::vector<int> &leaves)
{
    if (Cudd_IsConstant(node)) {
        // At a leaf: this leaf covers 2^(numVars - depth) assignments
        int count = 1 << (numVars - depth);
        int val = static_cast<int>(Cudd_V(node));
        for (int i = 0; i < count; i++)
            leaves.push_back(val);
        return;
    }

    int varIdx = Cudd_NodeReadIndex(node);

    if (varIdx > depth) {
        // Skipped variables: the node doesn't depend on variable 'depth'.
        // Both values of variable 'depth' lead to the same subtree (this node).
        // Recurse twice with depth+1, keeping the same node.
        enumerateADD(node, numVars, depth + 1, leaves);  // depth-variable = 0
        enumerateADD(node, numVars, depth + 1, leaves);  // depth-variable = 1
    } else {
        // No skip: this node tests variable 'depth'
        enumerateADD(Cudd_E(node), numVars, depth + 1, leaves);  // variable = 0
        enumerateADD(Cudd_T(node), numVars, depth + 1, leaves);  // variable = 1
    }
}

static bool verifyADD(Cudd &mgr, const Sequitur<int> &seq,
                      const GrammarADDResult &result)
{
    int numVars = result.maxVar + 1;
    if (numVars > 25) {
        std::cout << "  (skipping verification: too many variables)" << std::endl;
        return true;
    }

    // Enumerate all leaf values from the ADD decision tree
    std::vector<int> allLeaves;
    enumerateADD(result.add.getNode(), numVars, 0, allLeaves);

    // Squeeze out -1 entries
    std::vector<int> yield;
    for (int v : allLeaves) {
        if (v != -1) yield.push_back(v);
    }

    // Get expected string from SEQUITUR iteration
    std::vector<int> expected;
    for (auto it = seq.begin(); it != seq.end(); ++it) {
        expected.push_back(*it);
    }

    std::cout << "  ADD decision tree leaves: " << allLeaves.size()
              << " (2^" << numVars << " = " << (1 << numVars) << ")" << std::endl;
    std::cout << "  Non-(-1) leaves: " << yield.size() << std::endl;
    std::cout << "  Expected length: " << expected.size() << std::endl;

    if (yield.size() != expected.size()) {
        std::cout << "  Verification FAILED: length mismatch" << std::endl;
        return false;
    }

    int errors = 0;
    for (size_t i = 0; i < expected.size(); i++) {
        if (yield[i] != expected[i]) {
            if (errors < 10) {
                std::cout << "  MISMATCH at position " << i
                          << ": expected " << expected[i]
                          << " ('" << (char)expected[i] << "')"
                          << ", got " << yield[i]
                          << " ('" << (char)yield[i] << "')" << std::endl;
            }
            errors++;
        }
    }

    if (errors == 0) {
        std::cout << "  Verification PASSED" << std::endl;
    } else {
        std::cout << "  Verification FAILED: " << errors << " mismatches" << std::endl;
    }
    return errors == 0;
}

// =========================================================================
// Main
// =========================================================================
int main(int argc, char *argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: testSequitur <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::ifstream input(filename, std::ios::binary);
    if (!input.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return 1;
    }

    // ---- Phase 1: Build SEQUITUR grammar ----
    std::cout << "Building SEQUITUR grammar from: " << filename << std::endl;
    auto t0 = high_resolution_clock::now();

    Sequitur<int> seq;
    char ch;
    unsigned long count = 0;
    while (input.get(ch)) {
        seq.push_back(static_cast<int>(ch));
        count++;
    }

    auto t1 = high_resolution_clock::now();
    double seqTime = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;
    std::cout << "  Input length: " << count << " symbols" << std::endl;
    std::cout << "  SEQUITUR time: " << seqTime << "s" << std::endl;

    // Print grammar
    std::cout << "\nGrammar rules:" << std::endl;
    seq.printRules();
    std::cout << std::endl;

    // ---- Phase 2: Convert grammar to ADD ----
    Cudd mgr(0, 0, 256, 262144, 0);

    auto t2 = high_resolution_clock::now();
    GrammarADDResult result = sequiturToADD(mgr, seq);
    auto t3 = high_resolution_clock::now();
    double addTime = duration_cast<milliseconds>(t3 - t2).count() / 1000.0;

    printSequiturADDStats(seq, result);
    std::cout << "  Grammar->ADD time: " << addTime << "s" << std::endl;

    // ---- Phase 3: Verify ----
    std::cout << "\nVerifying ADD against SEQUITUR iteration..." << std::endl;
    verifyADD(mgr, seq, result);

    // ---- Phase 4: Convert ADD to CFLOBDD ----
    // Pad ADD to power-of-2 variables: add variables on top, left child = existing ADD, right = -1
    int nv = result.maxVar + 1;
    ADD paddedADD = result.add;
    while (nv > 0 && (nv & (nv - 1)) != 0) {
        paddedADD = mgr.addVar(nv).Ite(mgr.constant(-1), paddedADD);
        nv++;
    }
    if (nv == 0) nv = 1;  // degenerate: constant ADD needs at least 1 variable
    std::cout << "\n  Padded to " << nv << " variables (2^"
              << static_cast<int>(std::log2(nv)) << ")" << std::endl;

    if (nv <= (1 << CFL_OBDD::CFLOBDDMaxLevel)) {
        InitCFLOBDD();

        auto t4 = high_resolution_clock::now();
        CFLOBDD cf = ADD_to_CFLOBDD<int>(mgr, paddedADD);
        auto t5 = high_resolution_clock::now();
        double cfTime = duration_cast<milliseconds>(t5 - t4).count() / 1000.0;

        unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
        cf.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);

        std::cout << "\n=== CFLOBDD ===" << std::endl;
        std::cout << "  ADD->CFLOBDD time: " << cfTime << "s" << std::endl;
        std::cout << "  Nodes: " << nodeCount << std::endl;
        std::cout << "  Edges: " << edgeCount << std::endl;
        std::cout << "  Return map entries: " << returnEdgeCount << std::endl;
        std::cout << "  Total: " << (nodeCount + edgeCount) << std::endl;
    } else {
        std::cout << "\n(Skipping CFLOBDD conversion: too many variables for CFLOBDDMaxLevel)"
                  << std::endl;
    }

    return 0;
}
