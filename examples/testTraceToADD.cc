/**
  @file testTraceToADD.cc

  @brief Test driver for direct trace -> ADD -> CFLOBDD conversion.

  Reads a file as a stream of byte-valued tokens, builds an ADD directly
  using the binary-counter merge algorithm (no SEQUITUR), and converts
  to CFLOBDD.

  Usage: testTraceToADD <filename>
*/

#include "traceToAdd.h"
#include "addToCflobdd.hh"
#include "../CFLOBDD/cflobdd_int.h"
#include "../CFLOBDD/cflobdd_top_node_int.h"

#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <vector>

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
// Enumerate ADD decision tree leaves
// =========================================================================
static void enumerateADD(DdNode *node, unsigned int numVars, unsigned int depth,
                         std::vector<int> &leaves)
{
    if (Cudd_IsConstant(node)) {
        unsigned int count = 1u << (numVars - depth);
        int val = static_cast<int>(Cudd_V(node));
        for (unsigned int i = 0; i < count; i++)
            leaves.push_back(val);
        return;
    }

    unsigned int varIdx = Cudd_NodeReadIndex(node);

    if (varIdx > depth) {
        enumerateADD(node, numVars, depth + 1, leaves);
        enumerateADD(node, numVars, depth + 1, leaves);
    } else {
        enumerateADD(Cudd_E(node), numVars, depth + 1, leaves);
        enumerateADD(Cudd_T(node), numVars, depth + 1, leaves);
    }
}

// =========================================================================
// Verify ADD against file contents
// =========================================================================
static bool verifyADD(const std::string &filename,
                      const TraceADDResult &result)
{
    if (result.numVars > 25) {
        std::cout << "  (skipping verification: too many variables)" << std::endl;
        return true;
    }

    // Read expected values from file
    std::ifstream input(filename, std::ios::binary);
    std::vector<int> expected;
    char ch;
    while (input.get(ch)) {
        expected.push_back(static_cast<int>(static_cast<unsigned char>(ch)));
    }

    // Enumerate ADD decision tree
    std::vector<int> allLeaves;
    enumerateADD(result.add.getNode(), result.numVars, 0, allLeaves);

    // Squeeze out -1 entries
    std::vector<int> yield;
    for (int v : allLeaves) {
        if (v != -1) yield.push_back(v);
    }

    std::cout << "  ADD decision tree leaves: " << allLeaves.size()
              << " (2^" << result.numVars << ")" << std::endl;
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
        std::cerr << "Usage: testTraceToADD <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    // ---- Phase 1: Build ADD directly from trace ----
    Cudd mgr(0, 0, 256, 262144, 0);

    std::cout << "Building ADD directly from: " << filename << std::endl;
    auto t0 = high_resolution_clock::now();
    TraceADDResult result = traceFileToADD(mgr, filename);
    auto t1 = high_resolution_clock::now();
    double traceTime = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

    std::cout << "  Trace length: " << result.traceLength << " tokens" << std::endl;
    std::cout << "  Variables: " << result.numVars << std::endl;
    std::cout << "  ADD nodes: " << result.add.nodeCount() << std::endl;
    std::cout << "  Trace->ADD time: " << traceTime << "s" << std::endl;

    // ---- Phase 2: Verify ----
    std::cout << "\nVerifying ADD against file contents..." << std::endl;
    verifyADD(filename, result);

    // ---- Phase 3: Convert ADD to CFLOBDD ----
    // numVars is already a power of 2 (by construction in traceToADD)
    if (result.numVars > 0 &&
        result.numVars <= (1u << CFL_OBDD::CFLOBDDMaxLevel)) {
        InitCFLOBDD();

        auto t2 = high_resolution_clock::now();
        CFLOBDD cf = ADD_to_CFLOBDD<int>(mgr, result.add);
        auto t3 = high_resolution_clock::now();
        double cfTime = duration_cast<milliseconds>(t3 - t2).count() / 1000.0;

        unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
        cf.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);

        std::cout << "\n=== CFLOBDD ===" << std::endl;
        std::cout << "  ADD->CFLOBDD time: " << cfTime << "s" << std::endl;
        std::cout << "  Nodes: " << nodeCount << std::endl;
        std::cout << "  Edges: " << edgeCount << std::endl;
        std::cout << "  Return map entries: " << returnEdgeCount << std::endl;
        std::cout << "  Total: " << (nodeCount + edgeCount) << std::endl;

        // ---- Phase 4: Round-trip ----
        auto t4 = high_resolution_clock::now();
        ADD roundTrip = CFLOBDD_to_ADD<int>(mgr, cf);
        auto t5 = high_resolution_clock::now();
        double rtTime = duration_cast<milliseconds>(t5 - t4).count() / 1000.0;

        bool match = (roundTrip == result.add);
        std::cout << "\n=== Round-trip check ===" << std::endl;
        std::cout << "  CFLOBDD->ADD time: " << rtTime << "s" << std::endl;
        std::cout << "  ADD == original ADD: " << (match ? "PASSED" : "FAILED") << std::endl;
    } else {
        std::cout << "\n(Skipping CFLOBDD conversion)" << std::endl;
    }

    return 0;
}
