/**
  @file testConversion.cc

  @brief Test driver for bidirectional ADD <-> CFLOBDD conversion.

  Tests use structural equality (pointer equality) first.  Exhaustive
  evaluation on all assignments is performed only when structural equality
  fails or is not applicable.  This allows tests with many more variables
  than would be feasible with brute-force enumeration.

  CRT tests (NumsModK, MultModK) require -DNUM_BITS=32.
  Build with: make testConversion (from examples/ directory)
  Requires lld linker: pacman -S mingw-w64-ucrt-x86_64-lld
*/

#include "addToCflobdd.hh"
#include "../CFLOBDD/cflobdd_int.h"
#include "../CFLOBDD/cflobdd_top_node_int.h"

#ifdef NUM_BITS
#include "../CFLOBDD/multiplication_crt.h"
#define USE_CRT_TESTS
#endif

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <chrono>

using namespace CFL_OBDD;
using namespace std::chrono;

// Provide CFLTests::verbose (defined in tests_cfl.cpp, which we don't link)
#ifdef USE_CRT_TESTS
#include "../CFLOBDD/tests_cfl.h"
bool CFLTests::verbose = false;
#endif

// =========================================================================
// CFLOBDD subsystem initialization (mirrors CFLOBDD_module_init)
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
// Evaluate a CUDD ADD on a given bit-vector assignment
// =========================================================================
static double EvalADD(const ADD &f, unsigned int numVars, unsigned int assignment) {
    int* inputs = new int[numVars];
    for (unsigned int i = 0; i < numVars; i++) {
        inputs[i] = (assignment >> i) & 1;
    }
    ADD result = f.Eval(inputs);
    delete[] inputs;
    return cuddV(Cudd_Regular(result.getNode()));
}

// =========================================================================
// Evaluate a CFLOBDD on a given bit-vector assignment
// =========================================================================
static int EvalCFLOBDD(const CFLOBDD &f, unsigned int numVars, unsigned int assignment) {
    SH_OBDD::Assignment a(numVars);
    for (unsigned int i = 0; i < numVars; i++) {
        a[i] = (assignment >> i) & 1;
    }
    return f.root->EvaluateIteratively(a);
}

// =========================================================================
// Exhaustive comparison of ADD and CFLOBDD (used only as fallback)
// =========================================================================
static bool ExhaustiveCompare(const ADD &add, const CFLOBDD &cf, unsigned int numVars,
                               const char* label) {
    unsigned int total = 1u << numVars;
    for (unsigned int a = 0; a < total; a++) {
        int cfVal = EvalCFLOBDD(cf, numVars, a);
        double addVal = EvalADD(add, numVars, a);
        if (cfVal != (int)addVal) {
            std::cout << "FAIL [" << label << "]: assignment " << a
                      << " CFLOBDD=" << cfVal << " ADD=" << addVal << std::endl;
            return false;
        }
    }
    return true;
}

// =========================================================================
// Test 1: Round-trip ADD -> CFLOBDD -> ADD (structural equality on ADDs)
// =========================================================================
static bool testRoundTrip_ADD() {
    std::cout << "=== Test 1: Round-trip ADD -> CFLOBDD -> ADD ===" << std::endl;

    // 2 variables
    {
        Cudd mgr(0, 0);
        ADD x0 = mgr.addVar(0), x1 = mgr.addVar(1);
        ADD fs[] = {
            mgr.constant(42.0),      // constant
            x0,                       // projection
            x0 + x1,                 // sum
            x0 * x1,                 // product (AND)
        };
        const char* names[] = { "const(42)", "x0", "x0+x1", "x0*x1" };

        for (int t = 0; t < 4; t++) {
            CFLOBDD cf = ADD_to_CFLOBDD<int>(mgr, fs[t]);
            ADD rt = CFLOBDD_to_ADD<int>(mgr, cf);
            if (fs[t] == rt) {
                std::cout << "  1." << t << " " << names[t] << " (2 vars): PASS (structural)" << std::endl;
            } else {
                std::cout << "  1." << t << " " << names[t] << ": structural FAIL, checking exhaustively..." << std::endl;
                if (!ExhaustiveCompare(fs[t], cf, 2, names[t])) return false;
                std::cout << "  1." << t << " " << names[t] << ": PASS (semantic)" << std::endl;
            }
        }
    }

    // 4 variables
    {
        Cudd mgr(0, 0);
        ADD v0 = mgr.addVar(0), v1 = mgr.addVar(1);
        ADD v2 = mgr.addVar(2), v3 = mgr.addVar(3);
        ADD f = v0 * mgr.constant(8.0) + v1 * mgr.constant(4.0) +
                v2 * mgr.constant(2.0) + v3;
        CFLOBDD cf = ADD_to_CFLOBDD<int>(mgr, f);
        ADD rt = CFLOBDD_to_ADD<int>(mgr, cf);
        if (f == rt) {
            std::cout << "  1.4 weighted-sum (4 vars): PASS (structural)" << std::endl;
        } else {
            std::cout << "  1.4 weighted-sum: structural FAIL, checking exhaustively..." << std::endl;
            if (!ExhaustiveCompare(f, cf, 4, "weighted-sum-4")) return false;
            std::cout << "  1.4 weighted-sum: PASS (semantic)" << std::endl;
        }
    }

    // 16 variables — structural only (2^16 exhaustive not needed if structural passes)
    {
        Cudd mgr(0, 0);
        unsigned int n = 16;
        for (unsigned int i = 0; i < n; i++) mgr.addVar(i);

        // Parity function via XOR arithmetic
        ADD f = mgr.constant(0.0);
        for (unsigned int i = 0; i < n; i++) {
            ADD xi = mgr.addVar(i);
            f = f + xi - mgr.constant(2.0) * f * xi;
        }
        CFLOBDD cf = ADD_to_CFLOBDD<int>(mgr, f);
        ADD rt = CFLOBDD_to_ADD<int>(mgr, cf);
        if (f == rt) {
            std::cout << "  1.5 parity (16 vars): PASS (structural)" << std::endl;
        } else {
            std::cout << "  1.5 parity (16 vars): structural FAIL, checking exhaustively..." << std::endl;
            if (!ExhaustiveCompare(f, cf, n, "parity-16")) return false;
            std::cout << "  1.5 parity: PASS (semantic)" << std::endl;
        }
    }

    return true;
}

// =========================================================================
// Test 2: Round-trip CFLOBDD -> ADD -> CFLOBDD (structural equality on CFLOBDDs)
// =========================================================================
static bool testRoundTrip_CFLOBDD() {
    std::cout << "=== Test 2: Round-trip CFLOBDD -> ADD -> CFLOBDD ===" << std::endl;

    // Projection functions at various levels
    for (unsigned int level = 1; level <= 4; level++) {
        unsigned int numVars = 1u << level;
        Cudd mgr(0, 0);
        for (unsigned int i = 0; i < numVars; i++) mgr.addVar(i);

        // Test projection on each variable
        bool allStructural = true;
        for (unsigned int v = 0; v < numVars; v++) {
            CFLOBDD orig = MkProjection(v, level);
            ADD add = CFLOBDD_to_ADD<int>(mgr, orig);
            CFLOBDD rt = ADD_to_CFLOBDD<int>(mgr, add);

            if (!(orig == rt)) {
                allStructural = false;
                std::cout << "  Structural FAIL for Projection(" << v
                          << ") at level " << level << ", checking exhaustively..." << std::endl;
                if (!ExhaustiveCompare(add, orig, numVars, "projection"))
                    return false;
            }
        }
        if (allStructural) {
            std::cout << "  2." << level << " All " << numVars << " projections (level "
                      << level << "): PASS (structural)" << std::endl;
        }
    }

    // Parity at level 4 (16 vars) via XOR of projections
    {
        unsigned int level = 4;
        unsigned int numVars = 1u << level;
        Cudd mgr(0, 0);
        for (unsigned int i = 0; i < numVars; i++) mgr.addVar(i);

        CFLOBDD orig = MkProjection(0, level);
        for (unsigned int i = 1; i < numVars; i++) {
            orig = MkExclusiveOr(orig, MkProjection(i, level));
        }
        ADD add = CFLOBDD_to_ADD<int>(mgr, orig);
        CFLOBDD rt = ADD_to_CFLOBDD<int>(mgr, add);
        if (orig == rt) {
            std::cout << "  2.5 Parity-XOR (level " << level << ", "
                      << numVars << " vars): PASS (structural)" << std::endl;
        } else {
            std::cout << "  2.5 Parity-XOR (level " << level << "): structural FAIL" << std::endl;
            return false;
        }
    }

    return true;
}

#ifdef USE_CRT_TESTS
// =========================================================================
// Test 3: NumsModK round-trip (CRT primitive)
//
// Build NumsModK as a CFLOBDD, convert to ADD and back, check structural
// equality on the CFLOBDD side.
// =========================================================================
static bool testNumsModK() {
    std::cout << "=== Test 3: NumsModK round-trip ===" << std::endl;

    // NUM_BITS=32 => virtualMaxLevel=6, CFLOBDDMaxLevel vars
    unsigned int numVars = 1u << CFLOBDDMaxLevel;

    unsigned int testModuli[] = { 3, 5, 7, 11, 59 };

    for (unsigned int k : testModuli) {
        // Build CFLOBDD NumsModK (position A = first half of variables)
        CFLOBDD orig = NumsModK(k, CFL_OBDD::A);

        // Convert CFLOBDD -> ADD -> CFLOBDD
        Cudd mgr(0, 0);
        for (unsigned int i = 0; i < numVars; i++) mgr.addVar(i);

        auto start = high_resolution_clock::now();
        ADD add = CFLOBDD_to_ADD<int>(mgr, orig);
        auto mid = high_resolution_clock::now();
        CFLOBDD rt = ADD_to_CFLOBDD<int>(mgr, add);
        auto end = high_resolution_clock::now();

        auto cfToAdd = duration_cast<milliseconds>(mid - start).count();
        auto addToCf = duration_cast<milliseconds>(end - mid).count();

        // Structural equality check on CFLOBDDs
        if (orig == rt) {
            std::cout << "  3. NumsModK(k=" << k << "): PASS (structural)"
                      << "  [CF->ADD:" << cfToAdd << "ms, ADD->CF:" << addToCf << "ms]" << std::endl;
        } else {
            std::cout << "  3. NumsModK(k=" << k << "): structural FAIL, checking exhaustively..." << std::endl;
            if (!ExhaustiveCompare(add, orig, numVars, "NumsModK")) return false;
            std::cout << "  3. NumsModK(k=" << k << "): PASS (semantic)" << std::endl;
        }
    }

    return true;
}

// =========================================================================
// Test 4: MultModK round-trip (CRT primitive)
//
// Build MultModK as a CFLOBDD, convert to ADD and back, check structural
// equality on the CFLOBDD side.
// =========================================================================
static bool testMultModK() {
    std::cout << "=== Test 4: MultModK round-trip ===" << std::endl;

    unsigned int numVars = 1u << CFLOBDDMaxLevel;

    // Small moduli — MultModK is more expensive than NumsModK
    unsigned int testModuli[] = { 3, 5, 7 };

    for (unsigned int k : testModuli) {
        // Build CFLOBDD MultModK
        auto start = high_resolution_clock::now();
        CFLOBDD orig = MultModK(k);
        auto mid1 = high_resolution_clock::now();

        // Convert CFLOBDD -> ADD -> CFLOBDD
        Cudd mgr(0, 0);
        for (unsigned int i = 0; i < numVars; i++) mgr.addVar(i);

        ADD add = CFLOBDD_to_ADD<int>(mgr, orig);
        auto mid2 = high_resolution_clock::now();
        CFLOBDD rt = ADD_to_CFLOBDD<int>(mgr, add);
        auto end = high_resolution_clock::now();

        auto buildTime = duration_cast<milliseconds>(mid1 - start).count();
        auto cfToAdd = duration_cast<milliseconds>(mid2 - mid1).count();
        auto addToCf = duration_cast<milliseconds>(end - mid2).count();

        // Structural equality check on CFLOBDDs
        if (orig == rt) {
            std::cout << "  4. MultModK(k=" << k << "): PASS (structural)"
                      << "  [build:" << buildTime << "ms, CF->ADD:" << cfToAdd
                      << "ms, ADD->CF:" << addToCf << "ms]" << std::endl;
        } else {
            std::cout << "  4. MultModK(k=" << k << "): structural FAIL" << std::endl;
            return false;
        }
    }

    return true;
}
#endif // USE_CRT_TESTS

// =========================================================================
// Test 5: Round-trip with larger variable counts (structural only)
// =========================================================================
static bool testScaling() {
    std::cout << "=== Test 5: Scaling (structural round-trip) ===" << std::endl;

    // Test round-trip ADD -> CFLOBDD -> ADD at increasing sizes
    unsigned int sizes[] = { 8, 16, 32, 64 };

    for (unsigned int n : sizes) {
        Cudd mgr(0, 0);
        for (unsigned int i = 0; i < n; i++) mgr.addVar(i);

        // Build ADD: weighted sum (each var weighted by its index mod 13)
        ADD f = mgr.constant(0.0);
        for (unsigned int i = 0; i < n; i++) {
            f = f + mgr.addVar(i) * mgr.constant((double)(i % 13));
        }

        auto start = high_resolution_clock::now();
        CFLOBDD cf = ADD_to_CFLOBDD<int>(mgr, f);
        auto mid = high_resolution_clock::now();
        ADD rt = CFLOBDD_to_ADD<int>(mgr, cf);
        auto end = high_resolution_clock::now();

        auto addToCf = duration_cast<milliseconds>(mid - start).count();
        auto cfToAdd = duration_cast<milliseconds>(end - mid).count();

        if (f == rt) {
            std::cout << "  5. " << n << " vars: PASS (structural)"
                      << "  [ADD->CF:" << addToCf << "ms, CF->ADD:" << cfToAdd << "ms]" << std::endl;
        } else {
            std::cout << "  5. " << n << " vars: structural FAIL" << std::endl;
            return false;
        }
    }

    return true;
}

// =========================================================================
// main
// =========================================================================
int main() {
    InitCFLOBDD();

    bool allPassed = true;

    allPassed &= testRoundTrip_ADD();
    allPassed &= testRoundTrip_CFLOBDD();
#ifdef USE_CRT_TESTS
    allPassed &= testNumsModK();
    allPassed &= testMultModK();
#endif
    allPassed &= testScaling();

    if (allPassed) {
        std::cout << "\n*** ALL TESTS PASSED ***" << std::endl;
        return 0;
    } else {
        std::cout << "\n*** SOME TESTS FAILED ***" << std::endl;
        return 1;
    }
}
