/**
  @file addMultiplication.cc

  @brief Implementation of ADD-based multiplication modulo k.

  @details This file implements the functions declared in addMultiplication.hh
  for building ADDs that represent multiplication modulo k.

  The implementation uses CUDD's addResidue function to build NumsModK, and
  Cudd_addApply with a custom operator for the multiplication mod k combination.

  This is the ADD equivalent of the CFLOBDD implementation in multiplication_crt.cpp,
  which uses ApplyAndReduce with MultiplyModKFunc.

*/

#include "addMultiplication.hh"
#include "../cudd-3.0.0/cudd/cuddInt.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <chrono>
#include <string>
#include <iomanip>
#include <cctype>

using namespace std::chrono;


// =============================================================================
// Global variable for the apply operator
// =============================================================================

// Verbose flag: controls printing of intermediate structure sizes and
// cache diagnostics.  Set via --verbose / -v on the command line.
static bool verbose = false;

// Global variable to store the current modulus for the apply operator.
// This is necessary because CUDD's apply interface doesn't allow passing
// additional parameters to the operator function.
// This mirrors the CFLOBDD approach in multiplication_crt.cpp which uses
// "static unsigned int currentModulus".
static int currentModulusK = 0;

// Track the last modulus used - only flush cache when modulus changes.
// Flushing on every operation causes exponential blowup because Cudd_addApply
// loses its memoization and recomputes subproblems repeatedly.
static int lastCachedModulusK = -1;

/**
  @brief Sets the current modulus and flushes cache only if modulus changed.

  @details The CUDD cache is keyed on (operator, node1, node2), but our apply
  operators depend on the global currentModulusK. When the modulus changes,
  cached results become invalid. However, flushing on EVERY operation causes
  exponential blowup because Cudd_addApply relies on caching for efficiency.

  This function only flushes when the modulus actually changes.
*/
static void setModulusK(DdManager* ddMgr, int k)
{
    currentModulusK = k;
    if (k != lastCachedModulusK) {
        cuddCacheFlush(ddMgr);
        lastCachedModulusK = k;
    }
}


// =============================================================================
// Custom apply operator for multiplication mod k
// =============================================================================

/**
  @brief Apply operator for multiplication modulo k.

  @details This is used with Cudd_addApply to combine two ADDs element-wise
  with the operation (f * g) mod k, where k is stored in currentModulusK.

  This is analogous to MultiplyModKFunc in the CFLOBDD implementation.

  @param dd       CUDD manager
  @param f        Pointer to first operand node
  @param g        Pointer to second operand node

  @return Result node if both operands are constants; NULL otherwise.

*/
static DdNode* addTimesModK(DdManager* dd, DdNode** f, DdNode** g)
{
    DdNode* F = *f;
    DdNode* G = *g;

    // Short-circuit: 0 * x = 0
    if (F == DD_ZERO(dd)) return F;
    if (G == DD_ZERO(dd)) return G;

    // Short-circuit: 1 * x = x (only if 1 is a valid value, i.e., k > 1)
    if (cuddIsConstant(F) && cuddV(F) == 1.0) return G;
    if (cuddIsConstant(G) && cuddV(G) == 1.0) return F;

    // Terminal case: both operands are constants
    if (cuddIsConstant(F) && cuddIsConstant(G)) {
        CUDD_VALUE_TYPE valueF = cuddV(F);
        CUDD_VALUE_TYPE valueG = cuddV(G);
        // Compute (f * g) mod k
        int intF = (int)valueF;
        int intG = (int)valueG;
        int result = (intF * intG) % currentModulusK;
        return cuddUniqueConst(dd, (CUDD_VALUE_TYPE)result);
    }

    // For commutativity: ensure consistent ordering for cache efficiency
    if (F > G) {
        *f = G;
        *g = F;
    }

    return NULL;  // Not a terminal case; CUDD will recurse
}


/**
  @brief Apply operator for addition modulo k.
*/
static DdNode* addPlusModK(DdManager* dd, DdNode** f, DdNode** g)
{
    DdNode* F = *f;
    DdNode* G = *g;

    // Short-circuit: 0 + x = x (values are already in [0, k-1])
    if (F == DD_ZERO(dd)) return G;
    if (G == DD_ZERO(dd)) return F;

    // Terminal case: both operands are constants
    if (cuddIsConstant(F) && cuddIsConstant(G)) {
        CUDD_VALUE_TYPE valueF = cuddV(F);
        CUDD_VALUE_TYPE valueG = cuddV(G);
        int intF = (int)valueF;
        int intG = (int)valueG;
        int result = (intF + intG) % currentModulusK;
        return cuddUniqueConst(dd, (CUDD_VALUE_TYPE)result);
    }

    // For commutativity: ensure consistent ordering for cache efficiency
    if (F > G) {
        *f = G;
        *g = F;
    }

    return NULL;
}


/**
  @brief Apply operator for subtraction modulo k.
*/
static DdNode* addMinusModK(DdManager* dd, DdNode** f, DdNode** g)
{
    DdNode* F = *f;
    DdNode* G = *g;

    // Short-circuit: x - 0 = x
    if (G == DD_ZERO(dd)) return F;

    // Short-circuit: x - x = 0
    if (F == G) return DD_ZERO(dd);

    // Terminal case: both operands are constants
    if (cuddIsConstant(F) && cuddIsConstant(G)) {
        CUDD_VALUE_TYPE valueF = cuddV(F);
        CUDD_VALUE_TYPE valueG = cuddV(G);
        int intF = (int)valueF;
        int intG = (int)valueG;
        int result = (intF + currentModulusK - intG) % currentModulusK;
        return cuddUniqueConst(dd, (CUDD_VALUE_TYPE)result);
    }

    return NULL;
}


// =============================================================================
// Main functions
// =============================================================================

ADD ADD_NumsModK(Cudd& mgr, int n, int k, int top)
{
    assert(k >= 2);
    assert(n >= 1);

    // Use CUDD's built-in addResidue function.
    // CUDD_RESIDUE_MSB means MSB is at the top of the ADD (high-order to low-order).
    return mgr.addResidue(n, k, CUDD_RESIDUE_MSB, top);
}


ADD ADD_MultModK(Cudd& mgr, int n, int k)
{
    assert(k >= 2);
    assert(n >= 1);

    // Build NumsModK for x (variables 0..n-1)
    ADD numsModK_x = ADD_NumsModK(mgr, n, k, 0);

    // Build NumsModK for y (variables n..2n-1)
    ADD numsModK_y = ADD_NumsModK(mgr, n, k, n);

    // Set the global modulus for the apply operator (flushes cache only if k changed)
    DdManager* ddMgr = mgr.getManager();
    setModulusK(ddMgr, k);

    // Combine using multiplication mod k via Cudd_addApply
    // This is analogous to CFLOBDD's:
    //   return CFLOBDD(ApplyAndReduce<int>(cA.root, cB.root, MultiplyModKFunc));
    DdNode* result = Cudd_addApply(ddMgr, addTimesModK,
                                   numsModK_x.getNode(),
                                   numsModK_y.getNode());

    if (result == NULL) {
        std::cerr << "Error: Cudd_addApply returned NULL" << std::endl;
        return mgr.addZero();
    }

    return ADD(mgr, result);
}


CUDD_VALUE_TYPE ADD_Evaluate(Cudd& mgr, const ADD& f, int n,
                              INPUT_TYPE x, INPUT_TYPE y)
{
    (void)mgr;  // Unused, but kept for API consistency

    // Create an input array for CUDD's Eval function
    // The array has one entry per variable: 0 or 1
    int numVars = 2 * n;
    int* inputs = new int[numVars];

    // Fill in the x bits (variables 0..n-1)
    // High-order bit at lower index (MSB at variable 0)
    for (int i = 0; i < n; i++) {
        int bitPos = n - 1 - i;  // High-order to low-order
        inputs[i] = static_cast<int>((x >> bitPos) & 1);
    }

    // Fill in the y bits (variables n..2n-1)
    for (int i = 0; i < n; i++) {
        int bitPos = n - 1 - i;  // High-order to low-order
        inputs[n + i] = static_cast<int>((y >> bitPos) & 1);
    }

    // Evaluate the ADD - this returns the terminal node for this path
    ADD result = f.Eval(inputs);
    delete[] inputs;

    // Extract the constant value from the result
    // The result of Eval on a valid path should be a constant
    DdNode* node = result.getNode();
    return cuddV(Cudd_Regular(node));
}


bool verifyMultModK(Cudd& mgr, int n, int k)
{
    std::cout << "Building ADD_MultModK(n=" << n << ", k=" << k << ")..." << std::endl;

    high_resolution_clock::time_point start = high_resolution_clock::now();
    ADD multModK = ADD_MultModK(mgr, n, k);
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> buildTime = duration_cast<duration<double>>(end - start);

    std::cout << "Build time: " << buildTime.count() << " seconds" << std::endl;
    std::cout << "Node count: " << multModK.nodeCount() << std::endl;
    std::cout << "Leaf count: " << multModK.CountLeaves() << std::endl;

    // Test a subset of inputs
    INPUT_TYPE maxVal = INPUT_TYPE(1) << n;
    unsigned long long testCount = 0;
    unsigned long long errorCount = 0;

    // For small n, test all inputs; for large n, test a sample
    INPUT_TYPE step = (maxVal <= 256) ? INPUT_TYPE(1) : maxVal / 64;

    std::cout << "Verifying results..." << std::endl;

    start = high_resolution_clock::now();

    for (INPUT_TYPE x = 0; x < maxVal; x += step) {
        for (INPUT_TYPE y = 0; y < maxVal; y += step) {
            CUDD_VALUE_TYPE addResult = ADD_Evaluate(mgr, multModK, n, x, y);
            unsigned long long expected = static_cast<unsigned long long>(
                (OUTPUT_TYPE(x) * OUTPUT_TYPE(y)) % k);

            testCount++;

            if ((unsigned long long)addResult != expected) {
                errorCount++;
                if (errorCount <= 10) {
                    std::cout << "ERROR: x=" << x << ", y=" << y
                              << ", ADD result=" << addResult
                              << ", expected=" << expected << std::endl;
                }
            }
        }
    }

    end = high_resolution_clock::now();
    duration<double> verifyTime = duration_cast<duration<double>>(end - start);

    std::cout << "Verify time: " << verifyTime.count() << " seconds" << std::endl;
    std::cout << "Tested " << testCount << " input pairs" << std::endl;

    if (errorCount == 0) {
        std::cout << "All tests PASSED!" << std::endl;
        return true;
    } else {
        std::cout << errorCount << " tests FAILED!" << std::endl;
        return false;
    }
}


void printADDStats(const ADD& f, const char* name)
{
    std::cout << name << ": "
              << "nodes=" << f.nodeCount()
              << ", leaves=" << f.CountLeaves()
              << std::endl;
}


// =============================================================================
// Modular arithmetic functions
// =============================================================================

ADD ADD_AddModK(Cudd& mgr, const ADD& a, const ADD& b, int k)
{
    DdManager* ddMgr = mgr.getManager();
    setModulusK(ddMgr, k);

    DdNode* result = Cudd_addApply(ddMgr, addPlusModK, a.getNode(), b.getNode());
    if (result == NULL) {
        std::cerr << "Error: ADD_AddModK failed" << std::endl;
        return mgr.addZero();
    }

    return ADD(mgr, result);
}


ADD ADD_SubtractModK(Cudd& mgr, const ADD& a, const ADD& b, int k)
{
    DdManager* ddMgr = mgr.getManager();
    setModulusK(ddMgr, k);

    DdNode* result = Cudd_addApply(ddMgr, addMinusModK, a.getNode(), b.getNode());
    if (result == NULL) {
        std::cerr << "Error: ADD_SubtractModK failed" << std::endl;
        return mgr.addZero();
    }

    return ADD(mgr, result);
}


ADD ADD_MultiplyModK(Cudd& mgr, const ADD& a, const ADD& b, int k)
{
    DdManager* ddMgr = mgr.getManager();
    setModulusK(ddMgr, k);

    DdNode* result = Cudd_addApply(ddMgr, addTimesModK, a.getNode(), b.getNode());
    if (result == NULL) {
        std::cerr << "Error: ADD_MultiplyModK failed" << std::endl;
        return mgr.addZero();
    }

    return ADD(mgr, result);
}


// =============================================================================
// Positioned NumsModK for Karatsuba
// =============================================================================

ADD ADD_NumsModKPositioned(Cudd& mgr, int n, int k, Position p)
{
    assert(k >= 2);
    assert(p != TopLevel);

    int halfN = n / 2;

    switch (p) {
        case A:
            // First input: variables 0..n-1
            return ADD_NumsModK(mgr, n, k, 0);

        case B:
            // Second input: variables n..2n-1
            return ADD_NumsModK(mgr, n, k, n);

        case AA:
            // High half of first input: variables 0..halfN-1
            return ADD_NumsModK(mgr, halfN, k, 0);

        case AB:
            // Low half of first input: variables halfN..n-1
            return ADD_NumsModK(mgr, halfN, k, halfN);

        case BA:
            // High half of second input: variables n..n+halfN-1
            return ADD_NumsModK(mgr, halfN, k, n);

        case BB:
            // Low half of second input: variables n+halfN..2n-1
            return ADD_NumsModK(mgr, halfN, k, n + halfN);

        default:
            assert(false);
            return mgr.addZero();
    }
}


// =============================================================================
// Shift-and-Add Multiplication
// =============================================================================

ADD ADD_ShiftAndAddMultiplicationModK(Cudd& mgr, int n, int k)
{
    // Note: currentModulusK is set by the first call to ADD_MultiplyModK below

    // Get the second input (b): val(y) mod k for variables n..2n-1
    ADD b = ADD_NumsModKPositioned(mgr, n, k, B);

    // Initialize result to 0
    ADD result = mgr.constant(0);
    ADD two = mgr.constant(2 % k);

    // Process each bit of first input from high-order to low-order
    // i = 0 is MSB, i = n-1 is LSB
    for (int i = 0; i < n; i++) {
        // Get projection of bit i: returns 0 or 1
        ADD p = mgr.addVar(i);

        // addend = p * b mod k (either 0 or b, depending on bit value)
        ADD addend = ADD_MultiplyModK(mgr, p, b, k);

        // result = result * 2 mod k (shift left)
        ADD temp = ADD_MultiplyModK(mgr, two, result, k);

        // result = temp + addend mod k
        result = ADD_AddModK(mgr, temp, addend, k);
    }

    return result;
}


// =============================================================================
// Subtractive Karatsuba One Level
// =============================================================================

ADD ADD_SubtractiveKaratsubaOneLevel(Cudd& mgr, int n, int k)
{
    assert(n % 2 == 0);  // n must be even for Karatsuba

    // m = half of the number of bits in an input
    int m = n / 2;

    // Compute shift constants: 2^m mod k and 2^(2m) mod k
    unsigned int shift_m = 1;
    for (int i = 0; i < m; i++) {
        shift_m = (shift_m * 2) % k;
    }
    ADD const_m = mgr.constant(shift_m);

    unsigned int shift_2m = (shift_m * shift_m) % k;
    ADD const_2m = mgr.constant(shift_2m);

    // Get half-size number representations
    // X = X1*2^m + X0 (first input split into high/low halves)
    // Y = Y1*2^m + Y0 (second input split into high/low halves)
    if (verbose) std::cerr << "  [Karatsuba] building X1..." << std::flush;
    ADD X1 = ADD_NumsModKPositioned(mgr, n, k, AA);  // High bits of X
    if (verbose) std::cerr << " X0..." << std::flush;
    ADD X0 = ADD_NumsModKPositioned(mgr, n, k, AB);  // Low bits of X
    if (verbose) std::cerr << " Y1..." << std::flush;
    ADD Y1 = ADD_NumsModKPositioned(mgr, n, k, BA);  // High bits of Y
    if (verbose) std::cerr << " Y0..." << std::flush;
    ADD Y0 = ADD_NumsModKPositioned(mgr, n, k, BB);  // Low bits of Y
    if (verbose) std::cerr << " done" << std::endl;

    // Compute three products using subtractive Karatsuba formulation
    if (verbose) std::cerr << "  [Karatsuba] Z0 = X0*Y0..." << std::flush;
    ADD Z0 = ADD_MultiplyModK(mgr, X0, Y0, k);
    if (verbose) std::cerr << " done (" << Z0.nodeCount() << " nodes)" << std::endl;

    if (verbose) std::cerr << "  [Karatsuba] Z2 = X1*Y1..." << std::flush;
    ADD Z2 = ADD_MultiplyModK(mgr, X1, Y1, k);
    if (verbose) std::cerr << " done (" << Z2.nodeCount() << " nodes)" << std::endl;

    if (verbose) std::cerr << "  [Karatsuba] X_diff = X1-X0..." << std::flush;
    ADD X_diff = ADD_SubtractModK(mgr, X1, X0, k);
    if (verbose) std::cerr << " done (" << X_diff.nodeCount() << " nodes)" << std::endl;

    if (verbose) std::cerr << "  [Karatsuba] Y_diff = Y0-Y1..." << std::flush;
    ADD Y_diff = ADD_SubtractModK(mgr, Y0, Y1, k);
    if (verbose) std::cerr << " done (" << Y_diff.nodeCount() << " nodes)" << std::endl;

    if (verbose) std::cerr << "  [Karatsuba] Z1 = X_diff*Y_diff..." << std::flush;
    ADD Z1 = ADD_MultiplyModK(mgr, X_diff, Y_diff, k);
    if (verbose) std::cerr << " done (" << Z1.nodeCount() << " nodes)" << std::endl;

    // Combine: result = Z2*2^(2m) + (Z1 + Z2 + Z0)*2^m + Z0
    DdManager* ddMgr = mgr.getManager();
    if (verbose) std::cerr << "  [Karatsuba] middle = Z1+Z2..." << std::flush;
    ADD middle = ADD_AddModK(mgr, Z1, Z2, k);
    if (verbose) std::cerr << " done (" << middle.nodeCount() << " nodes)"
                           << " [cache: " << Cudd_ReadCacheSlots(ddMgr) << " slots]"
                           << std::endl;
    if (verbose) std::cerr << "  [Karatsuba] middle += Z0 (" << Z0.nodeCount()
                           << " nodes)..." << std::flush;
    middle = ADD_AddModK(mgr, middle, Z0, k);
    if (verbose) std::cerr << " done (" << middle.nodeCount() << " nodes)" << std::endl;

    if (verbose) std::cerr << "  [Karatsuba] term_Z2 = Z2*2^(2m)..." << std::flush;
    ADD term_Z2 = ADD_MultiplyModK(mgr, Z2, const_2m, k);
    if (verbose) std::cerr << " done" << std::endl;

    if (verbose) std::cerr << "  [Karatsuba] term_middle = middle*2^m..." << std::flush;
    ADD term_middle = ADD_MultiplyModK(mgr, middle, const_m, k);
    if (verbose) std::cerr << " done" << std::endl;

    // Final result: term_Z2 + term_middle + Z0
    if (verbose) std::cerr << "  [Karatsuba] final = term_Z2+term_middle..." << std::flush;
    ADD karatsuba = ADD_AddModK(mgr, term_Z2, term_middle, k);
    if (verbose) std::cerr << " done (" << karatsuba.nodeCount() << " nodes)" << std::endl;
    if (verbose) std::cerr << "  [Karatsuba] final += Z0..." << std::flush;
    karatsuba = ADD_AddModK(mgr, karatsuba, Z0, k);
    if (verbose) std::cerr << " done (" << karatsuba.nodeCount() << " nodes)" << std::endl;

    return karatsuba;
}


// =============================================================================
// Verification functions
// =============================================================================

bool ADD_VerifyShiftAndAddMultiplicationModK(Cudd& mgr, int n, int k)
{
    std::cout << "Verifying shift-and-add multiplication mod " << k << "..." << std::endl;

    auto start = high_resolution_clock::now();

    // Build specification (direct MultModK)
    ADD specification = ADD_MultModK(mgr, n, k);

    // Build shift-and-add result
    ADD result = ADD_ShiftAndAddMultiplicationModK(mgr, n, k);

    // Compare
    bool equal = (specification == result);

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);

    std::cout << "  Result: " << (equal ? "PASSED" : "FAILED")
              << " (" << duration.count() << " ms)" << std::endl;

    return equal;
}


void ADD_BuildMultiplicationSpecsModuliwise(Cudd& mgr, int n)
{
    std::cout << "Building specification ADDs for all " << numberOfMultRelations
              << " moduli for " << n << "-bit operands..." << std::endl;

    auto totalStart = high_resolution_clock::now();

    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        std::cout << "  Modulus " << Moduli[i] << " (" << i+1 << "/"
                  << numberOfMultRelations << ")..." << std::flush;
        if (i == numberOfMultRelations - 1) {
            auto start = high_resolution_clock::now();
            ADD spec = ADD_MultModK(mgr, n, Moduli[i]);
            auto end = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end - start);
            std::cout << " nodes=" << spec.nodeCount()
                      << ", leaves=" << spec.CountLeaves()
                      << " (" << duration.count() << " ms)" << std::endl;
        } else {
            ADD spec = ADD_MultModK(mgr, n, Moduli[i]);
            std::cout << " nodes=" << spec.nodeCount()
                      << ", leaves=" << spec.CountLeaves() << std::endl;
        }
    }

    auto totalEnd = high_resolution_clock::now();
    auto totalDuration = duration_cast<milliseconds>(totalEnd - totalStart);
    std::cout << "ADD_BuildMultiplicationSpecsModuliwise took " << totalDuration.count() << " ms" << std::endl;
}


bool ADD_VerifyShiftAndAddMultiplicationModuliwise(Cudd& mgr, int n)
{
    std::cout << "Verifying shift-and-add across all " << numberOfMultRelations
              << " moduli for " << n << "-bit operands..." << std::endl;

    auto start = high_resolution_clock::now();

    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        std::cout << "  Modulus " << Moduli[i] << " (" << i+1 << "/"
                  << numberOfMultRelations << ")..." << std::flush;

        if (i == numberOfMultRelations - 1) {
            auto lastStart = high_resolution_clock::now();
            ADD spec = ADD_MultModK(mgr, n, Moduli[i]);
            ADD result = ADD_ShiftAndAddMultiplicationModK(mgr, n, Moduli[i]);
            auto lastEnd = high_resolution_clock::now();
            if (!(spec == result)) {
                std::cout << " FAILED" << std::endl;
                return false;
            }
            auto lastDuration = duration_cast<milliseconds>(lastEnd - lastStart);
            std::cout << " OK (" << lastDuration.count() << " ms)" << std::endl;
        } else {
            ADD spec = ADD_MultModK(mgr, n, Moduli[i]);
            ADD result = ADD_ShiftAndAddMultiplicationModK(mgr, n, Moduli[i]);
            if (!(spec == result)) {
                std::cout << " FAILED" << std::endl;
                return false;
            }
            std::cout << " OK" << std::endl;
        }
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);

    std::cout << "SUCCESS! All moduli verified in " << duration.count() << " ms" << std::endl;
    return true;
}


bool ADD_VerifySubtractiveKaratsubaOneLevel(Cudd& mgr, int n, int k)
{
    std::cout << "Verifying Karatsuba one-level mod " << k << "..." << std::endl;

    DdManager* ddMgr = mgr.getManager();

    double lookupsBefore = 0, hitsBefore = 0;
    int gcBefore = 0;
    if (verbose) {
        std::cout << "  Cache slots: " << Cudd_ReadCacheSlots(ddMgr)
                  << ", max cache hard: " << Cudd_ReadMaxCacheHard(ddMgr)
                  << std::endl;
        lookupsBefore = Cudd_ReadCacheLookUps(ddMgr);
        hitsBefore = Cudd_ReadCacheHits(ddMgr);
        gcBefore = Cudd_ReadGarbageCollections(ddMgr);
    }

    auto start = high_resolution_clock::now();

    // Build specification (direct MultModK)
    ADD specification = ADD_MultModK(mgr, n, k);

    auto midpoint = high_resolution_clock::now();

    // Build Karatsuba result
    ADD karatsuba = ADD_SubtractiveKaratsubaOneLevel(mgr, n, k);

    auto end = high_resolution_clock::now();

    // Compare
    bool equal = (specification == karatsuba);

    auto specTime = duration_cast<milliseconds>(midpoint - start);
    auto karatsubaTime = duration_cast<milliseconds>(end - midpoint);
    auto totalTime = duration_cast<milliseconds>(end - start);

    if (verbose) {
        double lookupsAfter = Cudd_ReadCacheLookUps(ddMgr);
        double hitsAfter = Cudd_ReadCacheHits(ddMgr);
        int gcAfter = Cudd_ReadGarbageCollections(ddMgr);

        double lookups = lookupsAfter - lookupsBefore;
        double hits = hitsAfter - hitsBefore;
        double hitRate = (lookups > 0) ? (hits / lookups * 100.0) : 0.0;

        std::cout << "  MultModK time: " << specTime.count() << " ms"
                  << ", Karatsuba time: " << karatsubaTime.count() << " ms" << std::endl;
        std::cout << "  Cache lookups: " << (long long)lookups
                  << ", hits: " << (long long)hits
                  << ", hit rate: " << std::fixed << std::setprecision(1) << hitRate << "%"
                  << std::defaultfloat << std::endl;
        std::cout << "  GC collections: " << (gcAfter - gcBefore)
                  << ", dead nodes: " << Cudd_ReadDead(ddMgr) << std::endl;
    }

    std::cout << "  Result: " << (equal ? "PASSED" : "FAILED")
              << " (" << totalTime.count() << " ms)" << std::endl;

    return equal;
}


bool ADD_VerifySubtractiveKaratsubaOneLevelModuliwise(Cudd& mgr, int n)
{
    std::cout << "Verifying Karatsuba across all " << numberOfMultRelations
              << " moduli for " << n << "-bit operands..." << std::endl;

    auto start = high_resolution_clock::now();

    // Process in forward order (smallest primes first).
    // ADDs blow up for larger k (intermediate results grow as O(k^4*n)),
    // so do small moduli first to get results before hitting the limit.
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        std::cout << "  Modulus " << Moduli[i] << " (" << (i + 1) << "/"
                  << numberOfMultRelations << ")..." << std::flush;

        if (i == numberOfMultRelations - 1) {
            auto lastStart = high_resolution_clock::now();
            ADD spec = ADD_MultModK(mgr, n, Moduli[i]);
            ADD karatsuba = ADD_SubtractiveKaratsubaOneLevel(mgr, n, Moduli[i]);
            auto lastEnd = high_resolution_clock::now();
            if (!(spec == karatsuba)) {
                std::cout << " FAILED" << std::endl;
                return false;
            }
            auto lastDuration = duration_cast<milliseconds>(lastEnd - lastStart);
            std::cout << " OK (" << lastDuration.count() << " ms)" << std::endl;
        } else {
            ADD spec = ADD_MultModK(mgr, n, Moduli[i]);
            ADD karatsuba = ADD_SubtractiveKaratsubaOneLevel(mgr, n, Moduli[i]);
            if (!(spec == karatsuba)) {
                std::cout << " FAILED" << std::endl;
                return false;
            }
            std::cout << " OK" << std::endl;
        }
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);

    std::cout << "SUCCESS! All moduli verified in " << duration.count() << " ms" << std::endl;
    return true;
}


// =============================================================================
// First 999 odd primes (from multiplication_crt.cpp)
// =============================================================================

const unsigned int AllModuli[numberOfOddPrimes] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
    359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
    479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661,
    673, 677, 683, 691, 701, 709, 719, 727, 733, 739,
    743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
    821, 823, 827, 829, 839, 853, 857, 859, 863, 877,
    881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
    953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019,
    1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087,
    1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153,
    1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229,
    1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297,
    1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381,
    1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453,
    1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523,
    1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597,
    1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663,
    1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741,
    1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823,
    1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901,
    1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993,
    1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063,
    2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131,
    2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221,
    2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293,
    2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371,
    2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437,
    2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539,
    2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621,
    2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689,
    2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749,
    2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833,
    2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909,
    2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001,
    3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083,
    3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187,
    3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259,
    3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343,
    3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433,
    3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517,
    3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581,
    3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659,
    3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733,
    3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823,
    3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911,
    3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001,
    4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073,
    4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153,
    4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241,
    4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327,
    4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421,
    4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507,
    4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591,
    4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663,
    4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759,
    4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861,
    4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943,
    4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009,
    5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099,
    5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189,
    5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281,
    5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393,
    5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449,
    5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527,
    5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, 5641,
    5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701,
    5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801,
    5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861,
    5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953,
    5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067,
    6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143,
    6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229,
    6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311,
    6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373,
    6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481,
    6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577,
    6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679,
    6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763,
    6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841,
    6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947,
    6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001,
    7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109,
    7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211,
    7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307,
    7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417,
    7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507,
    7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573,
    7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649,
    7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727,
    7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841,
    7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919
};

/// Pointer to moduli array (for backward compatibility)
const unsigned int* Moduli = AllModuli;


// =============================================================================
// MultRelation implementation
// =============================================================================

MultRelation::MultRelation(Cudd& manager)
    : mgr(&manager)
{
    std::cout << "Building MultRelation with " << numberOfMultRelations
              << " primes for " << numBits << "-bit operands..." << std::endl;

    auto start = high_resolution_clock::now();

    // Initialize ModuliArray with first 26 odd primes
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        ModuliArray[i] = Moduli[i];
    }

    // Build ADDs for each prime
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        ModularMultRelations[i] = ADD_MultModK(*mgr, numBits, ModuliArray[i]);
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);

    std::cout << "MultRelation built in " << duration.count() << " ms" << std::endl;
}


MultRelation::MultRelation()
    : mgr(new Cudd(0, 0))
{
    std::cout << "Building MultRelation with " << numberOfMultRelations
              << " primes for " << numBits << "-bit operands..." << std::endl;

    auto start = high_resolution_clock::now();

    // Initialize ModuliArray with first 26 odd primes
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        ModuliArray[i] = Moduli[i];
    }

    // Build ADDs for each prime
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        ModularMultRelations[i] = ADD_MultModK(*mgr, numBits, ModuliArray[i]);
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);

    std::cout << "MultRelation built in " << duration.count() << " ms" << std::endl;
}


MultRelation::~MultRelation()
{
    // Note: We don't delete mgr here because it might be shared
    // In a production implementation, we'd use shared_ptr
}


bool MultRelation::operator==(const MultRelation& other) const
{
    // Compare each entry in ModuliArray
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        if (ModuliArray[i] != other.ModuliArray[i]) {
            return false;
        }
    }

    // Compare each ADD
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        if (ModularMultRelations[i] != other.ModularMultRelations[i]) {
            return false;
        }
    }

    return true;
}


// Helper: Extended Euclidean Algorithm for computing modular inverse
static unsigned int ModularInverse(unsigned int a, unsigned int m)
{
    int m0 = m;
    int y = 0, x = 1;

    if (m == 1) return 0;

    while (a > 1) {
        int q = a / m;
        int t = m;

        m = a % m;
        a = t;
        t = y;

        y = x - q * y;
        x = t;
    }

    if (x < 0) x += m0;

    return x;
}


OUTPUT_TYPE MultRelation::Multiply(INPUT_TYPE a, INPUT_TYPE b) const
{
    // Evaluate each ADD to get remainders modulo each prime
    unsigned int remainders[numberOfMultRelations];
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        remainders[i] = (unsigned int)ADD_Evaluate(*mgr, ModularMultRelations[i],
                                                    numBits, a, b);
    }

    // Apply Garner's algorithm to reconstruct the result
    // Garner's algorithm computes the mixed radix representation:
    // x = u[0] + u[1]*m[0] + u[2]*m[0]*m[1] + ... + u[n-1]*m[0]*...*m[n-2]

    unsigned int u[numberOfMultRelations];

    // Compute u values
    u[0] = remainders[0];
    for (unsigned int i = 1; i < numberOfMultRelations; i++) {
        u[i] = remainders[i];
        for (unsigned int j = 0; j < i; j++) {
            // u[i] = (u[i] - u[j]) * inverse(m[j], m[i]) mod m[i]
            int diff = u[i] - u[j];
            if (diff < 0) diff += ModuliArray[i];
            unsigned int inv = ModularInverse(ModuliArray[j], ModuliArray[i]);
            u[i] = (diff * inv) % ModuliArray[i];
        }
    }

    // Reconstruct the result using mixed radix representation
    OUTPUT_TYPE result = u[0];
    OUTPUT_TYPE product = ModuliArray[0];

    for (unsigned int i = 1; i < numberOfMultRelations; i++) {
        result += OUTPUT_TYPE(u[i]) * product;
        product *= ModuliArray[i];
    }

    return result;
}


void MultRelation::PrintStats() const
{
    std::cout << "MultRelation Statistics:" << std::endl;
    std::cout << "  Number of relations: " << numberOfMultRelations << std::endl;
    std::cout << "  Bits per operand: " << numBits << std::endl;
    std::cout << std::endl;

    unsigned long totalNodes = 0;
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        int nodes = ModularMultRelations[i].nodeCount();
        int leaves = ModularMultRelations[i].CountLeaves();
        totalNodes += nodes;
        std::cout << "  Modulus " << ModuliArray[i] << ": "
                  << nodes << " nodes, " << leaves << " leaves" << std::endl;
    }
    std::cout << std::endl;
    std::cout << "  Total nodes: " << totalNodes << std::endl;
}


// =============================================================================
// Main function for testing
// =============================================================================

void testMultRelation()
{
    std::cout << "=== Testing MultRelation (CRT-based " << NUM_BITS << "-bit multiplication) ===" << std::endl;
    std::cout << std::endl;

    // Build the MultRelation with all primes
    MultRelation R;

    std::cout << std::endl;
    R.PrintStats();
    std::cout << std::endl;

    // Test multiplication with some sample values
    std::cout << "--- Verifying multiplication ---" << std::endl;

    struct TestCase {
        INPUT_TYPE a;
        INPUT_TYPE b;
    };

    TestCase tests[] = {
        {INPUT_TYPE(0), INPUT_TYPE(0)},
        {INPUT_TYPE(1), INPUT_TYPE(1)},
        {INPUT_TYPE(2), INPUT_TYPE(3)},
        {INPUT_TYPE(10), INPUT_TYPE(20)},
        {INPUT_TYPE(100), INPUT_TYPE(100)},
        {INPUT_TYPE(12345), INPUT_TYPE(67890)},
        {INPUT_TYPE(0xFFFF), INPUT_TYPE(0xFFFF)},                // 65535 * 65535
        {INPUT_TYPE(0xFFFFFFFF), INPUT_TYPE(2)},                 // 2^32 - 1 times 2
        {INPUT_TYPE(1000000), INPUT_TYPE(1000000)},              // 10^6 * 10^6
        {INPUT_TYPE(6000000000ULL), INPUT_TYPE(5000000000ULL)},  // 6*10^9 * 5*10^9 = 3*10^19
        {INPUT_TYPE(0x123456789ABCDEFULL), INPUT_TYPE(2)},       // Large number times 2
    };

    int numTests = sizeof(tests) / sizeof(tests[0]);
    int passed = 0;

    for (int i = 0; i < numTests; i++) {
        INPUT_TYPE a = tests[i].a;
        INPUT_TYPE b = tests[i].b;

        OUTPUT_TYPE result = R.Multiply(a, b);

        // Compute expected result using OUTPUT_TYPE arithmetic
        OUTPUT_TYPE expected = OUTPUT_TYPE(a) * OUTPUT_TYPE(b);

        bool correct = (result == expected);

        if (correct) {
            passed++;
            std::cout << "  PASS: " << a << " * " << b << " = "
                      << result << std::endl;
        } else {
            std::cout << "  FAIL: " << a << " * " << b << std::endl;
            std::cout << "    Got:      " << result << std::endl;
            std::cout << "    Expected: " << expected << std::endl;
        }
    }

    std::cout << std::endl;
    std::cout << "Passed " << passed << "/" << numTests << " tests" << std::endl;
    std::cout << std::endl;
}


int main(int argc, char* argv[])
{
    std::cout << "=== ADD-based Multiplication Verification ===" << std::endl;
    std::cout << "Configured for " << NUM_BITS << "-bit operands, "
              << numberOfMultRelations << " moduli" << std::endl;
    std::cout << std::endl;

    // Bit-width from compile-time configuration
    int n = NUM_BITS;

    // Parse ALL arguments BEFORE creating CUDD manager, because cache size
    // can only be set at construction time.
    bool enableReorder = false;
    bool bigCache = false;
    bool noGC = false;
    unsigned int explicitCacheSize = 0;
    std::string testName;
    unsigned int testK = 0;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--verbose" || arg == "-v") {
            verbose = true;
        } else if (arg == "--reorder" || arg == "-r") {
            enableReorder = true;
        } else if (arg == "--big-cache") {
            bigCache = true;
        } else if (arg == "--no-gc") {
            noGC = true;
        } else if (arg.substr(0, 14) == "--cache-slots=") {
            explicitCacheSize = (unsigned int)atol(arg.substr(14).c_str());
        } else if (testName.empty()) {
            testName = arg;
        } else if (testK == 0 && isdigit(arg[0])) {
            testK = (unsigned int)atoi(arg.c_str());
        }
    }

    bool isKaratsubaTest = (testName == "karatsuba" || testName == "karatsuba-all");

    // =========================================================================
    // Compute initial cache size for CUDD's computed table.
    //
    // CUDD's computed table is direct-mapped: collisions EVICT cached results.
    // The hash (ddCHash2) depends on pointer addresses, so ASLR causes
    // non-deterministic collision patterns — some runs fast, some hang.
    //
    // For non-Karatsuba Apply (e.g., shift-and-add): the per-Apply working
    // set is O(k^2 * n).  We use 2x headroom: k^2 * 2n.
    //
    // For Karatsuba: the (Z1+Z2)+Z0 step creates up to O(k^4 * n/2)
    // subproblems (k^4 per ply in the n/2 worst plys), so we need a much
    // larger cache.
    //
    // IMPORTANT: don't over-allocate the cache!  A CUDD cache that exceeds
    // the CPU's L3 cache causes every hash lookup to miss in hardware,
    // which can slow things down by 5-10x.
    //
    // Auto-resize: DD_MAX_CACHE_TO_SLOTS_RATIO=1024 (in cuddInt.h) decouples
    // cache from unique table.  We only enable aggressive auto-resize
    // (minHit=2) when the initial cache was capped below our estimate,
    // i.e., when the cache may be insufficient.  Otherwise we leave CUDD's
    // default minHit=30, which effectively prevents spurious resizing.
    // =========================================================================
    unsigned int cacheSize;
    bool cacheCapped = false;
    if (explicitCacheSize > 0) {
        // User-specified cache size (must be power of 2)
        cacheSize = explicitCacheSize;
    } else {
        // Determine the largest modulus k to size for.
        unsigned int kEffective;
        if (testK > 0) {
            kEffective = testK;
        } else {
            kEffective = AllModuli[numberOfMultRelations - 1];
        }

        unsigned long long estimated;
        if (isKaratsubaTest) {
            // Karatsuba (Z1+Z2)+Z0 step: O(k^4 * n/2) subproblems.
            // Size the cache to hold these with load factor < 1.
            estimated = (unsigned long long)(n / 2)
                        * kEffective * kEffective * kEffective * kEffective;
        } else {
            // Per-Apply working set is O(k^2 * n), with 2x headroom.
            estimated = (unsigned long long)kEffective * kEffective * (2 * n);
        }

        // Round up to next power of 2, clamped to [CUDD_CACHE_SLOTS, maxDefault].
        // Default ceiling 2^26 = 64M slots (2 GB); --big-cache raises to
        // 2^28 = 256M slots (8 GB).
        unsigned int maxDefault = bigCache ? (1u << 28) : (1u << 26);
        cacheSize = CUDD_CACHE_SLOTS;  // minimum = CUDD default (262144)
        while (cacheSize < estimated && cacheSize < maxDefault) {
            cacheSize <<= 1;
        }

        cacheCapped = (estimated > cacheSize);
        if (cacheCapped) {
            std::cout << "  NOTE: estimated cache need (" << estimated
                      << " slots) exceeds cap (" << maxDefault << ")."
                      << std::endl;
            if (!bigCache) {
                std::cout << "  Use --big-cache to raise the cap, or "
                          << "--cache-slots=N for explicit control." << std::endl;
            }
        }
    }

    // When auto-resize is needed (cacheCapped), scale the unique table's
    // initial slots-per-subtable so that cacheSlack >= 0 from the start.
    // Formula: DD_MAX_CACHE_TO_SLOTS_RATIO * totalSlots >= 2 * cacheSize
    //          totalSlots = (2n + 1) * numSlots
    //          numSlots >= 2 * cacheSize / (ratio * (2n + 1))
    // When the cache is already adequate, the default 256 is fine.
    unsigned int numSubtables = 2 * n + 1;  // variables + constant
    unsigned int numUniqueSlots = CUDD_UNIQUE_SLOTS;  // default 256
    if (cacheCapped) {
        unsigned long long needed = 2ULL * cacheSize
                                    / ((unsigned long long)DD_MAX_CACHE_TO_SLOTS_RATIO * numSubtables);
        if (needed > numUniqueSlots) {
            // Round up to next power of 2
            unsigned int s = CUDD_UNIQUE_SLOTS;
            while (s < needed) s <<= 1;
            numUniqueSlots = s;
        }
    }

    Cudd mgr(0, 0, numUniqueSlots, cacheSize, 0);
    DdManager* ddMgr = mgr.getManager();

    if (cacheCapped) {
        // Cache was capped below our estimate — enable aggressive auto-resize.
        // Lower minHit from 30 to 2: cache grows when hit rate > 66.7%.
        // Reset cacheMisses to match, so the threshold is reachable.
        Cudd_SetMinHit(ddMgr, 2);
        ddMgr->cacheMisses = (double)((int)(cacheSize * 2 + 1));
    } else {
        // Cache is already adequate.  Cap maxCacheHard to prevent growth.
        // This also forces cacheSlack negative (= maxCacheHard - 2*cacheSlots
        // < 0), which short-circuits the resize check in cuddCacheLookup2's
        // hot path, avoiding a double-precision multiply on every cache hit.
        // (With DD_MAX_CACHE_TO_SLOTS_RATIO=1024, cacheSlack would otherwise
        // become positive once variables are created, adding measurable overhead.)
        ddMgr->maxCacheHard = cacheSize;
    }

    unsigned int minHit = (unsigned int)Cudd_ReadMinHit(ddMgr);
    std::cout << "  Cache: " << cacheSize << " slots ("
              << (cacheSize * sizeof(DdCache) / (1024*1024)) << " MB)"
              << ", unique table: " << numUniqueSlots << " slots/subtable"
              << std::endl;
    std::cout << "  Max cache hard: " << Cudd_ReadMaxCacheHard(ddMgr)
              << ", cacheSlack: " << ddMgr->cacheSlack
              << ", auto-resize: minHit=" << minHit
              << (cacheCapped ? " (aggressive)" : "") << std::endl;

    if (enableReorder) {
        std::cout << "  Dynamic variable reordering: ENABLED (CUDD_REORDER_SIFT)" << std::endl;
        mgr.AutodynEnable(CUDD_REORDER_SIFT);
    } else {
        std::cout << "  Dynamic variable reordering: DISABLED" << std::endl;
    }

    if (noGC) {
        Cudd_DisableGarbageCollection(ddMgr);
        std::cout << "  Garbage collection: DISABLED" << std::endl;
    }

    std::cout << std::endl;

    if (testName.empty()) {
        std::cout << "Usage: " << argv[0] << " <test> [options]" << std::endl;
        std::cout << std::endl;
        std::cout << "Tests:" << std::endl;
        std::cout << "  basic <k>          - Test MultModK with modulus k" << std::endl;
        std::cout << "  crt                - Test MultRelation with CRT" << std::endl;
        std::cout << "  shiftadd <k>       - Verify shift-and-add for modulus k" << std::endl;
        std::cout << "  shiftadd-all       - Verify shift-and-add for all moduli" << std::endl;
        std::cout << "  karatsuba <k>      - Verify Karatsuba for modulus k" << std::endl;
        std::cout << "  karatsuba-all      - Verify Karatsuba for all moduli" << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  --verbose, -v      - Print intermediate structure sizes and cache stats" << std::endl;
        std::cout << "  --reorder, -r      - Enable CUDD dynamic variable reordering" << std::endl;
        std::cout << "  --big-cache        - Raise cache size cap (for large Karatsuba)" << std::endl;
        std::cout << "  --cache-slots=N    - Set exact cache size (power of 2)" << std::endl;
        std::cout << "  --no-gc            - Disable CUDD garbage collection" << std::endl;
        std::cout << std::endl;
        return 1;
    }

    if (testName == "basic") {
        if (testK == 0) {
            std::cerr << "Usage: " << argv[0] << " basic <k> [--reorder]" << std::endl;
            return 1;
        }
        return verifyMultModK(mgr, n, testK) ? 0 : 1;
    }
    else if (testName == "crt") {
        testMultRelation();
        return 0;
    }
    else if (testName == "shiftadd") {
        if (testK == 0) {
            std::cerr << "Usage: " << argv[0] << " shiftadd <k> [--reorder]" << std::endl;
            return 1;
        }
        return ADD_VerifyShiftAndAddMultiplicationModK(mgr, n, testK) ? 0 : 1;
    }
    else if (testName == "spec-all") {
        ADD_BuildMultiplicationSpecsModuliwise(mgr, n);
        return 0;
    }
    else if (testName == "shiftadd-all") {
        return ADD_VerifyShiftAndAddMultiplicationModuliwise(mgr, n) ? 0 : 1;
    }
    else if (testName == "karatsuba") {
        if (testK == 0) {
            std::cerr << "Usage: " << argv[0] << " karatsuba <k> [--reorder]" << std::endl;
            return 1;
        }
        return ADD_VerifySubtractiveKaratsubaOneLevel(mgr, n, testK) ? 0 : 1;
    }
    else if (testName == "karatsuba-all") {
        return ADD_VerifySubtractiveKaratsubaOneLevelModuliwise(mgr, n) ? 0 : 1;
    }
    else {
        std::cerr << "Unknown test: " << testName << std::endl;
        std::cerr << "Run without arguments for usage." << std::endl;
        return 1;
    }
}
