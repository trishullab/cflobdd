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

using namespace std::chrono;


// =============================================================================
// Global variable for the apply operator
// =============================================================================

// Global variable to store the current modulus for the apply operator.
// This is necessary because CUDD's apply interface doesn't allow passing
// additional parameters to the operator function.
// This mirrors the CFLOBDD approach in multiplication_crt.cpp which uses
// "static unsigned int currentModulus".
static int currentModulusK = 0;


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

    // Set the global modulus for the apply operator
    // (mirrors CFLOBDD's "currentModulus = k;" in MultModK)
    currentModulusK = k;

    // IMPORTANT: Flush the cache before applying with a potentially different modulus.
    // CUDD caches apply results keyed by (function pointer, operand nodes), but doesn't
    // know about our global currentModulusK. Without this flush, we could get stale
    // cached results from a previous call with a different modulus.
    DdManager* ddMgr = mgr.getManager();
    cuddCacheFlush(ddMgr);

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

    // Wrap the result in an ADD object (with proper reference counting)
    Cudd_Ref(result);
    ADD resultADD(mgr, result);
    Cudd_RecursiveDeref(ddMgr, result);

    return resultADD;
}


CUDD_VALUE_TYPE ADD_Evaluate(Cudd& mgr, const ADD& f, int n,
                              unsigned long long x, unsigned long long y)
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
        inputs[i] = (x >> bitPos) & 1;
    }

    // Fill in the y bits (variables n..2n-1)
    for (int i = 0; i < n; i++) {
        int bitPos = n - 1 - i;  // High-order to low-order
        inputs[n + i] = (y >> bitPos) & 1;
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
    unsigned long long maxVal = (1ULL << n);
    unsigned long long testCount = 0;
    unsigned long long errorCount = 0;

    // For small n, test all inputs; for large n, test a sample
    unsigned long long step = (maxVal <= 256) ? 1 : maxVal / 64;

    std::cout << "Verifying results..." << std::endl;

    start = high_resolution_clock::now();

    for (unsigned long long x = 0; x < maxVal; x += step) {
        for (unsigned long long y = 0; y < maxVal; y += step) {
            CUDD_VALUE_TYPE addResult = ADD_Evaluate(mgr, multModK, n, x, y);
            unsigned long long expected = (x * y) % k;

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
// First 26 odd primes (Figure 12 in Multiplication_via_CRT.pdf)
// =============================================================================

const unsigned int Moduli[numberOfMultRelations] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103
};


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


void MultRelation::Multiply(unsigned long long a, unsigned long long b,
                            unsigned long long& resultHigh, unsigned long long& resultLow) const
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
            u[i] = ((unsigned int)diff * inv) % ModuliArray[i];
        }
    }

    // Reconstruct the result using mixed radix representation
    // We use 128-bit arithmetic via two 64-bit values
    // For simplicity, we'll compute this using unsigned __int128 if available,
    // otherwise we use a simple approach that works for our test cases

#ifdef __SIZEOF_INT128__
    unsigned __int128 result = u[0];
    unsigned __int128 product = ModuliArray[0];

    for (unsigned int i = 1; i < numberOfMultRelations; i++) {
        result += (unsigned __int128)u[i] * product;
        product *= ModuliArray[i];
    }

    resultLow = (unsigned long long)result;
    resultHigh = (unsigned long long)(result >> 64);
#else
    // Fallback: just compute low 64 bits (sufficient for verification with small inputs)
    unsigned long long result = u[0];
    unsigned long long product = ModuliArray[0];

    for (unsigned int i = 1; i < numberOfMultRelations; i++) {
        result += (unsigned long long)u[i] * product;
        product *= ModuliArray[i];
    }

    resultLow = result;
    resultHigh = 0;  // Not computed in fallback mode
#endif
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
    std::cout << "=== Testing MultRelation (CRT-based 64-bit multiplication) ===" << std::endl;
    std::cout << std::endl;

    // Build the MultRelation with all 26 primes
    MultRelation R;

    std::cout << std::endl;
    R.PrintStats();
    std::cout << std::endl;

    // Test multiplication with some sample values
    std::cout << "--- Verifying multiplication ---" << std::endl;

    struct TestCase {
        unsigned long long a;
        unsigned long long b;
    };

    TestCase tests[] = {
        {0, 0},
        {1, 1},
        {2, 3},
        {10, 20},
        {100, 100},
        {12345, 67890},
        {0xFFFF, 0xFFFF},                // 65535 * 65535
        {0xFFFFFFFF, 2},                 // 2^32 - 1 times 2
        {1000000, 1000000},              // 10^6 * 10^6
        {6000000000ULL, 5000000000ULL},  // 6*10^9 * 5*10^9 = 3*10^19 (overflows 64 bits)
        {0x123456789ABCDEFULL, 2},       // Large number times 2
    };

    int numTests = sizeof(tests) / sizeof(tests[0]);
    int passed = 0;

    for (int i = 0; i < numTests; i++) {
        unsigned long long a = tests[i].a;
        unsigned long long b = tests[i].b;

        unsigned long long resultHigh, resultLow;
        R.Multiply(a, b, resultHigh, resultLow);

        // Compute expected result using __int128 if available
#ifdef __SIZEOF_INT128__
        unsigned __int128 expected = (unsigned __int128)a * b;
        unsigned long long expectedLow = (unsigned long long)expected;
        unsigned long long expectedHigh = (unsigned long long)(expected >> 64);

        bool correct = (resultLow == expectedLow && resultHigh == expectedHigh);
#else
        // Fallback: only check low bits for small products
        unsigned long long expectedLow = a * b;  // May overflow
        bool correct = (resultLow == expectedLow);
#endif

        if (correct) {
            passed++;
            std::cout << "  PASS: " << a << " * " << b << " = ";
            if (resultHigh > 0) {
                std::cout << "0x" << std::hex << resultHigh << std::dec << ":";
            }
            std::cout << resultLow << std::endl;
        } else {
            std::cout << "  FAIL: " << a << " * " << b << std::endl;
            std::cout << "    Got:      " << resultHigh << ":" << resultLow << std::endl;
#ifdef __SIZEOF_INT128__
            std::cout << "    Expected: " << expectedHigh << ":" << expectedLow << std::endl;
#endif
        }
    }

    std::cout << std::endl;
    std::cout << "Passed " << passed << "/" << numTests << " tests" << std::endl;
    std::cout << std::endl;
}


int main(int argc, char* argv[])
{
    std::cout << "=== ADD-based MultModK Test ===" << std::endl;
    std::cout << std::endl;

    // Check for MultRelation test mode
    if (argc >= 2 && std::string(argv[1]) == "crt") {
        testMultRelation();
        return 0;
    }

    // Create CUDD manager
    Cudd mgr(0, 0);

    // Test parameters
    int n = 4;  // Number of bits per operand
    int k = 5;  // Modulus

    if (argc >= 3) {
        n = atoi(argv[1]);
        k = atoi(argv[2]);
    }

    std::cout << "Parameters: n=" << n << " bits, k=" << k << " (modulus)" << std::endl;
    std::cout << std::endl;

    // Test NumsModK
    std::cout << "--- Testing NumsModK ---" << std::endl;
    ADD numsModK_x = ADD_NumsModK(mgr, n, k, 0);
    printADDStats(numsModK_x, "NumsModK(x)");

    ADD numsModK_y = ADD_NumsModK(mgr, n, k, n);
    printADDStats(numsModK_y, "NumsModK(y)");

    std::cout << std::endl;

    // Test MultModK
    std::cout << "--- Testing MultModK ---" << std::endl;
    bool success = verifyMultModK(mgr, n, k);

    std::cout << std::endl;

    // Test with different moduli
    if (success && n <= 8) {
        std::cout << "--- Testing with other moduli ---" << std::endl;
        int moduli[] = {3, 7, 11, 13};
        for (int m : moduli) {
            std::cout << std::endl;
            verifyMultModK(mgr, n, m);
        }
    }

    std::cout << std::endl;
    std::cout << "Usage: " << argv[0] << " [n k]     - Test MultModK with n bits and modulus k" << std::endl;
    std::cout << "       " << argv[0] << " crt       - Test full MultRelation (26 primes, 64 bits)" << std::endl;
    std::cout << std::endl;
    std::cout << "=== Test Complete ===" << std::endl;

    return success ? 0 : 1;
}
