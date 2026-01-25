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

#include <cassert>
#include <algorithm>
#include <chrono>
#include "multiplication_crt.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"

namespace CFL_OBDD {

// First 26 odd primes (Figure 12)
const unsigned int Moduli[numberOfMultRelations] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103
};

// -----------------------------------------------------------------------------
// ProtoCFLOBDDNumsModK
//
// Utility procedure for construction of a Grouping that represents numbers mod k.
// (Figure 9 in the document)
//
// Translation from 1-indexed pseudo-code to 0-indexed C++:
// - Document arrays are 1-indexed; C++ uses 0-indexed
// - Document's AReturnTuple = [1, ..., k] becomes AReturnMap = {0, 1, ..., k-1}
// - Document's BReturnTuples[i] = [i] becomes BConnection[i-1].returnMapHandle = {i-1}
// -----------------------------------------------------------------------------
CFLOBDDNodeHandle ProtoCFLOBDDNumsModK(unsigned int lev, unsigned int k) {
    // Base case: level 0 returns a Fork node
    if (lev == 0) {
        return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
    }

    // Adjust for low-level Groupings that don't have k middle vertices
    // Document line [3]: numberOfMidVertices = ((lev-1) >= logLogOfMaxModulus) ? k : min(k, 2^(2^(lev-1)))
    unsigned int numberOfMidVertices;
    if ((lev - 1) >= logLogOfMaxModulus) {
        numberOfMidVertices = k;
    } else {
        unsigned long long maxMidVerts = 1ull << (1u << (lev - 1));  // 2^(2^(lev-1))
        numberOfMidVertices = std::min((unsigned long long)k, maxMidVerts);
    }

    // Adjust for low-level Groupings that don't have k exit vertices
    // Document line [4]: numberOfExits = (lev >= logLogOfMaxModulus) ? k : min(k, 2^(2^lev))
    unsigned int numberOfExits;
    if (lev >= logLogOfMaxModulus) {
        numberOfExits = k;
    } else {
        unsigned long long maxExits = 1ull << (1u << lev);  // 2^(2^lev)
        numberOfExits = std::min((unsigned long long)k, maxExits);
    }

    // Document line [5]: AConnectionPathsModK = 2^(2^(lev-1)) mod k
    unsigned long long AConnectionPathsModK = (1ull << (1u << (lev - 1))) % k;

    // Create internal node at this level
    CFLOBDDInternalNode* g = new CFLOBDDInternalNode(lev);

    // Document line [7]: Recursive call for A-connection
    CFLOBDDNodeHandle AConn = ProtoCFLOBDDNumsModK(lev - 1, k);

    // Document line [8]: AReturnTuple = [1, ..., numberOfMidVertices]
    // 0-indexed: AReturnMap = {0, 1, ..., numberOfMidVertices-1}
    CFLOBDDReturnMapHandle AReturnMap;
    for (unsigned int i = 0; i < numberOfMidVertices; i++) {
        AReturnMap.AddToEnd(i);
    }
    AReturnMap.Canonicalize();

    g->AConnection = Connection(AConn, AReturnMap);

    // Document line [9]: numberOfBConnections = numberOfMidVertices
    g->numBConnections = numberOfMidVertices;
    g->BConnection = new Connection[numberOfMidVertices];

    // Document lines [10-11]: BConnections[1] = AConnection, BReturnTuples[1] = [1, ..., numberOfMidVertices]
    // 0-indexed: BConnection[0] uses the same AReturnMap (identity map)
    g->BConnection[0] = Connection(AConn, AReturnMap);

    // Document lines [12-17]: For i = 2 to numberOfBConnections
    // 0-indexed: for i = 1 to numberOfMidVertices-1
    for (unsigned int i = 1; i < numberOfMidVertices; i++) {
        // Document line [14]: nextVal = ((i-1) * AConnectionPathsModK) mod k
        // In 0-indexed, our i corresponds to document's (i), so we use i directly
        unsigned long long nextVal = (i * AConnectionPathsModK) % k;

        // Document lines [15-16]: BReturnTuples[i] = [nextVal + 1, ..., ((nextVal + numberOfMidVertices - 1) mod numberOfExits) + 1]
        // 0-indexed: values are nextVal, (nextVal+1) mod numberOfExits, ..., (nextVal+numberOfMidVertices-1) mod numberOfExits
        CFLOBDDReturnMapHandle BReturnMap;
        for (unsigned int j = 0; j < numberOfMidVertices; j++) {
            unsigned int val = (nextVal + j) % numberOfExits;
            BReturnMap.AddToEnd(val);
        }
        BReturnMap.Canonicalize();

        // Document line [13]: BConnections[i] = AConnection (reuse the same grouping)
        g->BConnection[i] = Connection(AConn, BReturnMap);
    }

    // Document line [18]: numberOfExits
    g->numExits = numberOfExits;

#ifdef PATH_COUNTING_ENABLED
    g->InstallPathCounts();
#endif

    // Document line [19]: return RepresentativeGrouping(g)
    return CFLOBDDNodeHandle(g);
}

// -----------------------------------------------------------------------------
// NumsModK
//
// Construction of a CFLOBDD that represents numbers mod k (either in the
// topmost Grouping's A-connection or B-connection, as specified by the
// Position argument). (Figure 8 in the document)
// -----------------------------------------------------------------------------
CFLOBDD NumsModK(unsigned int k, Position p) {
    // Document line [4]: assert(2 <= k && k <= maxModulus)
    assert(2 <= k && k <= maxModulus);

    // Document line [5]: assert(p != TopLevel)
    assert(p != TopLevel);

    // Document line [6]: assert(maxLevel >= logLogOfMaxModulus + 1)
    assert(CFLOBDDMaxLevel >= logLogOfMaxModulus + 1);

    // Document line [7]: Create internal grouping at maxLevel
    CFLOBDDInternalNode* g = new CFLOBDDInternalNode(CFLOBDDMaxLevel);

    if (p == A) {
        // Document lines [8-15]: A-connection does the work, B-connections pass through

        // Document line [9]: AConnection = ProtoCFLOBDDNumsModK(maxLevel-1, k)
        CFLOBDDNodeHandle AConn = ProtoCFLOBDDNumsModK(CFLOBDDMaxLevel - 1, k);

        // Document line [10]: AReturnTuple = [1, ..., k]
        // 0-indexed: AReturnMap = {0, 1, ..., k-1}
        CFLOBDDReturnMapHandle AReturnMap;
        for (unsigned int i = 0; i < k; i++) {
            AReturnMap.AddToEnd(i);
        }
        AReturnMap.Canonicalize();

        g->AConnection = Connection(AConn, AReturnMap);

        // Document line [11]: numberOfBConnections = k
        g->numBConnections = k;
        g->BConnection = new Connection[k];

        // Document lines [12-15]: For each B-connection, use NoDistinction with single return value
        for (unsigned int i = 0; i < k; i++) {
            // Document line [13]: BConnections[i] = NoDistinctionProtoCFLOBDD(level-1)
            // Document line [14]: BReturnTuples[i] = [i]
            // 0-indexed: BReturnMap = {i}
            CFLOBDDReturnMapHandle BReturnMap;
            BReturnMap.AddToEnd(i);
            BReturnMap.Canonicalize();

            g->BConnection[i] = Connection(
                CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel - 1],
                BReturnMap
            );
        }
    }
    else if (p == B) {
        // Document lines [17-23]: A-connection is NoDistinction, B-connection does the work

        // Document lines [18-19]: AConnection = NoDistinctionProtoCFLOBDD(level-1), AReturnTuple = [1]
        // 0-indexed: AReturnMap = {0}
        CFLOBDDReturnMapHandle AReturnMap;
        AReturnMap.AddToEnd(0);
        AReturnMap.Canonicalize();

        g->AConnection = Connection(
            CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel - 1],
            AReturnMap
        );

        // Document line [20]: numberOfBConnections = 1
        g->numBConnections = 1;
        g->BConnection = new Connection[1];

        // Document line [21]: BConnections[1] = ProtoCFLOBDDNumsModK(maxLevel-1, k)
        CFLOBDDNodeHandle BConn = ProtoCFLOBDDNumsModK(CFLOBDDMaxLevel - 1, k);

        // Document line [22]: BReturnTuples[1] = [1, ..., k]
        // 0-indexed: BReturnMap = {0, 1, ..., k-1}
        CFLOBDDReturnMapHandle BReturnMap;
        for (unsigned int i = 0; i < k; i++) {
            BReturnMap.AddToEnd(i);
        }
        BReturnMap.Canonicalize();

        g->BConnection[0] = Connection(BConn, BReturnMap);
    }

    // Document line [24]: numberOfExits = k
    g->numExits = k;

#ifdef PATH_COUNTING_ENABLED
    g->InstallPathCounts();
#endif

    // Document line [25]: return RepresentativeCFLOBDD(g, [0, ..., k-1])
    // Create top-level CFLOBDD with valueTuple = [0, 1, ..., k-1]
    CFLOBDDNodeHandle nodeHandle(g);
    ReturnMapHandle<int> valueTuple;
    for (unsigned int i = 0; i < k; i++) {
        valueTuple.AddToEnd(i);
    }
    valueTuple.Canonicalize();

    return CFLOBDD(new CFLOBDDTopNodeT<int>(nodeHandle, valueTuple));
}

// -----------------------------------------------------------------------------
// Helper function for multiplication mod k
// -----------------------------------------------------------------------------
static unsigned int currentModulus = 0;

static int MultiplyModKFunc(int a, int b) {
    return (a * b) % currentModulus;
}

// -----------------------------------------------------------------------------
// Helper function for addition mod k
// -----------------------------------------------------------------------------
static int AddModKFunc(int a, int b) {
    return (a + b) % currentModulus;
}

// -----------------------------------------------------------------------------
// MultModK
//
// Construction of a CFLOBDD that represents the multiplication relation mod k.
// The bits of the two arguments are concatenated (not interleaved).
// (Figure 10 in the document)
// -----------------------------------------------------------------------------
CFLOBDD MultModK(unsigned int k) {
    // Document line [2]: cA = NumsModK(k, A)
    CFLOBDD cA = NumsModK(k, A);

    // Document line [3]: cB = NumsModK(k, B)
    CFLOBDD cB = NumsModK(k, B);

    // Document line [4]: return BinaryApplyAndReduce(cA, cB, λx, y.((x * y) mod k))
    // We use the function pointer version of ApplyAndReduce
    // Note: ApplyAndReduce takes CFLOBDDTopNodeTRefPtr, so we pass .root and wrap result
    currentModulus = k;
    return CFLOBDD(ApplyAndReduce<int>(cA.root, cB.root, MultiplyModKFunc));
}

// -----------------------------------------------------------------------------
// MultRelation Constructor
//
// Initializes the modular multiplication relations for all 26 primes.
// (Figure 12 in the document)
// -----------------------------------------------------------------------------
MultRelation::MultRelation() {
    // Initialize ModuliArray with first 26 odd primes
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        ModuliArray[i] = Moduli[i];
    }

    // Document line [32]: assert(n < 2^(2^(maxLevel-1)))
    // With 26 primes, product > 2^133, sufficient for 64-bit multiplication
    // Requires CFLOBDDMaxLevel >= 7
    assert(CFLOBDDMaxLevel >= 7);  // TWR TEST

    // Document lines [33-35]: Build CFLOBDDs for each prime
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        ModularMultRelations[i] = MultModK(ModuliArray[i]);
    }
}

// -----------------------------------------------------------------------------
// MultRelation operator==
//
// Compares two MultRelation objects for equality by comparing each entry
// in both ModuliArray and ModularMultRelations arrays.
// -----------------------------------------------------------------------------
bool MultRelation::operator==(const MultRelation& other) const {
    // Compare each entry in ModuliArray
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        if (ModuliArray[i] != other.ModuliArray[i]) {
            return false;
        }
    }

    // Compare each entry in ModularMultRelations
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        if (!(ModularMultRelations[i] == other.ModularMultRelations[i])) {
            return false;
        }
    }

    return true;
}

// -----------------------------------------------------------------------------
// Helper: Create a constant CFLOBDD that returns a specific value
// -----------------------------------------------------------------------------
static CFLOBDD MkConstantCFLOBDD(int value) {
    // Create a CFLOBDD with NoDistinction node that always returns 'value'
    ReturnMapHandle<int> valueTuple;
    valueTuple.AddToEnd(value);
    valueTuple.Canonicalize();

    return CFLOBDD(new CFLOBDDTopNodeT<int>(
        CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel],
        valueTuple
    ));
}

// -----------------------------------------------------------------------------
// Helper function for equality comparison
// -----------------------------------------------------------------------------
static int EqualityFunc(int a, int b) {
    return (a == b) ? 1 : 0;
}

// -----------------------------------------------------------------------------
// FactorViaCRT
//
// Operation to identify all factorizations of a number v.
// (Figure 13 in the document)
// -----------------------------------------------------------------------------
CFLOBDD FactorViaCRT(unsigned int v) {
    // Document line [2]: Create MultRelation
    MultRelation R;

    // Document line [3]: Create array of slices
    CFLOBDD slices[numberOfMultRelations];

    // Document lines [4-7]: For each modulus, constrain the result
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        // Document line [5]: P = ConstantCFLOBDD(v mod Moduli[i])
        CFLOBDD P = MkConstantCFLOBDD(v % R.ModuliArray[i]);

        // Document line [6]: slices[i] = BinaryApplyAndReduce(R.ModularMultRelations[i], P, λx, y.(x == y))
        // Note: ApplyAndReduce takes CFLOBDDTopNodeTRefPtr, so we pass .root and wrap result
        slices[i] = CFLOBDD(ApplyAndReduce<int>(R.ModularMultRelations[i].root, P.root, EqualityFunc));
    }

    // Document lines [8-11]: Conjoin all slices
    // Note: Tree reduction would be more efficient, but we use simple linear reduction here
    CFLOBDD ans = slices[0];
    for (unsigned int i = 1; i < numberOfMultRelations; i++) {
        ans = MkAnd(ans, slices[i]);
    }

    // Document line [12]: return ans
    return ans;
}

// -----------------------------------------------------------------------------
// ShiftAndAddMultiplicationModK
//
// Operation to build the multiplication relation mod k by a shift-and-add construction
// -----------------------------------------------------------------------------
CFLOBDD ShiftAndAddMultiplicationModK(unsigned int k) {
  currentModulus = k;

  // CFLOBDD interpretation of a shift-and-add multiplier
  CFLOBDD a = NumsModK(currentModulus, A);        // Set of all first inputs
  CFLOBDD b = NumsModK(currentModulus, B);        // Set of all second inputs
  CFLOBDD result = MkConstantCFLOBDD(0);
  CFLOBDD two = MkConstantCFLOBDD(2);
  for (int i = 0; i < MultRelation::numBits; i++) {
    // std::cout << "Iteration " << i << std::endl;
    CFLOBDD p = MkProjection(i);   // i = 0, ..., 63 corresponds to high-order bit to low-order bit
    CFLOBDD addend = CFLOBDD(ApplyAndReduce<int>(p.root, b.root, MultiplyModKFunc));
    CFLOBDD temp = CFLOBDD(ApplyAndReduce<int>(two.root, result.root, MultiplyModKFunc));
    result = CFLOBDD(ApplyAndReduce<int>(temp.root, addend.root, AddModKFunc));
  }

  return result;
}

// -----------------------------------------------------------------------------
// VerifyShiftAndAddMultiplicationModK
//
// Operation to check whether a shift-and-add multiplier produces the correct result
// Needs to be corrected
// -----------------------------------------------------------------------------
bool VerifyShiftAndAddMultiplicationModK(unsigned int k) {
  currentModulus = k;
  CFLOBDD result = ShiftAndAddMultiplicationModK(k);

  // Check whether result holds the correct value
  CFLOBDD specification = MultModK(currentModulus);
  std::cout << "Comparison of result with the specification of multiplication mod " << currentModulus << ": " << (result == specification) << std::endl;

  unsigned int nodeCount, edgeCount;
  unsigned int returnEdgesCount, returnEdgesObjCount;
  specification.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
  std::cout << nodeCount << ", " << edgeCount << std::endl;
  result.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
  std::cout << nodeCount << ", " << edgeCount << std::endl;

  	// Lookup a selection of products in specification and result
	// Create assignments for 128 variables from two 64-bit values
 	for(unsigned long long i = 20; i < 30; i++) {
		for(unsigned long long k = 20; k < 30; k++) {
			// Create a 128-bit assignment with i in the A position and k in the B position
			SH_OBDD::Assignment a(2*MultRelation::numBits);
			// Most significant bit first (high-order to low-order as in the CRT document)
			for (unsigned int j = 0; j < MultRelation::numBits; j++) {
				a[j] = (i >> (MultRelation::numBits - 1 - j)) & 1;
				a[MultRelation::numBits+j] = (k >> (MultRelation::numBits - 1 - j)) & 1;
			}
			// std::cout << "a: ";  a.print(std::cout);  std::cout << std::endl;
			int a_result = specification.root->Evaluate(a);
			std::cout << i << " * " << k << " specification = " << a_result << " Expected: " << i*k % 5 << std::endl;
			a_result = result.root->Evaluate(a);
			std::cout << i << " * " << k << " result = " << a_result << " Expected: " << i*k % 5  << std::endl;
		}
	}
    std::cout << std:: endl;

  return (result == specification);
}

// -----------------------------------------------------------------------------
// VerifyShiftAndAddMultiplication
//
// Operation to check whether a shift-and-add multiplier produces the correct result
// for all moduli
// -----------------------------------------------------------------------------
bool MultRelation::VerifyShiftAndAddMultiplication() {
    auto start = std::chrono::high_resolution_clock::now();

    // Create MultRelation
    MultRelation specification;

    MultRelation result;
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        std::cout << "ModularMultRelations[" << i << "]" << std::endl;
        result.ModularMultRelations[i] = ShiftAndAddMultiplicationModK(result.ModuliArray[i]);
    }

    // Compare result and specification
    bool equal = (result == specification);
    std::cout << "MultRelation::VerifyShiftAndAddMultiplication: " << equal << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "MultRelation::VerifyShiftAndAddMultiplication took " << duration.count() << " ms" << std::endl;

    // Emit the sizes of result.ModularMultRelations[]
    // If equal == false, emit the sizes of specification.ModularMultRelations[]
    unsigned int nodeCount, edgeCount;
    unsigned int returnEdgesCount, returnEdgesObjCount;
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        std::cout << "result.ModularMultRelations[" << i << "]: " << std::endl;
        result.ModularMultRelations[i].CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
        std::cout << nodeCount << ", " << edgeCount << std::endl;
        if (!equal) {
            std::cout << "specification.ModularMultRelations[" << i << "]: " << std::endl;
            specification.ModularMultRelations[i].CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
            std::cout << nodeCount << ", " << edgeCount << std::endl;
        }
    }


    return equal;
}


} // namespace CFL_OBDD
