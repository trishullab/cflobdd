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
#include "cross_product.h"
#include <boost/multiprecision/cpp_int.hpp>

namespace CFL_OBDD {

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

// -----------------------------------------------------------------------------
// PrintSize
//
// Utility procedure to print the number of groupings and edges
// -----------------------------------------------------------------------------
void PrintSize(CFLOBDD C) {
    unsigned int nodeCount, edgeCount;
    unsigned int returnEdgesCount, returnEdgesObjCount;
    C.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
    std::cout << nodeCount << ", " << edgeCount << std::endl;
}

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
    // Base case: at bottomLevel, return a "lifted fork" node
    // MkDistinction(bottomLevel, 0) is a level-bottomLevel node with 2 exits,
    // semantically equivalent to ForkNode but padded with NoDistinctionNodes below.
    if (lev == bottomLevel) {
        if (bottomLevel == 0)
            return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
        else
            return MkDistinction(bottomLevel, 0);
    }

    // virtualLevel = level in the non-embedded structure (pairs with virtualMaxLevel)
    // Used for modular arithmetic (how many bits this subtree handles: 2^virtualLevel)
    unsigned int virtualLevel = lev - bottomLevel;

    // Adjust for low-level Groupings that don't have k middle vertices
    // Document line [3]: numberOfMidVertices = ((virtualLevel-1) >= logLogOfMaxModulus) ? k : min(k, 2^(2^(virtualLevel-1)))
    unsigned int numberOfMidVertices;
    if ((virtualLevel - 1) >= logLogOfMaxModulus) {
        numberOfMidVertices = k;
    } else {
        INPUT_TYPE maxMidVerts = INPUT_TYPE(1) << (1u << (virtualLevel - 1));  // 2^(2^(virtualLevel-1))
        numberOfMidVertices = static_cast<unsigned int>(std::min((INPUT_TYPE)k, maxMidVerts));
    }

    // Adjust for low-level Groupings that don't have k exit vertices
    // Document line [4]: numberOfExits = (virtualLevel >= logLogOfMaxModulus) ? k : min(k, 2^(2^virtualLevel))
    unsigned int numberOfExits;
    if (virtualLevel >= logLogOfMaxModulus) {
        numberOfExits = k;
    } else {
        INPUT_TYPE maxExits = INPUT_TYPE(1) << (1u << virtualLevel);  // 2^(2^virtualLevel)
        numberOfExits = static_cast<unsigned int>(std::min((INPUT_TYPE)k, maxExits));
    }

    // Document line [5]: AConnectionPathsModK = 2^(2^(virtualLevel-1)) mod k
    // Computed via repeated squaring to avoid overflow
    unsigned long long AConnectionPathsModK = 2 % k;
    for (unsigned int s = 0; s < (virtualLevel - 1); s++) {
        AConnectionPathsModK = (AConnectionPathsModK * AConnectionPathsModK) % k;
    }

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
// Extended with AA, AB, BA, and BB positions to support Karatsuba multiplication
// -----------------------------------------------------------------------------
CFLOBDD NumsModK(unsigned int k, Position p) {
    // Document line [4]: assert(2 <= k && k <= maxModulus)
    assert(2 <= k && k <= maxModulus);

    // Document line [5]: assert(p != TopLevel)
    assert(p != TopLevel);

    // Document line [6]: assert(maxLevel >= logLogOfMaxModulus + 1)
    // Ensures that there are k exit vertices
    // Changed to support Karatsuba multiplication: ensures that A, B, AA, AB, BA, and BB variants all have k exit vertices
    assert(virtualMaxLevel - 2 >= logLogOfMaxModulus);

    CFLOBDDReturnMapHandle IdReturnMap;
    for (unsigned int i = 0; i < k; i++) {
        IdReturnMap.AddToEnd(i);
    }
    IdReturnMap.Canonicalize();

    // Document line [7]: Create internal grouping at maxLevel
    CFLOBDDInternalNode* g = new CFLOBDDInternalNode(CFLOBDDMaxLevel);

    if (p == A) {
        // Document lines [8-15]: A-connection does the work, B-connections pass through

        // Document line [9]: AConnection = ProtoCFLOBDDNumsModK(maxLevel-1, k)
        CFLOBDDNodeHandle AConn = ProtoCFLOBDDNumsModK(CFLOBDDMaxLevel - 1, k);

        // Document line [10]: AReturnTuple = [1, ..., k]
        // 0-indexed: AReturnMap = {0, 1, ..., k-1}
        g->AConnection = Connection(AConn, IdReturnMap);

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
        g->BConnection[0] = Connection(BConn, IdReturnMap);
    }
    else if (p == AA) {
        // Create the AConnection, with ProtoCFLOBDDNumsModK in the AA position and NoDistinction in the AB positions
        CFLOBDDInternalNode* gA = new CFLOBDDInternalNode(CFLOBDDMaxLevel-1);
        CFLOBDDNodeHandle AAConn = ProtoCFLOBDDNumsModK(CFLOBDDMaxLevel - 2, k);
        gA->AConnection = Connection(AAConn, IdReturnMap);
        gA->numBConnections = k;
        gA->BConnection = new Connection[k];
        for (unsigned int i = 0; i < k; i++) {
            CFLOBDDReturnMapHandle BReturnMap;
            BReturnMap.AddToEnd(i);
            BReturnMap.Canonicalize();

            gA->BConnection[i] = Connection(
                CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel - 2],
                BReturnMap
            );
        }       
        gA->numExits = k;
#ifdef PATH_COUNTING_ENABLED
        gA->InstallPathCounts();
#endif
    	CFLOBDDNodeHandle AConn(gA);

        g->AConnection = Connection(AConn, IdReturnMap);
        g->numBConnections = k;
        g->BConnection = new Connection[k];
        for (unsigned int i = 0; i < k; i++) {
            CFLOBDDReturnMapHandle BReturnMap;
            BReturnMap.AddToEnd(i);
            BReturnMap.Canonicalize();

            g->BConnection[i] = Connection(
                CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel - 1],
                BReturnMap
            );
        }
    }
    else if (p == AB) {
        // Create the AConnection, with NoDistinction in the AA position and ProtoCFLOBDDNumsModK in the AB position
        CFLOBDDInternalNode* gA = new CFLOBDDInternalNode(CFLOBDDMaxLevel-1);

        CFLOBDDReturnMapHandle AAReturnMap;
        AAReturnMap.AddToEnd(0);
        AAReturnMap.Canonicalize();

        gA->AConnection = Connection(
            CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel - 2],
            AAReturnMap
        );

        gA->numBConnections = 1;
        gA->BConnection = new Connection[1];
        CFLOBDDNodeHandle ABConn = ProtoCFLOBDDNumsModK(CFLOBDDMaxLevel - 2, k);
        gA->BConnection[0] = Connection(ABConn, IdReturnMap);
        gA->numExits = k;
#ifdef PATH_COUNTING_ENABLED
        gA->InstallPathCounts();
#endif
    	CFLOBDDNodeHandle AConn(gA);

        g->AConnection = Connection(AConn, IdReturnMap);
        g->numBConnections = k;
        g->BConnection = new Connection[k];
         for (unsigned int i = 0; i < k; i++) {
            CFLOBDDReturnMapHandle BReturnMap;
            BReturnMap.AddToEnd(i);
            BReturnMap.Canonicalize();

            g->BConnection[i] = Connection(
                CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel - 1],
                BReturnMap
            );
        }
    }
    else if (p == BA) {
        CFLOBDDReturnMapHandle AReturnMap;
        AReturnMap.AddToEnd(0);
        AReturnMap.Canonicalize();

        g->AConnection = Connection(
            CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel - 1],
            AReturnMap
        );
        g->numBConnections = 1;
        g->BConnection = new Connection[1];

        // Create the BConnection, with ProtoCFLOBDDNumsModK in the BA position and NoDistinction in the BB position
        CFLOBDDInternalNode* gB = new CFLOBDDInternalNode(CFLOBDDMaxLevel-1);
        CFLOBDDNodeHandle BAConn = ProtoCFLOBDDNumsModK(CFLOBDDMaxLevel - 2, k);
        gB->AConnection = Connection(BAConn, IdReturnMap);
        gB->numBConnections = k;
        gB->BConnection = new Connection[k];
        for (unsigned int i = 0; i < k; i++) {
            CFLOBDDReturnMapHandle BReturnMap;
            BReturnMap.AddToEnd(i);
            BReturnMap.Canonicalize();

            gB->BConnection[i] = Connection(
                CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel - 2],
                BReturnMap
            );
        }       
        gB->numExits = k;
#ifdef PATH_COUNTING_ENABLED
        gB->InstallPathCounts();
#endif
    	CFLOBDDNodeHandle BConn(gB);

        g->BConnection[0] = Connection(BConn, IdReturnMap);
    }
    else if (p == BB) {
        CFLOBDDReturnMapHandle AReturnMap;
        AReturnMap.AddToEnd(0);
        AReturnMap.Canonicalize();

        g->AConnection = Connection(
            CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel - 1],
            AReturnMap
        );
        g->numBConnections = 1;
        g->BConnection = new Connection[1];

        // Create the BConnection, with NoDistinction in the BA position and ProtoCFLOBDDNumsModK in the BB position
        CFLOBDDInternalNode* gB = new CFLOBDDInternalNode(CFLOBDDMaxLevel-1);
        gB->AConnection = Connection(
            CFLOBDDNodeHandle::NoDistinctionNode[CFLOBDDMaxLevel - 2],
            AReturnMap
        );
        gB->numBConnections = 1;
        gB->BConnection = new Connection[1];
        CFLOBDDNodeHandle BBConn = ProtoCFLOBDDNumsModK(CFLOBDDMaxLevel - 2, k);
        gB->BConnection[0] = Connection(BBConn, IdReturnMap);
        gB->numExits = k;
#ifdef PATH_COUNTING_ENABLED
        gB->InstallPathCounts();
#endif
    	CFLOBDDNodeHandle BConn(gB);

        g->BConnection[0] = Connection(BConn, IdReturnMap);
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
    return static_cast<int>(((long long)a * b) % currentModulus);
}

// -----------------------------------------------------------------------------
// Helper function for addition mod k
// -----------------------------------------------------------------------------
static int AddModKFunc(int a, int b) {
    return (a + b) % currentModulus;
}

// -----------------------------------------------------------------------------
// Helper function for (2*a + b) mod k — combines doubling and addition
// -----------------------------------------------------------------------------
static int TwoAPlusBModKFunc(int a, int b) {
    return (2 * a + b) % currentModulus;
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
// Initializes the modular multiplication relations for all numberOfMultRelations primes.
// (Figure 12 in the document)
// -----------------------------------------------------------------------------
MultRelation::MultRelation() {
    // Initialize ModuliArray with the appropriate set of odd primes
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        ModuliArray[i] = Moduli[i];
    }

    // Ensures that the AA, AB, BA, and BB variants of NumsModK all have k exit vertices
    assert(virtualMaxLevel - 2 >= logLogOfMaxModulus);

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
  for (int i = 0; i < MultRelation::numBits; i++) {
    // std::cout << "Iteration " << i << std::endl;
    CFLOBDD p = MkProjection(i * stride);   // i = 0, 1, ..., numBits corresponds to high-order bit to low-order bit
    CFLOBDD addend = CFLOBDD(ApplyAndReduce<int>(p.root, b.root, MultiplyModKFunc));
    result = CFLOBDD(ApplyAndReduce<int>(result.root, addend.root, TwoAPlusBModKFunc));
    // CFLOBDD two = MkConstantCFLOBDD(2);
    // CFLOBDD temp = CFLOBDD(ApplyAndReduce<int>(two.root, result.root, MultiplyModKFunc));
    // result = CFLOBDD(ApplyAndReduce<int>(temp.root, addend.root, AddModKFunc));
  }

  return result;
}

// -----------------------------------------------------------------------------
// VerifyShiftAndAddMultiplicationModK
//
// Operation to check whether a shift-and-add multiplier produces the correct result
// -----------------------------------------------------------------------------
bool VerifyShiftAndAddMultiplicationModK(unsigned int k) {
    // Create the representation of the multiplication relation for modulus k
    auto start = std::chrono::high_resolution_clock::now();
    CFLOBDD specification = MultModK(k);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "MultModK(" << k << ") took " << duration.count() << " ms" << std::endl;
    std::cout << "Size of the multiplication relation for modulus " << k << std::endl;
    PrintSize(specification);

#ifdef TEMPORARY
    // Create the representation of the multiplication relation for modulus k via shift-and-add
    currentModulus = k;
    CFLOBDD result = ShiftAndAddMultiplicationModK(k);

    // Check whether result holds the correct value
    bool resultOfComparison = (result == specification);
    std::cout << "Comparison of result with the specification of multiplication mod " << currentModulus << ": " << resultOfComparison << std::endl;

    if (resultOfComparison == false) {
        std::cout << "Size of the multiplication relation for modulus " << k << "via shift-and-add" << std::endl;
        PrintSize(specification);
    }

  	// Lookup a selection of products in specification and result
 	for(INPUT_TYPE i = 20; i < 30; i++) {
		for(INPUT_TYPE m = 20; m < 30; m++) {
			SH_OBDD::Assignment a = MultRelation::MakeAssignment(i, m);
			int a_result = specification.root->Evaluate(a);
			std::cout << i << " * " << m << " specification = " << a_result << " Expected: " << i*m % k << std::endl;
			a_result = result.root->Evaluate(a);
			std::cout << i << " * " << m << " result = " << a_result << " Expected: " << i*m % k << std::endl;
		}
	}
    std::cout << std:: endl;
#else
    bool resultOfComparison = true;
#endif

  return resultOfComparison;
}

// -----------------------------------------------------------------------------
// BuildMultiplicationSpecsModuliwise
//
// Build the specification CFLOBDDs for all moduli and report sizes and timing
// -----------------------------------------------------------------------------
void BuildMultiplicationSpecsModuliwise() {
    auto totalStart = std::chrono::high_resolution_clock::now();

    std::cout << "Sizes of the specification's CFLOBDDs" << std::endl;
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        std::cout << "Size of the multiplication relation for the " << i+1 << "th odd prime: " << Moduli[i] << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        CFLOBDD curSpec = MultModK(Moduli[i]);
        auto end = std::chrono::high_resolution_clock::now();
        PrintSize(curSpec);
        if (i == numberOfMultRelations - 1) {
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "MultModK(" << Moduli[i] << ") took " << duration.count() << " ms" << std::endl;
        }
    }

    auto totalEnd = std::chrono::high_resolution_clock::now();
    auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(totalEnd - totalStart);
    std::cout << "BuildMultiplicationSpecsModuliwise took " << totalDuration.count() << " ms" << std::endl;
}

// -----------------------------------------------------------------------------
// VerifyShiftAndAddMultiplicationModuliwise
//
// Iterate through moduli, and check whether a shift-and-add multiplier produces
// the correct result for all moduli
// -----------------------------------------------------------------------------
bool VerifyShiftAndAddMultiplicationModuliwise() {
    auto start = std::chrono::high_resolution_clock::now();

    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        std::cout << "Testing multiplication modulo the " << i+1 << "th odd prime: " << Moduli[i] << std::endl;
        CFLOBDD curSpec = MultModK(Moduli[i]);
        CFLOBDD curShiftAndAddResult = ShiftAndAddMultiplicationModK(Moduli[i]);
        bool equal = (curSpec == curShiftAndAddResult);
        if (!equal) {
            std::cout << "Failure" << std::endl;
            return false;
        }
        // std::cout << "reduceCache size: " << CFLOBDDNodeHandle::ReduceCacheSize() << std::endl;
        // std::cout << "pairProductCache size: " << PairProductCacheSize() << std::endl;
        // std::cout << "tripleProductCache size: " << TripleProductCacheSize() << std::endl;
    }
    std::cout << "Success" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "VerifyShiftAndAddMultiplicationModuliwise took " << duration.count() << " ms" << std::endl;

    return true;
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
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        std::cout << "result.ModularMultRelations[" << i << "]: " << std::endl;
        PrintSize(result.ModularMultRelations[i]);
        if (!equal) {
            std::cout << "specification.ModularMultRelations[" << i << "]: " << std::endl;
            PrintSize(specification.ModularMultRelations[i]);
        }
    }

    return equal;
}

int SubtractModKFunc(int a_val, int b_val) {
    return (a_val + currentModulus - b_val) % currentModulus;
}

CFLOBDD SubtractViaTerminalNegation(CFLOBDD A, CFLOBDD B, unsigned int k) {
    currentModulus = k;
    return CFLOBDD(ApplyAndReduce<int>(A.root, B.root, SubtractModKFunc));
}

// -----------------------------------------------------------------------------
// SubtractiveKaratsubaOneLevel
//
// Create the CFLOBDD for the multiplication relation mod k by a one-level
// subtractive Karatsuba multiplier produces
// -----------------------------------------------------------------------------
CFLOBDD SubtractiveKaratsubaOneLevel(unsigned int k) {
    unsigned int m = MultRelation::numBits / 2;  // half of the number of bits in an input number

    // Create constant CFLOBDDs for const_m = \x.2^m mod k and const_2m = \x.2^{2*m} mod k
    // for multiplying a CFLOBDD by these powers of 2 mod k
    unsigned int shift_m = 1;
    for (unsigned int i = 0; i < m; i++) {
        shift_m = (shift_m * 2) % k;
    }
    CFLOBDD const_m = MkConstantCFLOBDD(shift_m);
    unsigned int shift_2m = (shift_m * shift_m) % k;
    CFLOBDD const_2m = MkConstantCFLOBDD(shift_2m);
    // std::cout << "shift_m = " << shift_m << "; shift_2m = " << shift_2m << std::endl;

    // Create representations of the sets of half-size numbers
    CFLOBDD X1 = NumsModK(k, AA);      // High bits of first input X
    CFLOBDD X0 = NumsModK(k, AB);      // Low bits of first input X
    CFLOBDD Y1 = NumsModK(k, BA);      // High bits of second input Y
    CFLOBDD Y0 = NumsModK(k, BB);      // Low bits of second input Y

// #define DebugSizes 1
#ifdef DebugSizes
    std::cout << "Sizes of X1, X0, Y1, and Y0" << std::endl;
    PrintSize(X1);
    PrintSize(X0);
    PrintSize(Y1);
    PrintSize(Y0);
#endif

#ifdef DebugArguments
    // For debugging
    CFLOBDD X = NumsModK(k, A);        // Set of all first inputs
    CFLOBDD Y = NumsModK(k, B);        // Set of all second inputs
    CFLOBDD X01 = CFLOBDD(ApplyAndReduce<int>(X1.root, const_m.root, MultiplyModKFunc));
    X01 = CFLOBDD(ApplyAndReduce<int>(X01.root, X0.root, AddModKFunc));
    CFLOBDD Y01 = CFLOBDD(ApplyAndReduce<int>(Y1.root, const_m.root, MultiplyModKFunc));
    Y01 = CFLOBDD(ApplyAndReduce<int>(Y01.root, Y0.root, AddModKFunc)); 
    std::cout << "X == X01: " << (X == X01) << std::endl;
    std::cout << "Y == Y01: " << (Y == Y01) << std::endl;
#endif
    
    // Compute three products: Z0 = X0 * Y0; Z2 = X1 * Y1; Z1 = (X1 - X0)(Y0 - Y1)
    CFLOBDD X_diff = SubtractViaTerminalNegation(X1, X0, k);
    CFLOBDD Y_diff = SubtractViaTerminalNegation(Y0, Y1, k);  // Y0 - Y1 !
    currentModulus = k;
    CFLOBDD Z0 = CFLOBDD(ApplyAndReduce<int>(X0.root, Y0.root, MultiplyModKFunc));
    CFLOBDD Z2 = CFLOBDD(ApplyAndReduce<int>(X1.root, Y1.root, MultiplyModKFunc));
    CFLOBDD Z1 = CFLOBDD(ApplyAndReduce<int>(X_diff.root, Y_diff.root, MultiplyModKFunc));

#ifdef DebugSizes
    std::cout << "Sizes of X_diff, Y_diff, Z0, Z2, Z1" << std::endl;
    PrintSize(X_diff);
    PrintSize(Y_diff);
    PrintSize(Z0);
    PrintSize(Z2);
    PrintSize(Z1);
#endif
    
    // Middle term: Z1 + Z2 + Z0
    CFLOBDD middle = CFLOBDD(ApplyAndReduce<int>(Z1.root, Z2.root, AddModKFunc));
#ifdef DebugSizes
    std::cout << "Size of middle = Z1 + Z2" << std::endl;
    PrintSize(middle);
#endif
    middle = CFLOBDD(ApplyAndReduce<int>(middle.root, Z0.root, AddModKFunc));
#ifdef DebugSizes
    std::cout << "Size of middle = Z1 + Z2 + Z0" << std::endl;
    PrintSize(middle);
#endif
   
    // Set karatsuba = Z2*2^(2m) + middle*2^m + Z0
    CFLOBDD term_Z2 = CFLOBDD(ApplyAndReduce<int>(
        Z2.root, const_2m.root, MultiplyModKFunc
    ));
#ifdef DebugSizes
    std::cout << "Size of term_Z2 = Z2*2^(2m)" << std::endl;
    PrintSize(term_Z2);
#endif
    CFLOBDD term_middle = CFLOBDD(ApplyAndReduce<int>(
        middle.root, const_m.root, MultiplyModKFunc
    ));
#ifdef DebugSizes
    std::cout << "Size of term_middle = middle*2^m" << std::endl;
    PrintSize(term_middle);
#endif
    CFLOBDD karatsuba = CFLOBDD(ApplyAndReduce<int>(
        term_Z2.root, term_middle.root, AddModKFunc
    ));
#ifdef DebugSizes
    std::cout << "Size of karatsuba = Z2*2^(2m) + middle*2^m" << std::endl;
    PrintSize(karatsuba);
#endif
    karatsuba = CFLOBDD(ApplyAndReduce<int>(
        karatsuba.root, Z0.root, AddModKFunc
    ));
#ifdef DebugSizes
    std::cout << "Size of karatsuba = Z2*2^(2m) + middle*2^m + Z0" << std::endl;
    PrintSize(karatsuba);
#endif

    return karatsuba;    
}

// -----------------------------------------------------------------------------
// VerifySubtractiveKaratsubaOneLevel
//
// Check whether one level of a subtractive Karatsuba multiplier produces 
// the correct result for input modulus k
// -----------------------------------------------------------------------------
bool VerifySubtractiveKaratsubaOneLevel(unsigned int k) {
    std::cout << "Verifying subtractive one-level Karatsuba mod " << k << " : ";
  auto start = std::chrono::high_resolution_clock::now();

    // Build specification
    CFLOBDD spec = MultModK(k);

    CFLOBDD karatsuba = SubtractiveKaratsubaOneLevel(k);
    
    // Compare spec and karatuba
    bool passed = (spec == karatsuba);
    
    if (passed) {
        std::cout << "  PASSED" << std::endl;
    } else {
        std::cout << "  FAILED" << std::endl;
    }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "VerifySubtractiveKaratsubaOneLevel took " << duration.count() << " ms" << std::endl;
    
    return passed;
}

// -----------------------------------------------------------------------------
// VerifySubtractiveKaratsubaOneLevelModuliwise
//
// Iterate through moduli, and check whether a one-level Karatsuba multiplier produces
// the correct result for all moduli
// -----------------------------------------------------------------------------
bool VerifySubtractiveKaratsubaOneLevelModuliwise() {
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < numberOfMultRelations; i++) {
        std::cout << "Testing multiplication modulo the " << i+1 << "th odd prime: " << Moduli[i] << std::endl;
        CFLOBDD curSpec = MultModK(Moduli[i]);
        CFLOBDD curResult = SubtractiveKaratsubaOneLevel(Moduli[i]);
        bool equal = (curSpec == curResult);
        if (!equal) {
            std::cout << "Failure" << std::endl;
            return false;
        }
    }
    std::cout << "Success" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "VerifySubtractiveKaratsubaOneLevelModuliwise took " << duration.count() << " ms" << std::endl;

    return true;    
}


// -----------------------------------------------------------------------------
// Helper function: Extended Euclidean Algorithm for computing modular inverse
// Returns the modular inverse of a modulo m
// -----------------------------------------------------------------------------
static unsigned int ModularInverse(unsigned int a, unsigned int m) {
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

// -----------------------------------------------------------------------------
// MultRelation::Multiply
//
// Multiply two m-bit values using the CRT representation and Garner's algorithm
// to reconstruct the 2*m-bit result without overflow
// -----------------------------------------------------------------------------
OUTPUT_TYPE MultRelation::Multiply(INPUT_TYPE a, INPUT_TYPE b) const {
    using namespace boost::multiprecision;

    // Create assignment with a in the A position and b in the B position
    SH_OBDD::Assignment assignment = MakeAssignment(a, b);

    // Evaluate each CFLOBDD to get remainders modulo each prime
    unsigned int remainders[numberOfMultRelations];
    for (unsigned int i = 0; i < numberOfMultRelations; i++) {
        remainders[i] = ModularMultRelations[i].root->Evaluate(assignment);
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

// -----------------------------------------------------------------------------
// MultRelation::MakeAssignment
//
// Create an assignment mapping a_val to the A-position variables and b_val to
// the B-position variables, using stride-based embedding.
// A bit j maps to variable j * stride; B bit j maps to halfVars + j * stride.
// -----------------------------------------------------------------------------
SH_OBDD::Assignment MultRelation::MakeAssignment(INPUT_TYPE a_val, INPUT_TYPE b_val) {
    SH_OBDD::Assignment assignment(totalVars);
    // Zero-initialize all variables (the don't-care bits at lower levels
    // pass through NoDistinctionNodes, so their values don't matter,
    // but we initialize for cleanliness)
    for (unsigned int j = 0; j < totalVars; j++) {
        assignment[j] = 0;
    }
    // Map input bits to their stride-based positions
    for (unsigned int j = 0; j < MultRelation::numBits; j++) {
        assignment[j * stride] = static_cast<int>((a_val >> (MultRelation::numBits - 1 - j)) & 1);
        assignment[halfVars + j * stride] = static_cast<int>((b_val >> (MultRelation::numBits - 1 - j)) & 1);
    }
    return assignment;
}

// -----------------------------------------------------------------------------
// ShiftAndAddMultiplication
//
// Multiply two m-bit values using shift-and-add multiplication
// Return the 2*m-bit product
// -----------------------------------------------------------------------------
OUTPUT_TYPE ShiftAndAddMultiplication(INPUT_TYPE a, INPUT_TYPE b) {
    OUTPUT_TYPE result = 0;
    for (int i = MultRelation::numBits - 1; i >= 0; i--) {
        INPUT_TYPE addend = (a & (INPUT_TYPE(1) << i)) ? b : 0;
        result = 2 * result + addend;
    }
    return result;
}

} // namespace CFL_OBDD
