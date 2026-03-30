/**
  @file addToCflobdd.hh

  @brief Bidirectional conversion between CUDD ADDs and CFLOBDDs.

  @details Provides:
  - ADD_to_CFLOBDD<T>: Convert a CUDD ADD to a CFLOBDD_T<T>
  - CFLOBDD_to_ADD<T>: Convert a CFLOBDD_T<T> to a CUDD ADD

  Both functions are templatized on the terminal value type T (e.g., int).
  The ADD's terminal values (doubles) are cast to/from T at the boundary.

  Requirements:
  - The ADD must have 2^k variables for some k <= CFLOBDDMaxLevel.
  - CFLOBDD subsystem must be initialized before calling ADD_to_CFLOBDD.
  - The CUDD manager must remain alive during conversion.
*/

#ifndef ADD_TO_CFLOBDD_HH_
#define ADD_TO_CFLOBDD_HH_

#include "../cudd-3.0.0/cplusplus/cuddObj.hh"
#include "../cudd-3.0.0/cudd/cuddInt.h"
#include "../CFLOBDD/cflobdd_t.h"
#include "../CFLOBDD/cflobdd_node.h"

#include <vector>
#include <cassert>
#include <cmath>
#include <boost/unordered/unordered_flat_map.hpp>

using namespace CFL_OBDD;

// =========================================================================
// Memo key type: (pointer, unsigned int) pair
//
// Used to key memoization tables on (DdNode*, topPly) for ADD->CFLOBDD
// and (CFLOBDDNode*, topVar) for CFLOBDD->ADD.
// =========================================================================

struct MemoKey {
    uintptr_t ptr;
    unsigned int val;
    bool operator==(const MemoKey& o) const { return ptr == o.ptr && val == o.val; }
};

struct MemoKeyHash {
    size_t operator()(const MemoKey& k) const {
        size_t h = std::hash<uintptr_t>()(k.ptr);
        h ^= std::hash<unsigned int>()(k.val) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

inline MemoKey makeMemoKey(const void* ptr, unsigned int val) {
    return { reinterpret_cast<uintptr_t>(ptr), val };
}

// =========================================================================
// Shared exit-terminal buffer for ADD->CFLOBDD conversion
//
// Instead of allocating a separate std::vector<DdNode*> for each memo
// entry's exit-to-terminal mapping, all DdNode* pointers are stored
// contiguously in one shared buffer.  Each ADDConvertResult records
// an (offset, count) range into this buffer.  This eliminates millions
// of small heap allocations and their destructor overhead.
// =========================================================================

struct ExitTerminalBuffer {
    std::vector<DdNode*> data;

    // Append terminals and return the starting offset.
    unsigned int append(const std::vector<DdNode*>& terminals) {
        unsigned int start = data.size();
        data.insert(data.end(), terminals.begin(), terminals.end());
        return start;
    }

    // Append a single terminal and return the starting offset.
    unsigned int appendOne(DdNode* t) {
        unsigned int start = data.size();
        data.push_back(t);
        return start;
    }

    // Access terminal j of a result with the given start offset.
    DdNode* get(unsigned int start, unsigned int j) const {
        return data[start + j];
    }
};

// =========================================================================
// Result type for the inner ADD-to-CFLOBDD conversion
//
// Each result stores a CFLOBDD node handle and a range [exitStart,
// exitStart+exitCount) into a shared ExitTerminalBuffer.  The range
// maps exit index i to the DdNode* at buffer[exitStart + i].
// These DdNode* pointers may be real ADD terminals or "virtual terminals"
// (ADD internal nodes at a half-height boundary).
// =========================================================================

struct ADDConvertResult {
    CFLOBDDNodeHandle nodeHandle;
    unsigned int exitStart;   // offset into shared ExitTerminalBuffer
    unsigned int exitCount;   // number of exits
};

// Memo table type for ADD->CFLOBDD (uses boost::unordered_flat_map for
// cache-friendly open addressing and O(1) destruction).
typedef boost::unordered_flat_map<MemoKey, ADDConvertResult, MemoKeyHash>
    ADDConvertMemo;

// =========================================================================
// ADD_to_CFLOBDD: Convert a CUDD ADD to a CFLOBDD
// =========================================================================

// Inner recursive function (not part of public API)
ADDConvertResult ADDConvertInner(
    DdManager* mgr,
    DdNode* node,
    unsigned int level,
    unsigned int topPly,
    unsigned int bottomPly,
    ADDConvertMemo& memo,
    ExitTerminalBuffer& exitBuf
);

// Outer wrapper
template <typename T>
CFLOBDD_T<T> ADD_to_CFLOBDD(Cudd &mgr, const ADD &f);

// =========================================================================
// CFLOBDD_to_ADD: Convert a CFLOBDD to a CUDD ADD
// =========================================================================

// Proto-ADD memo entry: a proto-ADD (with exit-index terminals) and the
// topVar at which it was originally built.
struct ProtoADDEntry {
    ADD protoADD;
    unsigned int topVar;
};

// Proto-ADD memo table: keyed on CFLOBDDNode* (proto-CFLOBDD identity).
typedef std::unordered_map<CFLOBDDNode*, ProtoADDEntry> ProtoADDMemo;

// Fused variable-shift and leaf-splice: traverses a proto-ADD (whose
// terminals are exit indices 0, 1, ...) in a single pass, simultaneously
// shifting variable indices by `offset` and replacing each exit-index
// terminal i with leaves[i].  Uses its own memo (keyed on DdNode*) to
// handle DAG sharing within the proto-ADD.
ADD FusedVariableShiftAndSplice(
    Cudd &mgr,
    DdNode* protoNode,
    int offset,
    const std::vector<ADD> &leaves,
    std::unordered_map<DdNode*, ADD> &spliceMemo
);

// Inner recursive function (not part of public API).
// Converts a CFLOBDD node to a CUDD ADD bottom-up: B-connections (lower
// variables) are converted first, then their ADDs become the leaves for
// the A-connection conversion.  Return maps are consumed by indexing into
// the leaves vector, so no post-hoc substitution is needed.
//
// Proto-ADD memoization: the first time a CFLOBDD node is converted, a
// proto-ADD (with exit-index terminals) is built and cached in protoMemo.
// On subsequent encounters of the same node at a different topVar, the
// cached proto-ADD is reused via FusedVariableShiftAndSplice.
ADD CFLOBDDConvertNodeToADD(
    Cudd &mgr,
    CFLOBDDNodeHandle nodeHandle,
    unsigned int level,
    unsigned int topVar,
    const std::vector<ADD> &leaves,
    ProtoADDMemo &protoMemo
);

// Outer wrapper
template <typename T>
ADD CFLOBDD_to_ADD(Cudd &mgr, const CFLOBDD_T<T> &f);

// =========================================================================
// Template implementations (must be in header)
// =========================================================================

template <typename T>
CFLOBDD_T<T> ADD_to_CFLOBDD(Cudd &mgr, const ADD &f) {
    DdManager* ddMgr = mgr.getManager();
    int numVars = Cudd_ReadSize(ddMgr);

    // Check numVars is a power of 2
    assert(numVars > 0 && (numVars & (numVars - 1)) == 0
           && "Number of ADD variables must be a power of 2");

    unsigned int level = 0;
    { int nv = numVars; while (nv > 1) { nv >>= 1; level++; } }

    assert(level <= CFLOBDDMaxLevel
           && "ADD requires more levels than CFLOBDDMaxLevel");

    // The memo table uses boost::unordered_flat_map for cache-friendly
    // open addressing and fast O(1) destruction (no per-entry deallocation).
    ADDConvertMemo memo;

    // All exit-terminal DdNode* pointers are stored contiguously in this
    // shared buffer, avoiding per-entry vector allocations.  Each memo
    // entry records a (start, count) range into this buffer.
    ExitTerminalBuffer exitBuf;

    ADDConvertResult result = ADDConvertInner(
        ddMgr, f.getNode(), level, 0, numVars, memo, exitBuf);

    // Build top-level return map: exit index -> T value
    ReturnMapHandle<T> topReturnMap;
    for (unsigned int i = 0; i < result.exitCount; i++) {
        topReturnMap.AddToEnd(
            static_cast<T>(Cudd_V(exitBuf.get(result.exitStart, i))));
    }
    topReturnMap.Canonicalize();

    // Wrap in CFLOBDDTopNodeT<T>
    auto topNode = new CFLOBDDTopNodeT<T>(result.nodeHandle, topReturnMap);
    return CFLOBDD_T<T>(topNode);
}

template <typename T>
ADD CFLOBDD_to_ADD(Cudd &mgr, const CFLOBDD_T<T> &f) {
    unsigned int level = f.root->rootConnection.entryPointHandle.handleContents->Level();
    unsigned int numVars = 1u << level;

    // Ensure CUDD manager has enough variables
    while ((unsigned int)Cudd_ReadSize(mgr.getManager()) < numVars) {
        mgr.addVar();
    }

    // Build leaf ADDs from the top-level return map.
    // Each exit index i maps to a constant ADD with value T(returnMap[i]).
    unsigned int numExits = f.root->rootConnection.entryPointHandle.handleContents->numExits;
    CFLOBDDReturnMapHandle topRM = f.root->rootConnection.returnMapHandle;
    std::vector<ADD> leaves(numExits);
    for (unsigned int i = 0; i < numExits; i++) {
        leaves[i] = mgr.constant(static_cast<double>(topRM.Lookup(i)));
    }

    ProtoADDMemo protoMemo;
    return CFLOBDDConvertNodeToADD(
        mgr, f.root->rootConnection.entryPointHandle, level, 0, leaves, protoMemo);
}

#endif // ADD_TO_CFLOBDD_HH_
