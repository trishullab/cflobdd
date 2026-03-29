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
#include <unordered_map>
#include <cassert>
#include <cmath>

using namespace CFL_OBDD;

// =========================================================================
// Memo key type: (pointer, unsigned int) pair
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
// Result type for the inner ADD-to-CFLOBDD conversion
// =========================================================================

struct ADDConvertResult {
    CFLOBDDNodeHandle nodeHandle;
    std::vector<DdNode*> exitToTerminal;  // exit index -> DdNode* (real or virtual terminal)
};

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
    std::unordered_map<MemoKey, ADDConvertResult, MemoKeyHash>& memo
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

// Fused variable-shift and leaf-splice: traverses protoADD once,
// shifting variables by `offset` and replacing exit-index terminals
// with the corresponding leaf ADDs.
ADD FusedVariableShiftAndSplice(
    Cudd &mgr,
    DdNode* protoNode,
    int offset,
    const std::vector<ADD> &leaves,
    std::unordered_map<DdNode*, ADD> &spliceMemo
);

// Inner recursive function (not part of public API)
// `leaves` has length nodeHandle.numExits; leaves[i] is the ADD to place at exit i.
// `protoMemo` caches proto-ADDs keyed on CFLOBDD node identity.
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

    std::unordered_map<MemoKey, ADDConvertResult, MemoKeyHash> memo;

    ADDConvertResult result = ADDConvertInner(
        ddMgr, f.getNode(), level, 0, numVars, memo);

    // Build top-level return map: exit index -> T value
    ReturnMapHandle<T> topReturnMap;
    for (unsigned int i = 0; i < result.exitToTerminal.size(); i++) {
        topReturnMap.AddToEnd(static_cast<T>(Cudd_V(result.exitToTerminal[i])));
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

    // Build leaf ADDs from the top-level return map
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
