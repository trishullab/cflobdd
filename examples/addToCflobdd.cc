/**
  @file addToCflobdd.cc

  @brief Implementation of bidirectional ADD <-> CFLOBDD conversion.
*/

#include "addToCflobdd.hh"

// =========================================================================
// ADD_to_CFLOBDD inner implementation
//
// Recursively converts an ADD sub-DAG rooted at `node` into a CFLOBDD
// node at the given `level`.  The ply range [topPly, bottomPly) defines
// which ADD variables this CFLOBDD level covers:
//   - ADD nodes with variables in [topPly, bottomPly) are internal
//   - ADD nodes at ply >= bottomPly (or constants) are treated as
//     "virtual terminals" — they become exits of this CFLOBDD node
//
// At each recursive step, the range is split at midPly = topPly + 2^(level-1):
//   - The top half [topPly, midPly) produces the A-connection
//   - The bottom half [midPly, bottomPly) produces the B-connections
//
// Memoization on (DdNode*, topPly) ensures shared ADD sub-DAGs are
// converted only once.  Exit-terminal DdNode* pointers are stored in
// the shared ExitTerminalBuffer to avoid per-entry heap allocations.
// =========================================================================

ADDConvertResult ADDConvertInner(
    DdManager* mgr,
    DdNode* node,
    unsigned int level,
    unsigned int topPly,
    unsigned int bottomPly,
    ADDConvertMemo& memo,
    ExitTerminalBuffer& exitBuf)
{
    // Check memo
    MemoKey key = makeMemoKey(node, topPly);
    auto it = memo.find(key);
    if (it != memo.end()) {
        return it->second;
    }

    ADDConvertResult result;

    // --- BASE CASE: level == 0, covers one variable at ply topPly ---
    if (level == 0) {
        assert(bottomPly == topPly + 1);

        if (cuddIsConstant(node)) {
            // Terminal node — variable at topPly is irrelevant
            result.nodeHandle = CFLOBDDNodeHandle::NoDistinctionNode[0];
            result.exitStart = exitBuf.appendOne(node);
            result.exitCount = 1;
        }
        else {
            int nodePly = Cudd_ReadPerm(mgr, Cudd_NodeReadIndex(node));
            if ((unsigned int)nodePly == topPly) {
                // Node tests this exact variable
                DdNode* T = cuddT(node);
                DdNode* E = cuddE(node);
                if (T == E) {
                    // Both children identical — variable irrelevant
                    result.nodeHandle = CFLOBDDNodeHandle::NoDistinctionNode[0];
                    result.exitStart = exitBuf.appendOne(T);
                    result.exitCount = 1;
                }
                else {
                    // ForkNode: exit 0 = else (var=0), exit 1 = then (var=1)
                    result.nodeHandle = CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
                    std::vector<DdNode*> exits = { E, T };
                    result.exitStart = exitBuf.append(exits);
                    result.exitCount = 2;
                }
            }
            else {
                // Node's variable is below topPly (ply-skipping):
                // this variable is not tested — pass node through as virtual terminal
                result.nodeHandle = CFLOBDDNodeHandle::NoDistinctionNode[0];
                result.exitStart = exitBuf.appendOne(node);
                result.exitCount = 1;
            }
        }

        memo[key] = result;
        return result;
    }

    // --- RECURSIVE CASE: level > 0 ---
    unsigned int midPly = topPly + (1u << (level - 1));

    // Step 1: Convert top half (A-connection).
    // The ply range [topPly, midPly) means any ADD node at ply >= midPly
    // is treated as a virtual terminal — these are the "middle nodes".
    ADDConvertResult AResult = ADDConvertInner(
        mgr, node, level - 1, topPly, midPly, memo, exitBuf);
    unsigned int m = AResult.exitCount;  // number of middle nodes

    // Step 2: Convert bottom halves (B-connections).
    // For each middle node, recursively convert the sub-ADD rooted at it
    // over the ply range [midPly, bottomPly).
    std::vector<ADDConvertResult> BResults(m);

    // Collect the global set of distinct terminals reachable from all
    // B-connections.  These become the exits of the level-k CFLOBDD node.
    std::vector<DdNode*> globalTerminals;
    std::unordered_map<DdNode*, unsigned int> terminalToIndex;

    for (unsigned int i = 0; i < m; i++) {
        DdNode* middleNode = exitBuf.get(AResult.exitStart, i);
        BResults[i] = ADDConvertInner(
            mgr, middleNode, level - 1, midPly, bottomPly, memo, exitBuf);

        // Add B_i's terminals to the global set (deduplicating)
        for (unsigned int j = 0; j < BResults[i].exitCount; j++) {
            DdNode* t = exitBuf.get(BResults[i].exitStart, j);
            if (terminalToIndex.find(t) == terminalToIndex.end()) {
                terminalToIndex[t] = globalTerminals.size();
                globalTerminals.push_back(t);
            }
        }
    }

    unsigned int numGlobalExits = globalTerminals.size();

    // Step 3: Build the CFLOBDD node.
    CFLOBDDInternalNode* N = new CFLOBDDInternalNode(level);

    // AConnection: identity return map (A's exits correspond 1:1 to middle vertices)
    CFLOBDDReturnMapHandle AReturnMap = MakeIdentityReturnMap(m);
    N->AConnection = Connection(AResult.nodeHandle, AReturnMap);

    // BConnections: each B_i's return map maps its local exit indices
    // to the global exit numbering of this node.
    N->numBConnections = m;
    N->BConnection = new Connection[m];

    for (unsigned int i = 0; i < m; i++) {
        CFLOBDDReturnMapHandle BReturnMap;
        for (unsigned int j = 0; j < BResults[i].exitCount; j++) {
            DdNode* t = exitBuf.get(BResults[i].exitStart, j);
            BReturnMap.AddToEnd(terminalToIndex[t]);
        }
        BReturnMap.Canonicalize();
        N->BConnection[i] = Connection(BResults[i].nodeHandle, BReturnMap);
    }

    N->numExits = numGlobalExits;
#ifdef PATH_COUNTING_ENABLED
    N->InstallPathCounts();
#endif

    result.nodeHandle = CFLOBDDNodeHandle(N);  // auto-canonicalizes
    result.exitStart = exitBuf.append(globalTerminals);
    result.exitCount = numGlobalExits;

    memo[key] = result;
    return result;
}

// =========================================================================
// FusedVariableShiftAndSplice
//
// Traverses a proto-ADD (whose terminals are exit indices 0, 1, ...)
// in a single pass, simultaneously:
//   - shifting each internal node's variable index by `offset`
//   - replacing each exit-index terminal i with leaves[i]
//
// This is used when a previously-converted CFLOBDD node is encountered
// again at a different variable position (topVar).  Instead of re-doing
// the full CFLOBDD traversal, we reuse the cached proto-ADD structure
// and just adjust variable indices and splice in new leaf ADDs.
//
// Uses its own memo (keyed on DdNode*) to handle DAG sharing within
// the proto-ADD.  The offset and leaves are fixed for the entire call,
// so only the DdNode* is needed as the key.
// =========================================================================

ADD FusedVariableShiftAndSplice(
    Cudd &mgr,
    DdNode* protoNode,
    int offset,
    const std::vector<ADD> &leaves,
    std::unordered_map<DdNode*, ADD> &spliceMemo)
{
    auto it = spliceMemo.find(protoNode);
    if (it != spliceMemo.end()) {
        return it->second;
    }

    ADD result;
    if (cuddIsConstant(protoNode)) {
        // Terminal: the value is an exit index — replace with the leaf ADD
        unsigned int exitIndex = (unsigned int)cuddV(protoNode);
        assert(exitIndex < leaves.size());
        result = leaves[exitIndex];
    }
    else {
        // Internal node: shift the variable index and recurse
        int v = Cudd_NodeReadIndex(protoNode);
        ADD T = FusedVariableShiftAndSplice(mgr, cuddT(protoNode), offset, leaves, spliceMemo);
        ADD E = FusedVariableShiftAndSplice(mgr, cuddE(protoNode), offset, leaves, spliceMemo);
        result = mgr.addVar(v + offset).Ite(T, E);
    }

    spliceMemo[protoNode] = result;
    return result;
}

// =========================================================================
// CFLOBDD_to_ADD inner implementation (bottom-up with proto-ADD memoization)
//
// Converts a CFLOBDD node to a CUDD ADD by building bottom-up:
//   1. Convert each B-connection sub-node first (lower variables).
//      Each B_i's return map tells us which of this node's exits each
//      B_i exit maps to, so B_i's leaves are the parent's leaves
//      re-indexed through the return map.
//   2. Convert the A-connection sub-node (upper variables), using
//      the B-connection ADDs as its leaves.  The A return map maps
//      A's exits to middle vertex indices, selecting which B-ADD
//      each A exit leads to.
//
// Proto-ADD memoization:
//   On a cache miss, after building the actual ADD, we also build and
//   cache a "proto-ADD" — an ADD with the same structure but with
//   exit-index constants (0.0, 1.0, ...) as terminals instead of real
//   values.  On a subsequent cache hit for the same CFLOBDD node at a
//   different topVar, we reuse the proto-ADD via FusedVariableShiftAndSplice,
//   which shifts variable indices and splices in the new leaves in one pass.
//
// NoDistinctionNode handling:
//   A NoDistinctionNode at any level has exactly 1 exit.  Its
//   corresponding ADD is just leaves[0] (no variable structure at all,
//   analogous to ply-skipping in the ADD).  We short-circuit these
//   without recursion or proto-ADD construction.
// =========================================================================

ADD CFLOBDDConvertNodeToADD(
    Cudd &mgr,
    CFLOBDDNodeHandle nodeHandle,
    unsigned int level,
    unsigned int topVar,
    const std::vector<ADD> &leaves,
    ProtoADDMemo &protoMemo)
{
    CFLOBDDNode* cfNode = nodeHandle.handleContents;

    // Fast path: NoDistinctionNode at any level has 1 exit — just return
    // leaves[0].  This avoids recursing through NoDistinctionNode padding
    // levels in the topmost embedding.
    if (nodeHandle == CFLOBDDNodeHandle::NoDistinctionNode[level]) {
        return leaves[0];
    }

    // Check proto-ADD memo: have we converted this CFLOBDD node before?
    auto pit = protoMemo.find(cfNode);
    if (pit != protoMemo.end()) {
        // Hit: reuse the cached proto-ADD.  The proto-ADD was built at
        // cachedTopVar; we need it at topVar, so shift by the difference.
        int offset = (int)topVar - (int)pit->second.topVar;
        std::unordered_map<DdNode*, ADD> spliceMemo;
        return FusedVariableShiftAndSplice(
            mgr, pit->second.protoADD.getNode(), offset, leaves, spliceMemo);
    }

    // Cache miss: do the full bottom-up conversion.

    // --- BASE CASE: level == 0 ---
    if (level == 0) {
        ADD protoADD;
        ADD result;
        switch (cfNode->NodeKind()) {
            case CFLOBDD_DONTCARE:
                // 1 exit: proto is constant 0.0; result is leaves[0]
                protoADD = mgr.constant(0.0);
                result = leaves[0];
                break;
            case CFLOBDD_FORK:
                // 2 exits: variable selects between them
                // exit 0 = var is 0 (else); exit 1 = var is 1 (then)
                protoADD = mgr.addVar(topVar).Ite(mgr.constant(1.0), mgr.constant(0.0));
                result = mgr.addVar(topVar).Ite(leaves[1], leaves[0]);
                break;
            default:
                assert(false && "Unexpected node kind at level 0");
                return leaves[0];
        }
        protoMemo[cfNode] = { protoADD, topVar };
        return result;
    }

    // --- RECURSIVE CASE: level > 0 ---
    CFLOBDDInternalNode* node =
        dynamic_cast<CFLOBDDInternalNode*>(cfNode);
    assert(node != nullptr);

    unsigned int midVar = topVar + (1u << (level - 1));
    unsigned int m = node->numBConnections;
    unsigned int numExits = node->numExits;

    // Step 1: Convert each B sub-node (lower variables first).
    // Each B_i's return map maps B_i's exits to this node's exits,
    // so B_i's leaves are this node's leaves re-indexed through the map.
    std::vector<ADD> B_ADDs(m);
    for (unsigned int i = 0; i < m; i++) {
        unsigned int numBExits = node->BConnection[i].entryPointHandle.handleContents->numExits;
        std::vector<ADD> B_leaves(numBExits);
        for (unsigned int j = 0; j < numBExits; j++) {
            B_leaves[j] = leaves[node->BConnection[i].returnMapHandle.Lookup(j)];
        }
        B_ADDs[i] = CFLOBDDConvertNodeToADD(
            mgr, node->BConnection[i].entryPointHandle, level - 1, midVar, B_leaves, protoMemo);
    }

    // Step 2: Convert A sub-node (upper variables).
    // A's return map maps A's exits to middle vertex indices 0..m-1.
    // The leaves for A are the B-connection ADDs, selected through the map.
    unsigned int numAExits = node->AConnection.entryPointHandle.handleContents->numExits;
    std::vector<ADD> A_leaves(numAExits);
    for (unsigned int j = 0; j < numAExits; j++) {
        A_leaves[j] = B_ADDs[node->AConnection.returnMapHandle.Lookup(j)];
    }

    ADD result = CFLOBDDConvertNodeToADD(
        mgr, node->AConnection.entryPointHandle, level - 1, topVar, A_leaves, protoMemo);

    // --- Build and cache proto-ADD for this node ---
    // The proto-ADD has the same structure as the real ADD but uses
    // exit-index constants (0.0, 1.0, ...) as terminals.  This allows
    // reuse via FusedVariableShiftAndSplice when the same CFLOBDD node
    // is encountered at a different variable position.

    // Proto-leaves: constant ADDs for each exit index
    std::vector<ADD> exitLeaves(numExits);
    for (unsigned int i = 0; i < numExits; i++) {
        exitLeaves[i] = mgr.constant((double)i);
    }

    // Build proto B-ADDs from cached proto-ADDs of B sub-nodes.
    // NoDistinctionNodes have 1 exit and no variable structure, so their
    // proto is just the corresponding exit-index constant directly.
    std::vector<ADD> B_protos(m);
    for (unsigned int i = 0; i < m; i++) {
        CFLOBDDNodeHandle &bHandle = node->BConnection[i].entryPointHandle;
        unsigned int numBExits = bHandle.handleContents->numExits;
        std::vector<ADD> bExitLeaves(numBExits);
        for (unsigned int j = 0; j < numBExits; j++) {
            bExitLeaves[j] = exitLeaves[node->BConnection[i].returnMapHandle.Lookup(j)];
        }
        if (bHandle == CFLOBDDNodeHandle::NoDistinctionNode[level - 1]) {
            // NoDistinctionNode: 1 exit, no variables — proto is just the leaf
            B_protos[i] = bExitLeaves[0];
        }
        else {
            // Use cached proto-ADD with variable shift and exit-index splice
            assert(protoMemo.find(bHandle.handleContents) != protoMemo.end());
            auto &entry = protoMemo[bHandle.handleContents];
            int bOffset = (int)midVar - (int)entry.topVar;
            std::unordered_map<DdNode*, ADD> spliceMemo;
            B_protos[i] = FusedVariableShiftAndSplice(
                mgr, entry.protoADD.getNode(), bOffset, bExitLeaves, spliceMemo);
        }
    }

    // Build proto A-ADD similarly.
    CFLOBDDNodeHandle &aHandle = node->AConnection.entryPointHandle;
    std::vector<ADD> A_protoLeaves(numAExits);
    for (unsigned int j = 0; j < numAExits; j++) {
        A_protoLeaves[j] = B_protos[node->AConnection.returnMapHandle.Lookup(j)];
    }
    ADD protoADD;
    if (aHandle == CFLOBDDNodeHandle::NoDistinctionNode[level - 1]) {
        protoADD = A_protoLeaves[0];
    }
    else {
        assert(protoMemo.find(aHandle.handleContents) != protoMemo.end());
        auto &entry = protoMemo[aHandle.handleContents];
        int aOffset = (int)topVar - (int)entry.topVar;
        std::unordered_map<DdNode*, ADD> aSpliceMemo;
        protoADD = FusedVariableShiftAndSplice(
            mgr, entry.protoADD.getNode(), aOffset, A_protoLeaves, aSpliceMemo);
    }

    protoMemo[cfNode] = { protoADD, topVar };
    return result;
}
