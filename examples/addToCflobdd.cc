/**
  @file addToCflobdd.cc

  @brief Implementation of bidirectional ADD <-> CFLOBDD conversion.
*/

#include "addToCflobdd.hh"

// =========================================================================
// ADD_to_CFLOBDD inner implementation
// =========================================================================

ADDConvertResult ADDConvertInner(
    DdManager* mgr,
    DdNode* node,
    unsigned int level,
    unsigned int topPly,
    unsigned int bottomPly,
    std::unordered_map<MemoKey, ADDConvertResult, MemoKeyHash>& memo)
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
            result.exitToTerminal = { node };
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
                    result.exitToTerminal = { T };
                }
                else {
                    // ForkNode: exit 0 = else (var=0), exit 1 = then (var=1)
                    result.nodeHandle = CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
                    result.exitToTerminal = { E, T };
                }
            }
            else {
                // Node's variable is below topPly (ply-skipping)
                // This variable is not tested — pass node through as virtual terminal
                result.nodeHandle = CFLOBDDNodeHandle::NoDistinctionNode[0];
                result.exitToTerminal = { node };
            }
        }

        memo[key] = result;
        return result;
    }

    // --- RECURSIVE CASE: level > 0 ---
    unsigned int midPly = topPly + (1u << (level - 1));

    // Step 1: Convert top half (A-connection)
    // The ply range [topPly, midPly) means any node at ply >= midPly is a virtual terminal
    ADDConvertResult AResult = ADDConvertInner(mgr, node, level - 1, topPly, midPly, memo);
    unsigned int m = AResult.exitToTerminal.size();  // number of middle nodes

    // Step 2: Convert bottom halves (B-connections)
    std::vector<ADDConvertResult> BResults(m);

    // Collect global terminals (union of all B_i terminals)
    std::vector<DdNode*> globalTerminals;
    std::unordered_map<DdNode*, unsigned int> terminalToIndex;

    for (unsigned int i = 0; i < m; i++) {
        DdNode* middleNode = AResult.exitToTerminal[i];
        BResults[i] = ADDConvertInner(mgr, middleNode, level - 1, midPly, bottomPly, memo);

        for (DdNode* t : BResults[i].exitToTerminal) {
            if (terminalToIndex.find(t) == terminalToIndex.end()) {
                terminalToIndex[t] = globalTerminals.size();
                globalTerminals.push_back(t);
            }
        }
    }

    unsigned int numGlobalExits = globalTerminals.size();

    // Step 3: Build the CFLOBDD node
    CFLOBDDInternalNode* N = new CFLOBDDInternalNode(level);

    // AConnection: identity return map (A's exits map 1:1 to middle vertices)
    CFLOBDDReturnMapHandle AReturnMap = MakeIdentityReturnMap(m);
    N->AConnection = Connection(AResult.nodeHandle, AReturnMap);

    // BConnections
    N->numBConnections = m;
    N->BConnection = new Connection[m];

    for (unsigned int i = 0; i < m; i++) {
        CFLOBDDReturnMapHandle BReturnMap;
        for (unsigned int j = 0; j < BResults[i].exitToTerminal.size(); j++) {
            DdNode* t = BResults[i].exitToTerminal[j];
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
    result.exitToTerminal = globalTerminals;

    memo[key] = result;
    return result;
}

// =========================================================================
// FusedVariableShiftAndSplice
//
// Traverses a proto-ADD (whose terminals are exit indices 0, 1, ...)
// in a single pass, simultaneously shifting variable indices by `offset`
// and replacing each exit-index terminal i with leaves[i].
// Uses its own memo (keyed on DdNode*) to handle DAG sharing within
// the proto-ADD.
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
        unsigned int exitIndex = (unsigned int)cuddV(protoNode);
        assert(exitIndex < leaves.size());
        result = leaves[exitIndex];
    }
    else {
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
// On a cache miss: does the bottom-up conversion (B-connections first,
// then A-connection with B-ADDs as leaves). Also builds and caches a
// proto-ADD (with exit-index terminals) for this CFLOBDD node.
//
// On a cache hit: uses FusedVariableShiftAndSplice to reuse the cached
// proto-ADD, applying the variable offset and splicing in the new leaves
// in a single traversal.
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

    // Check proto-ADD memo
    auto pit = protoMemo.find(cfNode);
    if (pit != protoMemo.end()) {
        // Hit: reuse proto-ADD with variable shift and leaf splice
        int offset = (int)topVar - (int)pit->second.topVar;
        std::unordered_map<DdNode*, ADD> spliceMemo;
        return FusedVariableShiftAndSplice(
            mgr, pit->second.protoADD.getNode(), offset, leaves, spliceMemo);
    }

    // Miss: do bottom-up conversion and build + cache the proto-ADD

    // --- BASE CASE: level == 0 ---
    if (level == 0) {
        ADD protoADD;
        ADD result;
        switch (cfNode->NodeKind()) {
            case CFLOBDD_DONTCARE:
                protoADD = mgr.constant(0.0);
                result = leaves[0];
                break;
            case CFLOBDD_FORK:
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
    unsigned int numAExits = node->AConnection.entryPointHandle.handleContents->numExits;
    std::vector<ADD> A_leaves(numAExits);
    for (unsigned int j = 0; j < numAExits; j++) {
        A_leaves[j] = B_ADDs[node->AConnection.returnMapHandle.Lookup(j)];
    }

    ADD result = CFLOBDDConvertNodeToADD(
        mgr, node->AConnection.entryPointHandle, level - 1, topVar, A_leaves, protoMemo);

    // Build and cache proto-ADD for this node using exit-index leaves.
    // All sub-nodes are now in protoMemo (populated by the recursive calls above).
    std::vector<ADD> exitLeaves(numExits);
    for (unsigned int i = 0; i < numExits; i++) {
        exitLeaves[i] = mgr.constant((double)i);
    }

    // Build proto B-ADDs (using cached proto-ADDs of B sub-nodes)
    std::vector<ADD> B_protos(m);
    for (unsigned int i = 0; i < m; i++) {
        CFLOBDDNode* bNode = node->BConnection[i].entryPointHandle.handleContents;
        unsigned int numBExits = bNode->numExits;
        std::vector<ADD> bExitLeaves(numBExits);
        for (unsigned int j = 0; j < numBExits; j++) {
            bExitLeaves[j] = exitLeaves[node->BConnection[i].returnMapHandle.Lookup(j)];
        }
        assert(protoMemo.find(bNode) != protoMemo.end());
        int bOffset = (int)midVar - (int)protoMemo[bNode].topVar;
        std::unordered_map<DdNode*, ADD> spliceMemo;
        B_protos[i] = FusedVariableShiftAndSplice(
            mgr, protoMemo[bNode].protoADD.getNode(), bOffset, bExitLeaves, spliceMemo);
    }

    // Build proto A-ADD (using cached proto-ADD of A sub-node)
    CFLOBDDNode* aNode = node->AConnection.entryPointHandle.handleContents;
    std::vector<ADD> A_protoLeaves(numAExits);
    for (unsigned int j = 0; j < numAExits; j++) {
        A_protoLeaves[j] = B_protos[node->AConnection.returnMapHandle.Lookup(j)];
    }
    assert(protoMemo.find(aNode) != protoMemo.end());
    int aOffset = (int)topVar - (int)protoMemo[aNode].topVar;
    std::unordered_map<DdNode*, ADD> aSpliceMemo;
    ADD protoADD = FusedVariableShiftAndSplice(
        mgr, protoMemo[aNode].protoADD.getNode(), aOffset, A_protoLeaves, aSpliceMemo);

    protoMemo[cfNode] = { protoADD, topVar };
    return result;
}
