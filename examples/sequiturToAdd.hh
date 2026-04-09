/**
  @file sequiturToAdd.hh

  @brief Convert a SEQUITUR grammar (Sequitur<int>) to a CUDD ADD.

  Uses the non-dense selector-variable approach:
  - Each rule P -> Q1 Q2 ... Qk allocates ceil(log2(k)) selector variables
  - The selector encodes which child Qi; indices >= k map to constant -1
  - Each nonterminal's ADD is built once and memoized
  - Repeated children share the same sub-ADD (canonicalized by CUDD)

  The resulting ADD maps Boolean variable assignments to int terminal values,
  representing the string yielded by the grammar as a function of position.
*/

#ifndef SEQUITUR_TO_ADD_HH_
#define SEQUITUR_TO_ADD_HH_

#include "../cudd-3.0.0/cplusplus/cuddObj.hh"
#include "../sequitur/sequitur.hpp"

#include <unordered_map>
#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>

using namespace jw;

// Result of converting a grammar node to an ADD
struct GrammarADDResult {
    ADD add;           // The ADD for this node's yield
    int maxVar;        // Maximum variable index used (-1 for constants)
    unsigned long yieldLength;  // Length of this node's yield
};

/**
  Build a selector ADD over variables [baseVar, baseVar+numSelectorVars).

  For a rule with k children, numSelectorVars = ceil(log2(k)).
  The selector maps each binary assignment of the selector variables
  to the corresponding child's ADD.  Assignments >= k map to constant -1.

  @param mgr         CUDD manager
  @param childADDs   ADDs for each child (size k)
  @param baseVar     First selector variable index
  @param numSelectorVars  Number of selector variables (ceil(log2(k)))
*/
static ADD buildSelector(Cudd &mgr, const std::vector<ADD> &childADDs,
                         int baseVar, int numSelectorVars)
{
    int k = childADDs.size();
    int numSlots = 1 << numSelectorVars;
    ADD dontCare = mgr.constant(-1);

    // Build bottom-up: start with the leaf-level assignments,
    // then combine with ITE at each selector variable level.

    // Layer 0: the 2^numSelectorVars slots, each holding a child ADD or -1
    std::vector<ADD> current(numSlots);
    for (int i = 0; i < numSlots; i++) {
        current[i] = (i < k) ? childADDs[i] : dontCare;
    }

    // Build the selector tree bottom-up over selector variables.
    // Variable baseVar is the most significant selector bit,
    // baseVar + numSelectorVars - 1 is the least significant.
    for (int level = numSelectorVars - 1; level >= 0; level--) {
        int varIdx = baseVar + level;
        ADD var = mgr.addVar(varIdx);
        int stride = 1 << (numSelectorVars - 1 - level);
        std::vector<ADD> next(current.size() / 2);
        for (int i = 0; i < (int)next.size(); i++) {
            // var=0 selects current[2*i], var=1 selects current[2*i+1]
            next[i] = var.Ite(current[2 * i + 1], current[2 * i]);
        }
        current = std::move(next);
    }

    assert(current.size() == 1);
    return current[0];
}

/**
  Convert a SEQUITUR grammar to a CUDD ADD.

  Walks the grammar DAG bottom-up (memoized per rule ID), building an ADD
  for each nonterminal using the non-dense selector-variable approach.

  @param mgr  CUDD manager
  @param seq  The SEQUITUR grammar (Sequitur<int>)
  @return     GrammarADDResult for the start rule (rule 0)
*/
static GrammarADDResult sequiturToADD(Cudd &mgr, const Sequitur<int> &seq)
{
    const auto &rules = seq.getRules();

    // Type info for dispatching
    const std::type_info &RuleHeadType = typeid(RuleHead);
    const std::type_info &RuleTailType = typeid(RuleTail);
    const std::type_info &RuleSymbolType = typeid(RuleSymbol);
    const std::type_info &ValueType = typeid(ValueSymbol<int>);

    // Parse each rule's children
    struct Child {
        bool isTerminal;
        int terminalValue;       // if isTerminal
        unsigned int ruleID;     // if !isTerminal
    };
    std::unordered_map<unsigned int, std::vector<Child>> ruleChildren;

    for (const auto &pair : rules) {
        unsigned int ruleID = pair.first;
        Symbol *head = pair.second;
        std::vector<Child> children;
        Symbol *sym = head->next();
        while (typeid(*sym) != RuleTailType) {
            if (typeid(*sym) == ValueType) {
                children.push_back({true, static_cast<const ValueSymbol<int>*>(sym)->getValue(), 0});
            } else if (typeid(*sym) == RuleSymbolType) {
                children.push_back({false, 0, static_cast<const RuleSymbol*>(sym)->getID()});
            }
            sym = sym->next();
        }
        ruleChildren[ruleID] = std::move(children);
    }

    // Phase 1: Top-down variable allocation.
    // Each rule gets a baseVar; its selector variables are [baseVar, baseVar + numSelectorVars).
    // Children's variables start at baseVar + numSelectorVars.
    // Since all nonterminal children share the same variable space (they are memoized),
    // a child nonterminal's baseVar is baseVar + numSelectorVars.
    std::unordered_map<unsigned int, int> ruleBaseVar;   // first variable for this rule's selector
    std::unordered_map<unsigned int, int> ruleTotalVars;  // total variables used by this rule's subtree

    std::function<int(unsigned int, int)> assignVars;
    assignVars = [&](unsigned int ruleID, int base) -> int {
        // Returns total number of variables used by this rule's subtree
        auto it = ruleTotalVars.find(ruleID);
        if (it != ruleTotalVars.end()) {
            ruleBaseVar[ruleID] = base;  // update base for this context
            return it->second;
        }

        const auto &children = ruleChildren[ruleID];
        int k = children.size();

        if (k == 1 && children[0].isTerminal) {
            ruleBaseVar[ruleID] = base;
            ruleTotalVars[ruleID] = 0;
            return 0;
        }

        int numSelectorVars = (k <= 1) ? 0 : static_cast<int>(std::ceil(std::log2(k)));
        ruleBaseVar[ruleID] = base;

        // Find max child depth (all nonterminal children share same variable range)
        int maxChildVars = 0;
        for (int i = 0; i < k; i++) {
            if (!children[i].isTerminal) {
                int cv = assignVars(children[i].ruleID, base + numSelectorVars);
                if (cv > maxChildVars) maxChildVars = cv;
            }
        }

        int total = numSelectorVars + maxChildVars;
        ruleTotalVars[ruleID] = total;
        return total;
    };
    assignVars(0, 0);

    // Phase 2: Bottom-up ADD construction (memoized per rule ID).
    std::unordered_map<unsigned int, GrammarADDResult> memo;

    std::function<GrammarADDResult(unsigned int)> convertRule;
    convertRule = [&](unsigned int ruleID) -> GrammarADDResult {
        auto it = memo.find(ruleID);
        if (it != memo.end()) return it->second;

        const auto &children = ruleChildren[ruleID];
        int k = children.size();
        int base = ruleBaseVar[ruleID];
        int totalVars = ruleTotalVars[ruleID];

        // Base case: single terminal
        if (k == 1 && children[0].isTerminal) {
            GrammarADDResult result;
            result.add = mgr.constant(static_cast<double>(children[0].terminalValue));
            result.maxVar = base - 1;
            result.yieldLength = 1;
            memo[ruleID] = result;
            return result;
        }

        int numSelectorVars = (k <= 1) ? 0 : static_cast<int>(std::ceil(std::log2(k)));
        int childBase = base + numSelectorVars;
        int maxVar = base + totalVars - 1;

        // Build child ADDs
        // For terminals: narrow to all-zeros of vars [childBase..maxVar] → value, else → -1
        std::vector<ADD> childADDs(k);
        unsigned long totalYield = 0;

        for (int i = 0; i < k; i++) {
            if (children[i].isTerminal) {
                ADD termADD = mgr.constant(static_cast<double>(children[i].terminalValue));
                if (maxVar >= childBase) {
                    ADD dontCare = mgr.constant(-1);
                    for (int v = childBase; v <= maxVar; v++) {
                        termADD = mgr.addVar(v).Ite(dontCare, termADD);
                    }
                }
                childADDs[i] = termADD;
                totalYield += 1;
            } else {
                GrammarADDResult childResult = convertRule(children[i].ruleID);
                ADD childADD = childResult.add;
                // If this child uses fewer variables than maxVar, narrow the
                // extra lower variables to all-zeros → pass through, else → -1
                if (childResult.maxVar < maxVar) {
                    ADD dontCare = mgr.constant(-1);
                    for (int v = childResult.maxVar + 1; v <= maxVar; v++) {
                        childADD = mgr.addVar(v).Ite(dontCare, childADD);
                    }
                }
                childADDs[i] = childADD;
                totalYield += childResult.yieldLength;
            }
        }

        // Build selector
        ADD resultADD;
        if (numSelectorVars == 0) {
            resultADD = childADDs[0];
        } else {
            resultADD = buildSelector(mgr, childADDs, base, numSelectorVars);
        }

        GrammarADDResult result;
        result.add = resultADD;
        result.maxVar = maxVar;
        result.yieldLength = totalYield;
        memo[ruleID] = result;
        return result;
    };

    return convertRule(0);
}

/**
  Print statistics about the grammar and resulting ADD.
*/
static void printSequiturADDStats(const Sequitur<int> &seq,
                                  const GrammarADDResult &result)
{
    const auto &rules = seq.getRules();

    // Count grammar symbols
    const std::type_info &RuleHeadType = typeid(RuleHead);
    const std::type_info &RuleTailType = typeid(RuleTail);
    unsigned int totalSymbols = 0;
    unsigned int ruleCount = 0;
    for (const auto &pair : rules) {
        ruleCount++;
        Symbol *sym = pair.second->next();
        while (typeid(*sym) != RuleTailType) {
            totalSymbols++;
            sym = sym->next();
        }
    }

    std::cout << "=== SEQUITUR Grammar ===" << std::endl;
    std::cout << "  Rules: " << ruleCount << std::endl;
    std::cout << "  Total RHS symbols: " << totalSymbols << std::endl;
    std::cout << "  Yield length: " << result.yieldLength << std::endl;
    std::cout << "=== ADD ===" << std::endl;
    std::cout << "  Variables used: " << (result.maxVar + 1) << std::endl;
    std::cout << "  ADD nodes: " << result.add.nodeCount() << std::endl;
}

#endif // SEQUITUR_TO_ADD_HH_
