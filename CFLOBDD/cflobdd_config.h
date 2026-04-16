#ifndef CFLOBDD_CONFIG_H
#define CFLOBDD_CONFIG_H

#include <cstddef>

// Tunable runtime parameters for the CFLOBDD library.
// All fields have sensible defaults; override via command-line flags
// (see --help) or by assigning to cflobddConfig before calling any
// CFLOBDD operations.
struct CFLOBDDConfig {
    // --- Flat-array vs hash-map thresholds ---
    // When the pair/reduction space is <= threshold, use an O(1) flat array;
    // otherwise fall back to boost::unordered_flat_map.
    size_t flatLookupThreshold    = 33554432;  // 2^25 (PairProduct)
    size_t composeFlatThreshold   = 33554432;  // 2^25 (ComposeAndReduce)
    size_t identityMapArrayThreshold = 1024;   // MakeIdentityReturnMap

    // --- Object-pool freelist caps ---
    // 0 = unlimited (no eviction).
    size_t nodeFreelistCap           = 64;     // CFLOBDDInternalNode
    size_t pairProductFreelistCap    = 64;     // PairProductMapBody
    size_t tripleProductFreelistCap  = 64;     // TripleProductMapBody
    size_t reductionMapFreelistCap   = 64;     // ReductionMapBody
    size_t returnMapFreelistCap      = 0;      // ReturnMapBody (unlimited by default)
};

extern CFLOBDDConfig cflobddConfig;

void printCFLOBDDUsage(const char* progName);

#endif // CFLOBDD_CONFIG_H
