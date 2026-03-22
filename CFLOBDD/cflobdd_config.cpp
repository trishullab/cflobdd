#include <iostream>
#include <cstring>
#include "cflobdd_config.h"

CFLOBDDConfig cflobddConfig;

void printCFLOBDDUsage(const char* progName)
{
    std::cout
        << "Usage: " << progName << " [options] <test> [args...]\n"
        << "\n"
        << "Options:\n"
        << "  -h, --help                          Show this help message and exit\n"
        << "  -v, --verbose                       Enable verbose output\n"
        << "\n"
        << "Flat-array thresholds (0 = always use hash map):\n"
        << "  --flat-lookup-threshold=N            PairProduct flat-array cutoff\n"
        << "                                       (default: " << CFLOBDDConfig{}.flatLookupThreshold << ")\n"
        << "  --compose-flat-threshold=N           ComposeAndReduce flat-array cutoff\n"
        << "                                       (default: " << CFLOBDDConfig{}.composeFlatThreshold << ")\n"
        << "  --identity-map-threshold=N           MakeIdentityReturnMap array cutoff\n"
        << "                                       (default: " << CFLOBDDConfig{}.identityMapArrayThreshold << ")\n"
        << "\n"
        << "Freelist caps (0 = unlimited):\n"
        << "  --freelist-cap=N                     Set all freelist caps at once\n"
        << "  --node-freelist-cap=N                CFLOBDDInternalNode pool cap\n"
        << "                                       (default: " << CFLOBDDConfig{}.nodeFreelistCap << ")\n"
        << "  --pair-product-freelist-cap=N        PairProductMapBody pool cap\n"
        << "                                       (default: " << CFLOBDDConfig{}.pairProductFreelistCap << ")\n"
        << "  --triple-product-freelist-cap=N      TripleProductMapBody pool cap\n"
        << "                                       (default: " << CFLOBDDConfig{}.tripleProductFreelistCap << ")\n"
        << "  --reduction-map-freelist-cap=N       ReductionMapBody pool cap\n"
        << "                                       (default: " << CFLOBDDConfig{}.reductionMapFreelistCap << ")\n"
        << "  --return-map-freelist-cap=N          ReturnMapBody pool cap\n"
        << "                                       (default: " << CFLOBDDConfig{}.returnMapFreelistCap << ")\n"
        << std::endl;
}
