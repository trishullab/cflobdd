#ifndef CFLOBDD_C_H
#define CFLOBDD_C_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdint.h>
#include <stdlib.h>
#if defined(_WIN32)
    #ifdef CFLOBDD_C_EXPORTS
        #define CFLOBDD_C_API __declspec(dllexport)
    #else
        #define CFLOBDD_C_API __declspec(dllimport)
    #endif
#else
    #define CFLOBDD_C_API __attribute__((visibility("default")))
#endif

// Opaque pointer to the C++ CFLOBDD object
typedef size_t CFLOBDD;

// C-style functions
// Another "level" is needed since it do not have a "engine" concept.
CFLOBDD_C_API CFLOBDD CFLOBDD_createVar(uint32_t i, int level);
// destroy is not specifically needed since the destruction function is empty]
// copy constructor
CFLOBDD_C_API CFLOBDD CFLOBDD_copy(CFLOBDD cflobdd);

CFLOBDD_C_API CFLOBDD CFLOBDD_and(CFLOBDD a, CFLOBDD b);
CFLOBDD_C_API CFLOBDD CFLOBDD_or(CFLOBDD a, CFLOBDD b);
CFLOBDD_C_API CFLOBDD CFLOBDD_not(CFLOBDD a);
CFLOBDD_C_API CFLOBDD CFLOBDD_xor(CFLOBDD a, CFLOBDD b);
CFLOBDD_C_API CFLOBDD CFLOBDD_exists(CFLOBDD a, uint32_t i);
CFLOBDD_C_API CFLOBDD CFLOBDD_forall(CFLOBDD a, uint32_t i);
#ifdef __cplusplus
}
#endif

#endif // CFLOBDD_C_H