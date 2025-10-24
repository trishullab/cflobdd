#ifndef CFLOBDD_C_H
#define CFLOBDD_C_H
#ifdef __cplusplus
extern "C" {
#endif
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>
#include <stdbool.h>
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
// also use ref to indicate that the function will not destroy the input CFLOBDD
typedef size_t CFLOBDD_Ref;

// C-style functions
CFLOBDD_C_API void CFLOBDD_module_init();
CFLOBDD_C_API void CFLOBDD_module_dispose();
// Another "level" is needed since it do not have a "engine" concept.
CFLOBDD_C_API CFLOBDD CFLOBDD_createVar(uint32_t i, int level);
CFLOBDD_C_API CFLOBDD CFLOBDD_createTrue(int level);
CFLOBDD_C_API CFLOBDD CFLOBDD_createFalse(int level);
// destroy is not specifically needed since the destruction function is empty]
// copy constructor
CFLOBDD_C_API CFLOBDD CFLOBDD_copy(const CFLOBDD_Ref cflobdd);
// return the root address of CFLOBDD; can be used in hash
CFLOBDD_C_API void CFLOBDD_delete(CFLOBDD cflobdd);
CFLOBDD_C_API CFLOBDD CFLOBDD_and(const CFLOBDD_Ref a, const CFLOBDD_Ref b);
CFLOBDD_C_API CFLOBDD CFLOBDD_and_to(CFLOBDD product, const CFLOBDD_Ref b);
CFLOBDD_C_API CFLOBDD CFLOBDD_or(const CFLOBDD_Ref a, const CFLOBDD_Ref b);
CFLOBDD_C_API CFLOBDD CFLOBDD_or_to(CFLOBDD sum, const CFLOBDD_Ref b);
CFLOBDD_C_API CFLOBDD CFLOBDD_minus(const CFLOBDD_Ref a, const CFLOBDD_Ref b);
CFLOBDD_C_API CFLOBDD CFLOBDD_not(const CFLOBDD_Ref a);
CFLOBDD_C_API CFLOBDD CFLOBDD_xor(const CFLOBDD_Ref a, const CFLOBDD_Ref b);
CFLOBDD_C_API CFLOBDD CFLOBDD_xor_to(CFLOBDD sum, const CFLOBDD_Ref b);
CFLOBDD_C_API CFLOBDD CFLOBDD_implies(const CFLOBDD_Ref a, const CFLOBDD_Ref b);
CFLOBDD_C_API CFLOBDD CFLOBDD_exists(const CFLOBDD_Ref a, uint32_t i);
CFLOBDD_C_API CFLOBDD CFLOBDD_forall(const CFLOBDD_Ref a, uint32_t i);
CFLOBDD_C_API bool CFLOBDD_isValid(const CFLOBDD_Ref a);
CFLOBDD_C_API uint32_t CFLOBDD_hash(const CFLOBDD_Ref cflobdd);
CFLOBDD_C_API uint32_t CFLOBDD_countNodes(const CFLOBDD_Ref cflobdd);
CFLOBDD_C_API uint32_t CFLOBDD_getLevel(const CFLOBDD_Ref cflobdd);
CFLOBDD_C_API void *CFLOBDD_getRoot(const CFLOBDD_Ref cflobdd);
#ifdef PATH_COUNTING_ENABLED
CFLOBDD_C_API uint32_t CFLOBDD_numSatisfyingAssignments(const CFLOBDD_Ref cflobdd);
#endif

/// Get one satisfying assignment of the CFLOBDD.
/// @param cflobdd The CFLOBDD.
/// @param assignment_buffer The buffer to store the satisfying assignment.
/// @param assignment_buffer_size The size (both in byte and in number of booleans) of the assignment_buffer.
/// @return 0 if no satisfying assignment is found. If there is a satisfiable assignment,
/// The size of buffer used to store the satisfying assignment if success, -1 if the assignment_buffer is too small.
CFLOBDD_C_API ssize_t CFLOBDD_getOneSatisfyingAssignment(const CFLOBDD_Ref cflobdd, bool *assignment_buffer, size_t assignment_buffer_size);
CFLOBDD_C_API bool CFLOBDD_eq(const CFLOBDD_Ref a, const CFLOBDD_Ref b);
CFLOBDD_C_API void CFLOBDD_print(const CFLOBDD_Ref a);
#ifdef __cplusplus
}
#endif

#endif // CFLOBDD_C_H