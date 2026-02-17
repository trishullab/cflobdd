#include "cflobdd_c.h"
#include "cflobdd_t.h"
#include "cflobdd_int.h"
#include <utility>
#include "vector_float_boost.h"
// #include "vector_complex_float_boost.h"
#include "quantum_algos.h"
// #include "matrix1234_fourier.h"
// #include "Solver/uwr/matrix/HowellMatrix.h"
// #include "Solver/uwr/matrix/ModularSquareMatrix.h"
// #include "memory_check.h"
#include "matmult_map.h"
// #include "matrix1234_node.h"
#include "matrix1234_int.h"
#ifdef WCFLOBDD_SUPPORTED
#include "wmatrix1234_fb_mul.h"
#include "weighted_cross_product.h"
#include "wvector_fb_mul.h"
#include "weighted_quantum_algos.h"
#include "wvector_complex_fb_mul.h"
#include "wmatrix1234_fourier_mul.h"
#include "wvector_fourier_mul.h"
#include "weighted_bdd_node_t.h"
#include "weighted_cross_product_bdd.h"
#include "wvector_complex_fb_mul_bdd_node.h"
#endif
static_assert(sizeof(CFL_OBDD::CFLOBDD) == sizeof(CFLOBDD));
extern "C" {

void CFLOBDD_module_init() {
    
	CFL_OBDD::CFLOBDDNodeHandle::InitNoDistinctionTable();
	// CFLOBDDNodeHandle::InitAdditionInterleavedTable();
	CFL_OBDD::CFLOBDDNodeHandle::InitReduceCache();
	CFL_OBDD::InitPairProductCache();
	CFL_OBDD::InitTripleProductCache();
	CFL_OBDD::Matrix1234Int::Matrix1234Initializer();
	CFL_OBDD::VectorFloatBoost::VectorInitializer();

#ifdef WCFLOBDD_SUPPORTED
	// typedef double BIG_FLOAT;
	CFL_OBDD::WeightedCFLOBDDNodeHandleT<CFL_OBDD::BIG_FLOAT, std::multiplies<CFL_OBDD::BIG_FLOAT>>::InitNoDistinctionTable();
	CFL_OBDD::WeightedCFLOBDDNodeHandleT<CFL_OBDD::BIG_FLOAT, std::multiplies<CFL_OBDD::BIG_FLOAT>>::InitNoDistinctionTable_Ann();
	CFL_OBDD::WeightedCFLOBDDNodeHandleT<CFL_OBDD::BIG_FLOAT, std::multiplies<CFL_OBDD::BIG_FLOAT>>::InitIdentityNodeTable();
	CFL_OBDD::WeightedCFLOBDDNodeHandleT<CFL_OBDD::BIG_FLOAT, std::multiplies<CFL_OBDD::BIG_FLOAT>>::InitReduceCache();
	CFL_OBDD::WeightedMatrix1234FloatBoostMul::Matrix1234Initializer();
	CFL_OBDD::WeightedVectorFloatBoostMul::VectorInitializer();
	CFL_OBDD::InitWeightedPairProductCache<CFL_OBDD::BIG_FLOAT, std::multiplies<CFL_OBDD::BIG_FLOAT>>();

	CFL_OBDD::WeightedCFLOBDDNodeHandleT<CFL_OBDD::BIG_COMPLEX_FLOAT, std::multiplies<CFL_OBDD::BIG_COMPLEX_FLOAT>>::InitNoDistinctionTable();
	CFL_OBDD::WeightedCFLOBDDNodeHandleT<CFL_OBDD::BIG_COMPLEX_FLOAT, std::multiplies<CFL_OBDD::BIG_COMPLEX_FLOAT>>::InitNoDistinctionTable_Ann();
	CFL_OBDD::WeightedCFLOBDDNodeHandleT<CFL_OBDD::BIG_COMPLEX_FLOAT, std::multiplies<CFL_OBDD::BIG_COMPLEX_FLOAT>>::InitIdentityNodeTable();
	CFL_OBDD::WeightedCFLOBDDNodeHandleT<CFL_OBDD::BIG_COMPLEX_FLOAT, std::multiplies<CFL_OBDD::BIG_COMPLEX_FLOAT>>::InitReduceCache();
	CFL_OBDD::WeightedMatrix1234ComplexFloatBoostMul::Matrix1234Initializer();
	CFL_OBDD::WeightedVectorComplexFloatBoostMul::VectorInitializer();
	CFL_OBDD::InitWeightedPairProductCache<CFL_OBDD::BIG_COMPLEX_FLOAT, std::multiplies<CFL_OBDD::BIG_COMPLEX_FLOAT>>();

	CFL_OBDD::WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>>::InitNoDistinctionTable();
	CFL_OBDD::WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>>::InitNoDistinctionTable_Ann();
	CFL_OBDD::WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>>::InitIdentityNodeTable();
	CFL_OBDD::WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>>::InitReduceCache();
	CFL_OBDD::WeightedMatrix1234FourierMul::Matrix1234Initializer();
	CFL_OBDD::WeightedVectorFourierMul::VectorInitializer();
	CFL_OBDD::InitWeightedPairProductCache<fourierSemiring, std::multiplies<fourierSemiring>>();

	CFL_OBDD::WeightedBDDNodeHandle<CFL_OBDD::BIG_COMPLEX_FLOAT, std::multiplies<CFL_OBDD::BIG_COMPLEX_FLOAT>>::InitLeafNodes();
	CFL_OBDD::InitWeightedBDDPairProductCache<CFL_OBDD::BIG_COMPLEX_FLOAT, std::multiplies<CFL_OBDD::BIG_COMPLEX_FLOAT>>();
#endif
}

void CFLOBDD_module_dispose() {
	CFL_OBDD::DisposeOfTripleProductCache();
	CFL_OBDD::DisposeOfPairProductCache();
	CFL_OBDD::CFLOBDDNodeHandle::DisposeOfReduceCache();
#ifdef WCFLOBDD_SUPPORTED
	CFL_OBDD::DisposeOfWeightedPairProductCache<CFL_OBDD::BIG_FLOAT, std::multiplies<CFL_OBDD::BIG_FLOAT>>();
	CFL_OBDD::DisposeOfWeightedPairProductCache<CFL_OBDD::BIG_COMPLEX_FLOAT, std::multiplies<CFL_OBDD::BIG_COMPLEX_FLOAT>>();
	CFL_OBDD::DisposeOfWeightedPairProductCache<fourierSemiring, std::multiplies<fourierSemiring>>();

	CFL_OBDD::DisposeOfWeightedBDDPairProductCache<CFL_OBDD::BIG_COMPLEX_FLOAT, std::multiplies<CFL_OBDD::BIG_COMPLEX_FLOAT>>();
#endif
}

class CFLOBDD_Nodestruct {
private:
    union {
        CFL_OBDD::CFLOBDD cflobdd;
    };
public:
    CFLOBDD_Nodestruct(CFL_OBDD::CFLOBDD &&cflobdd) : cflobdd(std::move(cflobdd)) {}
    ~CFLOBDD_Nodestruct() {}
    CFLOBDD get_value() {
        return *(CFLOBDD *)&cflobdd;
    }
};

CFLOBDD cast(CFL_OBDD::CFLOBDD &&cflobdd) {
    CFLOBDD_Nodestruct a1(std::move(cflobdd));
    return a1.get_value();
}
CFL_OBDD::CFLOBDD retrieve(CFLOBDD cflobdd) {
    return CFL_OBDD::CFLOBDD(std::move(*(CFL_OBDD::CFLOBDD *)(&cflobdd)));
}

void forget(CFL_OBDD::CFLOBDD &&cflobdd) {
    CFLOBDD_Nodestruct a1(std::move(cflobdd));
}

CFLOBDD CFLOBDD_createVar(uint32_t i, int level) {
    return cast(std::move(CFL_OBDD::MkProjection(i, level)));
}
CFLOBDD CFLOBDD_createTrue(int level) {
    return cast(std::move(CFL_OBDD::MkTrue(level)));
}
CFLOBDD CFLOBDD_createFalse(int level) {
    return cast(std::move(CFL_OBDD::MkFalse(level)));
}

void CFLOBDD_delete(CFLOBDD cflobdd) {
    auto a1 = retrieve(cflobdd);
    // a1 will go out of scope and call destructor
}

// return the root address of CFLOBDD; can be used in hash
CFLOBDD_C_API void *CFLOBDD_getRoot(const CFLOBDD_Ref cflobdd) {
    auto a1 = retrieve(cflobdd);
    void *ans = a1.root.get_ptr();
    forget(std::move(a1));
    return ans;
}

CFLOBDD CFLOBDD_copy(const CFLOBDD_Ref cflobdd) {
    auto a1 = retrieve(cflobdd);
    auto ans = cast(CFL_OBDD::CFLOBDD(a1));
    forget(std::move(a1));
    return ans;
}



CFLOBDD CFLOBDD_and(const CFLOBDD_Ref a, const CFLOBDD_Ref b) {
    auto a1 = retrieve(a);
    auto b1 = retrieve(b);
    auto ans = cast(CFL_OBDD::MkAnd(a1, b1));
    forget(std::move(a1));
    forget(std::move(b1));
    return ans;
}

CFLOBDD CFLOBDD_and_to(CFLOBDD product, const CFLOBDD_Ref b) {
    auto a1 = retrieve(product);
    auto b1 = retrieve(b);
    auto ans = cast(CFL_OBDD::MkAnd(a1, b1));
    // a1 calls destructor when ans goes out of scope
    // forget(a1);
    forget(std::move(b1));
    return ans;
}

CFLOBDD CFLOBDD_or(const CFLOBDD_Ref a, const CFLOBDD_Ref b) {
    auto a1 = retrieve(a);
    auto b1 = retrieve(b);
    auto ans = cast(CFL_OBDD::MkOr(a1, b1));
    forget(std::move(a1));
    forget(std::move(b1));
    return ans;
}

CFLOBDD CFLOBDD_or_to(CFLOBDD sum, const CFLOBDD_Ref b) {
    auto a1 = retrieve(sum);
    auto b1 = retrieve(b);
    auto ans = cast(CFL_OBDD::MkOr(a1, b1));
    // a1 calls destructor when ans goes out of scope
    // forget(a1);
    forget(std::move(b1));
    return ans;
}

CFLOBDD CFLOBDD_minus(const CFLOBDD_Ref a, const CFLOBDD_Ref b) {
    auto a1 = retrieve(a);
    auto b1 = retrieve(b);
    auto ans = cast(CFL_OBDD::MkMinus(a1, b1));
    forget(std::move(a1));
    forget(std::move(b1));
    return ans;
}

CFLOBDD CFLOBDD_not(const CFLOBDD_Ref a) {
    auto a1 = retrieve(a);
    auto ans = cast(CFL_OBDD::MkNot(a1));
    forget(std::move(a1));
    return ans;
}

CFLOBDD CFLOBDD_xor(const CFLOBDD_Ref a, const CFLOBDD_Ref b) {
    auto a1 = retrieve(a);
    auto b1 = retrieve(b);
    auto ans = cast(CFL_OBDD::MkExclusiveOr(a1, b1));
    forget(std::move(a1));
    forget(std::move(b1));
    return ans;
}

CFLOBDD CFLOBDD_xor_to(CFLOBDD sum, const CFLOBDD_Ref b) {
    auto a1 = retrieve(sum);
    auto b1 = retrieve(b);
    auto ans = cast(CFL_OBDD::MkExclusiveOr(a1, b1));
    // a1 calls destructor when ans goes out of scope
    // forget(a1);
    forget(std::move(b1));
    return ans;
}

CFLOBDD CFLOBDD_implies(const CFLOBDD_Ref a, const CFLOBDD_Ref b) {
    auto a1 = retrieve(a);
    auto b1 = retrieve(b);
    auto ans = cast(CFL_OBDD::MkImplies(a1, b1));
    forget(std::move(a1));
    forget(std::move(b1));
    return ans;
}

CFLOBDD CFLOBDD_exists(const CFLOBDD_Ref a, uint32_t i) {
    auto a1 = retrieve(a);
    auto ans = cast(CFL_OBDD::MkExists(a1, i));
    forget(std::move(a1));
    return ans;
}

CFLOBDD CFLOBDD_forall(const CFLOBDD_Ref a, uint32_t i) {
    auto a1 = retrieve(a);
    auto ans = cast(CFL_OBDD::MkForall(a1, i));
    forget(std::move(a1));
    return ans;
}
bool CFLOBDD_isValid(const CFLOBDD_Ref a) {
    auto a1 = retrieve(a);
    auto ans = a1.IsValid();
    forget(std::move(a1));
    return ans;
}
CFLOBDD_C_API uint32_t CFLOBDD_getLevel(const CFLOBDD_Ref cflobdd) {
    auto a1 = retrieve(cflobdd);
    auto ans = (*a1.root).level;
    forget(std::move(a1));
    return ans;
}
#ifdef PATH_COUNTING_ENABLED
CFLOBDD_C_API uint32_t CFLOBDD_numSatisfyingAssignments(const CFLOBDD_Ref cflobdd) {
    auto a1 = retrieve(cflobdd);
    auto ans = a1.NumSatisfyingAssignments();
    forget(std::move(a1));
    return ans;
}
#endif

CFLOBDD_C_API ssize_t CFLOBDD_getOneSatisfyingAssignment(const CFLOBDD_Ref cflobdd, bool *assignment_buffer, size_t assignment_buffer_size) {
    auto a1 = retrieve(cflobdd);
    uint32_t level_size = 1 << ((*a1.root).level);
    uint32_t level_start = (1 << CFL_OBDD::CFLOBDD::maxLevel) - level_size;
    SH_OBDD::Assignment *assignment;
    bool ans = a1.FindOneSatisfyingAssignment(assignment);
    forget(std::move(a1));
    if (!ans) {
        return 0;
    }
    size_t copy_size = assignment_buffer_size < level_size ? assignment_buffer_size : level_size;
    memcpy(assignment_buffer, assignment->get_data() + level_start, copy_size);
    delete assignment;
    return copy_size;
}

bool CFLOBDD_eq(const CFLOBDD_Ref a, const CFLOBDD_Ref b) {
    auto a1 = retrieve(a);
    auto b1 = retrieve(b);
    auto ans = a1 == b1;
    forget(std::move(a1));
    forget(std::move(b1));
    return ans;
}

uint32_t CFLOBDD_hash(const CFLOBDD_Ref cflobdd) {
    auto a1 = retrieve(cflobdd);
    uint32_t ans = (uint32_t)a1.Hash();
    forget(std::move(a1));
    return ans;
}

uint32_t CFLOBDD_countNodes(const CFLOBDD_Ref cflobdd) {
    auto a1 = retrieve(cflobdd);
    uint32_t ans = 0;
    a1.CountNodes(ans);
    forget(std::move(a1));
    return ans;
}

void CFLOBDD_print(const CFLOBDD_Ref a) {
    auto a1 = retrieve(a);
    CFL_OBDD::PrintCFLOBDD(a1);
    forget(std::move(a1));
}
}