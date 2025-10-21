#include "cflobdd_c.h"
#include "cflobdd_int.h"

static_assert(sizeof(CFL_OBDD::CFLOBDD) == sizeof(CFLOBDD));
extern "C" {


CFLOBDD CFLOBDD_createVar(uint32_t i, int level) {
    auto ans = CFL_OBDD::MkProjection(i, level);
    return *(CFLOBDD *)(&ans);
}

CFLOBDD CFLOBDD_copy(CFLOBDD cflobdd) {
    auto a1 = *(CFL_OBDD::CFLOBDD *)(&cflobdd);
    auto ans = CFL_OBDD::CFLOBDD(a1);
    return *(CFLOBDD *)(&ans);
}



CFLOBDD CFLOBDD_and(CFLOBDD a, CFLOBDD b) {
    auto a1 = *(CFL_OBDD::CFLOBDD *)(&a);
    auto b1 = *(CFL_OBDD::CFLOBDD *)(&b);
    auto ans = CFL_OBDD::MkAnd(a1, b1);
    return *(CFLOBDD *)(&ans);
}

CFLOBDD CFLOBDD_or(CFLOBDD a, CFLOBDD b) {
    auto a1 = *(CFL_OBDD::CFLOBDD *)(&a);
    auto b1 = *(CFL_OBDD::CFLOBDD *)(&b);
    auto ans = CFL_OBDD::MkOr(a1, b1);
    return *(CFLOBDD *)(&ans);
}

CFLOBDD CFLOBDD_not(CFLOBDD a) {
    auto a1 = *(CFL_OBDD::CFLOBDD *)(&a);
    auto ans = CFL_OBDD::MkNot(a1);
    return *(CFLOBDD *)(&ans);
}

CFLOBDD CFLOBDD_xor(CFLOBDD a, CFLOBDD b) {
    auto a1 = *(CFL_OBDD::CFLOBDD *)(&a);
    auto b1 = *(CFL_OBDD::CFLOBDD *)(&b);
    auto ans = CFL_OBDD::MkExclusiveOr(a1, b1);
    return *(CFLOBDD *)(&ans);
}

CFLOBDD CFLOBDD_exists(CFLOBDD a, uint32_t i) {
    auto a1 = *(CFL_OBDD::CFLOBDD *)(&a);
    auto ans = CFL_OBDD::MkExists(a1, i);
    return *(CFLOBDD *)(&ans);
}

CFLOBDD CFLOBDD_forall(CFLOBDD a, uint32_t i) {
    auto a1 = *(CFL_OBDD::CFLOBDD *)(&a);
    auto ans = CFL_OBDD::MkForall(a1, i);
    return *(CFLOBDD *)(&ans);
}
}