#include "cuddAbsVal.h"
#include <mpfr.h>

void getAbsValue(CUDD_VALUE_TYPE v, CUDD_VALUE_TYPE* v_abs)
{
  mpfr_set_zero(v_abs->imag, MPFR_RNDN);
  mpfr_mul(v_abs->real, v.real, v.real, MPFR_RNDN);
  mpfr_t tmp; mpfr_init(tmp);
  mpfr_mul(tmp, v.imag, v.imag, MPFR_RNDN);
  mpfr_add(v_abs->real, v_abs->real, tmp, MPFR_RNDN);
  mpfr_clear(tmp);
}