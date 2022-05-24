/**
  @file

  @ingroup cudd

  @brief Apply functions for ADDs and their operators.

  @author Fabio Somenzi

  @copyright@parblock
  Copyright (c) 1995-2015, Regents of the University of Colorado

  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  Neither the name of the University of Colorado nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
  @endparblock

*/

#include "util.h"
#include "cuddInt.h"


/*---------------------------------------------------------------------------*/
/* Constant declarations                                                     */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Stucture declarations                                                     */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Type declarations                                                         */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Variable declarations                                                     */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/

/** \cond */

/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/

/** \endcond */


/*---------------------------------------------------------------------------*/
/* Definition of exported functions                                          */
/*---------------------------------------------------------------------------*/

/**
  @brief Applies op to the corresponding discriminants of f and g.

  @return a pointer to the result if succssful; NULL otherwise.

  @sideeffect None

  @see Cudd_addMonadicApply Cudd_addPlus Cudd_addTimes
  Cudd_addThreshold Cudd_addSetNZ Cudd_addDivide Cudd_addMinus Cudd_addMinimum
  Cudd_addMaximum Cudd_addOneZeroMaximum Cudd_addDiff Cudd_addAgreement
  Cudd_addOr Cudd_addNand Cudd_addNor Cudd_addXor Cudd_addXnor

*/
DdNode *
Cudd_addApply(
  DdManager * dd /**< manager */,
  DD_AOP op /**< operator */,
  DdNode * f /**< first operand */,
  DdNode * g /**< second operand */)
{
    DdNode *res;

    do {
    	dd->reordered = 0;
      // printf("apply %g %g\n", Cudd_V(f), Cudd_V(g));
      // mpfr_printf("%.16Rf\n", Cudd_V(f).real);
        res = cuddAddApplyRecur(dd,op,f,g);
  } while (dd->reordered == 1);
    if (dd->errorCode == CUDD_TIMEOUT_EXPIRED && dd->timeoutHandler) {
        dd->timeoutHandler(dd, dd->tohArg);
    }

    return(res);

} /* end of Cudd_addApply */


/**
  @brief Integer and floating point addition.

  @return NULL if not a terminal case; f+g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addPlus(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *res;
    DdNode *F, *G;
    CUDD_VALUE_TYPE value;
    int is_complex_assigned = cuddV(*f).is_complex_assigned;
    value.is_complex_assigned = is_complex_assigned;
    if (is_complex_assigned == 1){
      mpfr_init(value.real);
      mpfr_init(value.imag);
    }

    F = *f; G = *g;
    if (F == DD_ZERO(dd)){
      if (is_complex_assigned == 1){
        mpfr_clear(value.real); mpfr_clear(value.imag); 
      }
      return(G);
    }
    if (G == DD_ZERO(dd)) {
      if (is_complex_assigned == 1){
        mpfr_clear(value.real); mpfr_clear(value.imag);
      }
      return(F);
    }
    if (cuddIsConstant(F) && cuddIsConstant(G)) {
      if (is_complex_assigned == 1){
        mpfr_add(value.real, cuddV(F).real, cuddV(G).real, RND_TYPE);
        mpfr_add(value.imag, cuddV(F).imag, cuddV(G).imag, RND_TYPE);
        // value = cuddV(F)+cuddV(G);
        res = cuddUniqueConst(dd,value);
        mpfr_clear(value.real);
        mpfr_clear(value.imag);
        return(res);
      }
    }
    if (F > G) { /* swap f and g */
	*f = G;
	*g = F;
    }
    if (is_complex_assigned == 1){
      mpfr_clear(value.real);
      mpfr_clear(value.imag);
    }
    return(NULL);

} /* end of Cudd_addPlus */


/**
  @brief Integer and floating point multiplication.

  @details This function can be used also to take the AND of two 0-1
  ADDs.

  @return NULL if not a terminal case; f * g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addTimes(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *res;
    DdNode *F, *G;
    CUDD_VALUE_TYPE value;
    int is_complex_assigned = cuddV(*f).is_complex_assigned;
    value.is_complex_assigned = is_complex_assigned;
    if (is_complex_assigned == 1){
      mpfr_init(value.real);
      mpfr_init(value.imag);
    }

    F = *f; G = *g;
    // printf("0\n");
    if ((F == DD_ZERO(dd)) || (G == DD_ZERO(dd))){
      if (is_complex_assigned == 1){
        mpfr_clear(value.real);mpfr_clear(value.imag); 
      }
      return(DD_ZERO(dd));
    }
    // printf("1\n");
    if (F == DD_ONE(dd)){
      if (is_complex_assigned == 1){
        mpfr_clear(value.real);mpfr_clear(value.imag); 
      }
        return(G);
    }
    // printf("2\n");
    if (G == DD_ONE(dd)){
      if (is_complex_assigned == 1){
        mpfr_clear(value.real);mpfr_clear(value.imag); 
      }
      return(F);
    }
    // printf("3\n");
    if (cuddIsConstant(F) && cuddIsConstant(G)) {
      if (is_complex_assigned == 1){
        mpfr_mul(value.real, cuddV(F).real, cuddV(G).real, RND_TYPE);
        mpfr_t tmp_val; mpfr_init(tmp_val);
        mpfr_mul(tmp_val, cuddV(F).imag, cuddV(G).imag, RND_TYPE);
        mpfr_sub(value.real, value.real, tmp_val, RND_TYPE); // real

        mpfr_mul(value.imag, cuddV(F).real, cuddV(G).imag, RND_TYPE);
        mpfr_mul(tmp_val, cuddV(F).imag, cuddV(G).real, RND_TYPE);
        mpfr_add(value.imag, value.imag, tmp_val, RND_TYPE); // imag
        mpfr_clear(tmp_val);

        // value = cuddV(F)*cuddV(G);
      	res = cuddUniqueConst(dd,value);
        mpfr_clear(value.real);
        mpfr_clear(value.imag);
      }else{
        if (cuddV(F).base > cuddV(G).base){
          //printf("%d %d %d %d %d\n", cuddV(F).val, cuddV(F).base, cuddV(G).val, cuddV(G).base, (cuddV(F).val * (cuddV(G).base/cuddV(F).base) + cuddV(G).val));
          value.val = (cuddV(F).val * (cuddV(G).base/cuddV(F).base) + cuddV(G).val) % cuddV(G).base;
          value.base = cuddV(G).base;
          res = cuddUniqueConst(dd, value);
        }else if (cuddV(G).base > cuddV(F).base){
          value.val = (cuddV(G).val * (cuddV(F).base/cuddV(G).base) + cuddV(F).val) % cuddV(F).base;
          value.base = cuddV(F).base;
          res = cuddUniqueConst(dd, value);
        } else{
          value.val = (cuddV(F).val + cuddV(G).val) % cuddV(G).base;
          value.base = cuddV(F).base;
          res = cuddUniqueConst(dd, value);
        }
      }
    	return(res);
    }
    if (F > G) { /* swap f and g */
	*f = G;
	*g = F;
    }
    if (is_complex_assigned == 1){
      mpfr_clear(value.real);
      mpfr_clear(value.imag);
    }
    return(NULL);

} /* end of Cudd_addTimes */

/**
  @brief Converts -1 to 0
*/
DdNode *
Cudd_addConvertToBase(
  DdManager * dd,
  DdNode * f,
  unsigned int base)
{

    if (cuddIsConstant(f)) {
      CUDD_VALUE_TYPE value;
      value.val = (cuddV(f).val * (base/cuddV(f).base)) % base;
      value.base = base;
      value.is_complex_assigned = 0;
      DdNode *res = cuddUniqueConst(dd,value);
      return(res);
    }
    return(NULL);

} /* end of Cudd_addSimonsRemoveMinusOne */

    DdNode *
Cudd_addConvertToComplex(
  DdManager * dd,
  DdNode * f)
{

    if (cuddIsConstant(f)) {
      CUDD_VALUE_TYPE value;
      value.val = cuddV(f).val;
      value.base = cuddV(f).base;
      value.is_complex_assigned = 1;

      mpfr_init(value.real); mpfr_init(value.imag);
      mpfr_const_pi(value.real, RND_TYPE); mpfr_const_pi(value.imag, RND_TYPE);
      mpfr_mul_d(value.real, value.real, (2.0 * value.val)/(double)pow(2, value.base), RND_TYPE);
      mpfr_mul_d(value.imag, value.imag, (2.0 * value.val)/(double)pow(2, value.base), RND_TYPE);
      mpfr_cos(value.real, value.real, RND_TYPE);
      mpfr_sin(value.imag, value.imag, RND_TYPE);

      DdNode *res = cuddUniqueConst(dd,value);
      mpfr_clear(value.real); mpfr_clear(value.imag);
      return(res);
    }
    return(NULL);

} /* end of Cudd_addSimonsRemoveMinusOne */

    DdNode *
Cudd_addSetToComplex(
  DdManager * dd,
  DdNode * f)
{

    if (cuddIsConstant(f)) {
      CUDD_VALUE_TYPE value;
      value.is_complex_assigned = 1;

      mpfr_init(value.real); mpfr_init(value.imag);
      mpfr_set_si(value.real, cuddV(f).val, RND_TYPE);
      mpfr_set_si(value.imag, 0, RND_TYPE);
      DdNode *res = cuddUniqueConst(dd,value);
      mpfr_clear(value.real); mpfr_clear(value.imag);
      return(res);
    }
    return(NULL);

} /* end of Cudd_addSimonsRemoveMinusOne */

/**
  @brief Converts -1 to 0
*/
DdNode *
Cudd_addSimonsRemoveMinusOne(
  DdManager * dd,
  DdNode * f)
{

    if (cuddIsConstant(f)) {
      CUDD_VALUE_TYPE value;
      value.is_complex_assigned = 0;
      mpfr_init(value.real);
      mpfr_init(value.imag);
      if (mpfr_cmp_si(cuddV(f).real, -1) == 0){
        mpfr_set_si(value.real, 0, RND_TYPE); 
        mpfr_set_si(value.imag, 0, RND_TYPE); 
      }else{
        mpfr_set(value.real, cuddV(f).real, RND_TYPE);
        mpfr_set(value.imag, cuddV(f).imag, RND_TYPE);
      }
      DdNode *res = cuddUniqueConst(dd,value);
      mpfr_clear(value.real);
      mpfr_clear(value.imag);
      return(res);
    }
    return(NULL);

} /* end of Cudd_addSimonsRemoveMinusOne */


/**
  @brief Converts 1 to 0 and -1 to 1
*/
DdNode *
Cudd_addSimonsRemoveOne(
  DdManager * dd,
  DdNode * f)
{
    if (cuddIsConstant(f)) {
      CUDD_VALUE_TYPE value;
      value.is_complex_assigned = 0;
      mpfr_init(value.real);
      mpfr_init(value.imag);
      if (mpfr_cmp_si(cuddV(f).real, -1) == 0){
        mpfr_set_si(value.real, 1, RND_TYPE); 
        mpfr_set_si(value.imag, 0, RND_TYPE); 
      }
      else if (mpfr_cmp_si(cuddV(f).real, 1) == 0){
        mpfr_set_si(value.real, 0, RND_TYPE); 
        mpfr_set_si(value.imag, 0, RND_TYPE); 
      }
      else{
        mpfr_set(value.real, cuddV(f).real, RND_TYPE);
        mpfr_set(value.imag, cuddV(f).imag, RND_TYPE);
      }
      DdNode *res = cuddUniqueConst(dd,value);
      mpfr_clear(value.real);
      mpfr_clear(value.imag);
      return(res);
    }
    return(NULL);

} /* end of Cudd_addSimonsRemoveOne */

/**
  @brief Squares the terminal values
*/
DdNode *
Cudd_addSquareTerminalValues(
  DdManager * dd,
  DdNode * f)
{
    if (cuddIsConstant(f)) {
      CUDD_VALUE_TYPE value;
      value.is_complex_assigned = 0;
      mpfr_init_set(value.real, cuddV(f).real, RND_TYPE);
      mpfr_init_set(value.imag, cuddV(f).imag, RND_TYPE);
      mpfr_mul(value.real, value.real, value.real, RND_TYPE);
      mpfr_mul(value.imag, value.imag, value.imag, RND_TYPE);
      mpfr_sub(value.real, value.real, value.imag, RND_TYPE);
      mpfr_set_si(value.imag, 0, RND_TYPE);
      DdNode *res = cuddUniqueConst(dd,value);
      mpfr_clear(value.real);
      mpfr_clear(value.imag);
      return(res);
    }
    return(NULL);

} /* end of Cudd_addSimonsRemoveOne */


/**
  @brief f if f&ge;g; 0 if f&lt;g.

  @details Threshold operator for Apply (f if f &ge;g; 0 if f&lt;g).

  @return NULL if not a terminal case; f op g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addThreshold(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == G || F == DD_PLUS_INFINITY(dd)) return(F);
    if (cuddIsConstant(F) && cuddIsConstant(G)) {
	// if (cuddV(F) >= cuddV(G)) {

      mpfr_t f_val, g_val;
      mpfr_init(f_val); mpfr_init(g_val);
      mpfr_t tmp; mpfr_init(tmp);
      mpfr_mul(f_val, cuddV(F).real, cuddV(F).real, RND_TYPE);
      mpfr_mul(tmp, cuddV(F).imag, cuddV(F).imag, RND_TYPE);
      mpfr_sub(f_val, f_val, tmp, RND_TYPE);

      mpfr_mul(g_val, cuddV(G).real, cuddV(G).real, RND_TYPE);
      mpfr_mul(tmp, cuddV(G).imag, cuddV(G).imag, RND_TYPE);
      mpfr_sub(g_val, g_val, tmp, RND_TYPE);      
      mpfr_clear(tmp);


      if (mpfr_cmp(f_val,g_val) >= 0) {
        mpfr_clear(f_val); mpfr_clear(g_val);
        return(F);
      }
      else {
        mpfr_clear(f_val); mpfr_clear(g_val);
        return(DD_ZERO(dd));
      }
    }
    return(NULL);

} /* end of Cudd_addThreshold */


/**
  @brief This operator sets f to the value of g wherever g != 0.

  @return NULL if not a terminal case; f op g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addSetNZ(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == G) return(F);
    if (F == DD_ZERO(dd)) return(G);
    if (G == DD_ZERO(dd)) return(F);
    if (cuddIsConstant(G)) return(G);
    return(NULL);

} /* end of Cudd_addSetNZ */


/**
  @brief Integer and floating point division.

  @return NULL if not a terminal case; f / g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addDivide(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *res;
    DdNode *F, *G;
    CUDD_VALUE_TYPE value;
    value.is_complex_assigned = 1;
    mpfr_init(value.real);
    mpfr_init(value.imag);

    F = *f; G = *g;
    /* We would like to use F == G -> F/G == 1, but F and G may
    ** contain zeroes. */
    if (F == DD_ZERO(dd)){mpfr_clear(value.real);mpfr_clear(value.imag); return(DD_ZERO(dd));}
    if (G == DD_ONE(dd)) {mpfr_clear(value.real);mpfr_clear(value.imag);return(F);}
    if (cuddIsConstant(F) && cuddIsConstant(G)) {
      mpfr_div(value.real, cuddV(F).real, cuddV(G).real, RND_TYPE);
      mpfr_div(value.imag, cuddV(F).imag, cuddV(G).imag, RND_TYPE);
	// value = cuddV(F)/cuddV(G);
	res = cuddUniqueConst(dd,value);
  mpfr_clear(value.real);
  mpfr_clear(value.imag);
	return(res);
    }
    mpfr_clear(value.real);
    mpfr_clear(value.imag);
    return(NULL);

} /* end of Cudd_addDivide */


/**
  @brief Integer and floating point subtraction.

  @return NULL if not a terminal case; f - g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addMinus(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *res;
    DdNode *F, *G;
    CUDD_VALUE_TYPE value;
    mpfr_init(value.real);
    mpfr_init(value.imag);
    value.is_complex_assigned = 1;

    F = *f; G = *g;
    if (F == G) {mpfr_clear(value.real);mpfr_clear(value.imag);return(DD_ZERO(dd));}
    if (F == DD_ZERO(dd)) {mpfr_clear(value.real);mpfr_clear(value.imag);return(cuddAddNegateRecur(dd,G));}
    if (G == DD_ZERO(dd)) {mpfr_clear(value.real);mpfr_clear(value.imag);return(F);}
    if (cuddIsConstant(F) && cuddIsConstant(G)) {
      mpfr_sub(value.real, cuddV(F).real, cuddV(G).real, RND_TYPE);
      mpfr_sub(value.imag, cuddV(F).imag, cuddV(G).imag, RND_TYPE);
      // mpfr_printf("minus: %.128Rf %.128Rf %.128Rf\n", value.real, cuddV(F).real, cuddV(G).real);
	// value = cuddV(F)-cuddV(G);
	res = cuddUniqueConst(dd,value);
  mpfr_clear(value.real);
  mpfr_clear(value.imag);
	return(res);
    }
    mpfr_clear(value.real);
    mpfr_clear(value.imag);
    return(NULL);

} /* end of Cudd_addMinus */


/**
  @brief Integer and floating point min.

  @details Integer and floating point min for Cudd_addApply.
  
  @return NULL if not a terminal case; min(f,g) otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addMinimum(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == DD_PLUS_INFINITY(dd)) return(G);
    if (G == DD_PLUS_INFINITY(dd)) return(F);
    if (F == G) return(F);
#if 0
    /* These special cases probably do not pay off. */
    if (F == DD_MINUS_INFINITY(dd)) return(F);
    if (G == DD_MINUS_INFINITY(dd)) return(G);
#endif
    if (cuddIsConstant(F) && cuddIsConstant(G)) {

      mpfr_t f_val, g_val;
      mpfr_init(f_val); mpfr_init(g_val);
      mpfr_t tmp; mpfr_init(tmp);
      mpfr_mul(f_val, cuddV(F).real, cuddV(F).real, RND_TYPE);
      mpfr_mul(tmp, cuddV(F).imag, cuddV(F).imag, RND_TYPE);
      mpfr_sub(f_val, f_val, tmp, RND_TYPE);

      mpfr_mul(g_val, cuddV(G).real, cuddV(G).real, RND_TYPE);
      mpfr_mul(tmp, cuddV(G).imag, cuddV(G).imag, RND_TYPE);
      mpfr_sub(g_val, g_val, tmp, RND_TYPE);      
      mpfr_clear(tmp);


      if (mpfr_cmp(f_val,g_val) <= 0) { // real
	// if (cuddV(F) <= cuddV(G)) {
        mpfr_clear(f_val); mpfr_clear(g_val);
	    return(F);
	} else{
    mpfr_clear(f_val); mpfr_clear(g_val);
	    return(G);
	}
    }
    if (F > G) { /* swap f and g */
	*f = G;
	*g = F;
    }
    return(NULL);

} /* end of Cudd_addMinimum */


/**
  @brief Integer and floating point max.

  @details Integer and floating point max for Cudd_addApply.

  @return NULL if not a terminal case; max(f,g) otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addMaximum(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == G) return(F);
    if (F == DD_MINUS_INFINITY(dd)) return(G);
    if (G == DD_MINUS_INFINITY(dd)) return(F);
#if 0
    /* These special cases probably do not pay off. */
    if (F == DD_PLUS_INFINITY(dd)) return(F);
    if (G == DD_PLUS_INFINITY(dd)) return(G);
#endif
    if (cuddIsConstant(F) && cuddIsConstant(G)) {

      mpfr_t f_val, g_val;
      mpfr_init(f_val); mpfr_init(g_val);
      mpfr_t tmp; mpfr_init(tmp);
      mpfr_mul(f_val, cuddV(F).real, cuddV(F).real, RND_TYPE);
      mpfr_mul(tmp, cuddV(F).imag, cuddV(F).imag, RND_TYPE);
      mpfr_sub(f_val, f_val, tmp, RND_TYPE);

      mpfr_mul(g_val, cuddV(G).real, cuddV(G).real, RND_TYPE);
      mpfr_mul(tmp, cuddV(G).imag, cuddV(G).imag, RND_TYPE);
      mpfr_sub(g_val, g_val, tmp, RND_TYPE);      
      mpfr_clear(tmp);

      if (mpfr_cmp(f_val,g_val) >= 0) {
        mpfr_clear(f_val); mpfr_clear(g_val);
	// if (cuddV(F) >= cuddV(G)) {
	    return(F);
	} else {
    mpfr_clear(f_val); mpfr_clear(g_val);
	    return(G);
	}
    }
    if (F > G) { /* swap f and g */
	*f = G;
	*g = F;
    }
    return(NULL);

} /* end of Cudd_addMaximum */


/**
  @brief Returns 1 if f &gt; g and 0 otherwise.

  @details Used in conjunction with Cudd_addApply.

  @return NULL if not a terminal case.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addOneZeroMaximum(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{

    if (*f == *g) return(DD_ZERO(dd));
    if (*g == DD_PLUS_INFINITY(dd))
	return DD_ZERO(dd);
    if (cuddIsConstant(*f) && cuddIsConstant(*g)) {
      mpfr_t f_val, g_val;
      mpfr_init(f_val); mpfr_init(g_val);
      mpfr_t tmp; mpfr_init(tmp);
      mpfr_mul(f_val, cuddV(*f).real, cuddV(*f).real, RND_TYPE);
      mpfr_mul(tmp, cuddV(*f).imag, cuddV(*f).imag, RND_TYPE);
      mpfr_sub(f_val, f_val, tmp, RND_TYPE);

      mpfr_mul(g_val, cuddV(*g).real, cuddV(*g).real, RND_TYPE);
      mpfr_mul(tmp, cuddV(*g).imag, cuddV(*g).imag, RND_TYPE);
      mpfr_sub(g_val, g_val, tmp, RND_TYPE);      
      mpfr_clear(tmp);
      if (mpfr_cmp(f_val, g_val) > 0) {
	// if (cuddV(*f) > cuddV(*g)) {
        mpfr_clear(f_val); mpfr_clear(g_val);
	    return(DD_ONE(dd));
	} else {
    mpfr_clear(f_val); mpfr_clear(g_val);
	    return(DD_ZERO(dd));
	}
    }

    return(NULL);

} /* end of Cudd_addOneZeroMaximum */


/**
  @brief Returns plusinfinity if f=g; returns min(f,g) if f!=g.

  @return NULL if not a terminal case; f op g otherwise, where f op g
  is plusinfinity if f=g; min(f,g) if f!=g.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addDiff(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == G) return(DD_PLUS_INFINITY(dd));
    if (F == DD_PLUS_INFINITY(dd)) return(G);
    if (G == DD_PLUS_INFINITY(dd)) return(F);
    if (cuddIsConstant(F) && cuddIsConstant(G)) {

      mpfr_t f_val, g_val;
      mpfr_init(f_val); mpfr_init(g_val);
      mpfr_t tmp; mpfr_init(tmp);
      mpfr_mul(f_val, cuddV(F).real, cuddV(F).real, RND_TYPE);
      mpfr_mul(tmp, cuddV(F).imag, cuddV(F).imag, RND_TYPE);
      mpfr_sub(f_val, f_val, tmp, RND_TYPE);

      mpfr_mul(g_val, cuddV(G).real, cuddV(G).real, RND_TYPE);
      mpfr_mul(tmp, cuddV(G).imag, cuddV(G).imag, RND_TYPE);
      mpfr_sub(g_val, g_val, tmp, RND_TYPE);      
      mpfr_clear(tmp);

      if (mpfr_cmp(f_val, g_val) != 0) {
	// if (cuddV(F) != cuddV(G)) {

        if (mpfr_cmp(f_val, g_val) < 0) {
	    // if (cuddV(F) < cuddV(G)) {
          mpfr_clear(f_val); mpfr_clear(g_val);
		return(F);
	    } else {
        mpfr_clear(f_val); mpfr_clear(g_val);
		return(G);
	    }
	} else {
    mpfr_clear(f_val); mpfr_clear(g_val);
	    return(DD_PLUS_INFINITY(dd));
	}
    }
    return(NULL);

} /* end of Cudd_addDiff */


/**
  @brief f if f==g; background if f!=g.

  @return NULL if not a terminal case; f op g otherwise, where f op g
  is f if f==g; background if f!=g.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addAgreement(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == G) return(F);
    if (F == dd->background) return(F);
    if (G == dd->background) return(G);
    if (cuddIsConstant(F) && cuddIsConstant(G)) return(dd->background);
    return(NULL);

} /* end of Cudd_addAgreement */


/**
  @brief Disjunction of two 0-1 ADDs.

  @return NULL if not a terminal case; f OR g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addOr(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == DD_ONE(dd) || G == DD_ONE(dd)) return(DD_ONE(dd));
    if (cuddIsConstant(F)) return(G);
    if (cuddIsConstant(G)) return(F);
    if (F == G) return(F);
    if (F > G) { /* swap f and g */
	*f = G;
	*g = F;
    }
    return(NULL);

} /* end of Cudd_addOr */


/**
  @brief NAND of two 0-1 ADDs.

  @return NULL if not a terminal case; f NAND g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addNand(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == DD_ZERO(dd) || G == DD_ZERO(dd)) return(DD_ONE(dd));
    if (cuddIsConstant(F) && cuddIsConstant(G)) return(DD_ZERO(dd));
    if (F > G) { /* swap f and g */
	*f = G;
	*g = F;
    }
    return(NULL);

} /* end of Cudd_addNand */


/**
  @brief NOR of two 0-1 ADDs.

  @return NULL if not a terminal case; f NOR g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addNor(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == DD_ONE(dd) || G == DD_ONE(dd)) return(DD_ZERO(dd));
    if (cuddIsConstant(F) && cuddIsConstant(G)) return(DD_ONE(dd));
    if (F > G) { /* swap f and g */
	*f = G;
	*g = F;
    }
    return(NULL);

} /* end of Cudd_addNor */


/**
  @brief XOR of two 0-1 ADDs.

  @return NULL if not a terminal case; f XOR g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addXor(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == G) return(DD_ZERO(dd));
    if (F == DD_ONE(dd) && G == DD_ZERO(dd)) return(DD_ONE(dd));
    if (G == DD_ONE(dd) && F == DD_ZERO(dd)) return(DD_ONE(dd));
    if (cuddIsConstant(F) && cuddIsConstant(G)) return(DD_ZERO(dd));
    if (F > G) { /* swap f and g */
	*f = G;
	*g = F;
    }
    return(NULL);

} /* end of Cudd_addXor */


/**
  @brief XNOR of two 0-1 ADDs.

  @return NULL if not a terminal case; f XNOR g otherwise.

  @sideeffect None

  @see Cudd_addApply

*/
DdNode *
Cudd_addXnor(
  DdManager * dd,
  DdNode ** f,
  DdNode ** g)
{
    DdNode *F, *G;

    F = *f; G = *g;
    if (F == G) return(DD_ONE(dd));
    if (F == DD_ONE(dd) && G == DD_ONE(dd)) return(DD_ONE(dd));
    if (G == DD_ZERO(dd) && F == DD_ZERO(dd)) return(DD_ONE(dd));
    if (cuddIsConstant(F) && cuddIsConstant(G)) return(DD_ZERO(dd));
    if (F > G) { /* swap f and g */
	*f = G;
	*g = F;
    }
    return(NULL);

} /* end of Cudd_addXnor */


/**
  @brief Applies op to the discriminants of f.

  @return a pointer to the result if succssful; NULL otherwise.

  @sideeffect None

  @see Cudd_addApply Cudd_addLog

*/
DdNode *
Cudd_addMonadicApply(
  DdManager * dd,
  DD_MAOP op,
  DdNode * f)
{
    DdNode *res;

    do {
	dd->reordered = 0;
	res = cuddAddMonadicApplyRecur(dd,op,f);
    } while (dd->reordered == 1);
    if (dd->errorCode == CUDD_TIMEOUT_EXPIRED && dd->timeoutHandler) {
        dd->timeoutHandler(dd, dd->tohArg);
    }

    return(res);

} /* end of Cudd_addMonadicApply */


/**
  @brief Performs the recursive step of Cudd_addMonadicApply.

  @return a pointer to the result if successful; NULL otherwise.

  @sideeffect None

  @see cuddAddApplyRecur

*/
DdNode *
cuddAddMonadicWithArgApplyRecur(
  DdManager * dd,
  DD_MAOP2 op,
  unsigned int arg,
  DdNode * f)
{
    DdNode *res, *ft, *fe, *T, *E;
    unsigned int index;

    /* Check terminal cases. */
    statLine(dd);
    res = (*op)(dd,f, arg);
    if (res != NULL) return(res);

    /* Check cache. */
    res = cuddCacheLookup12(dd,op,f);
    if (res != NULL) return(res);

    checkWhetherToGiveUp(dd);

    /* Recursive step. */
    index = f->index;
    ft = cuddT(f);
    fe = cuddE(f);

    T = cuddAddMonadicWithArgApplyRecur(dd,op,arg, ft);
    if (T == NULL) return(NULL);
    cuddRef(T);

    E = cuddAddMonadicWithArgApplyRecur(dd,op, arg, fe);
    if (E == NULL) {
  Cudd_RecursiveDeref(dd,T);
  return(NULL);
    }
    cuddRef(E);

    res = (T == E) ? T : cuddUniqueInter(dd,(int)index,T,E);
    if (res == NULL) {
  Cudd_RecursiveDeref(dd, T);
  Cudd_RecursiveDeref(dd, E);
  return(NULL);
    }
    cuddDeref(T);
    cuddDeref(E);

    /* Store result. */
    cuddCacheInsert12(dd,op,f,res);

    return(res);

} /* end of cuddAddMonadicApplyRecur */

/**
  @brief Applies op to the discriminants of f.

  @return a pointer to the result if succssful; NULL otherwise.

  @sideeffect None

  @see Cudd_addApply Cudd_addLog

*/
DdNode *
Cudd_addMonadicWithArgApply(
  DdManager * dd,
  DD_MAOP2 op,
  unsigned int arg,
  DdNode * f)
{
    DdNode *res;

    do {
  dd->reordered = 0;
  res = cuddAddMonadicWithArgApplyRecur(dd,op, arg, f);
    } while (dd->reordered == 1);
    if (dd->errorCode == CUDD_TIMEOUT_EXPIRED && dd->timeoutHandler) {
        dd->timeoutHandler(dd, dd->tohArg);
    }

    return(res);

} /* end of Cudd_addMonadicApply */


/**
  @brief Natural logarithm of an %ADD.

  @details The discriminants of f must be positive double's.

  @return NULL if not a terminal case; log(f) otherwise.

  @sideeffect None

  @see Cudd_addMonadicApply

*/
DdNode *
Cudd_addLog(
  DdManager * dd,
  DdNode * f)
{
    if (cuddIsConstant(f)) {
	CUDD_VALUE_TYPE value;
  value.is_complex_assigned = 1;
  mpfr_init(value.real); mpfr_init(value.imag);
  mpfr_log(value.real, cuddV(f).real, RND_TYPE);
  mpfr_log(value.imag, cuddV(f).imag, RND_TYPE);
	DdNode *res = cuddUniqueConst(dd,value);
  mpfr_clear(value.real);
  mpfr_clear(value.imag);
	return(res);
    }
    return(NULL);

} /* end of Cudd_addLog */


/*---------------------------------------------------------------------------*/
/* Definition of internal functions                                          */
/*---------------------------------------------------------------------------*/


/**
  @brief Performs the recursive step of Cudd_addApply.

  @return a pointer to the result if successful; NULL otherwise.

  @sideeffect None

  @see cuddAddMonadicApplyRecur

*/
DdNode *
cuddAddApplyRecur(
  DdManager * dd,
  DD_AOP op,
  DdNode * f,
  DdNode * g)
{
    DdNode *res,
	   *fv, *fvn, *gv, *gvn,
	   *T, *E;
    int ford, gord;
    unsigned int index;
    DD_CTFP cacheOp;

    /* Check terminal cases. Op may swap f and g to increase the
     * cache hit rate.
     */
    statLine(dd);
    // printf("rec %g %g\n", Cudd_V(f), Cudd_V(g));
    res = (*op)(dd,&f,&g);
    if (res != NULL) return(res);

    /* Check cache. */
    cacheOp = (DD_CTFP) op;
    res = cuddCacheLookup2(dd,cacheOp,f,g);
    if (res != NULL) return(res);

    checkWhetherToGiveUp(dd);

    /* Recursive step. */
    ford = cuddI(dd,f->index);
    gord = cuddI(dd,g->index);
    if (ford <= gord) {
	index = f->index;
	fv = cuddT(f);
	fvn = cuddE(f);
    } else {
	index = g->index;
	fv = fvn = f;
    }
    if (gord <= ford) {
	gv = cuddT(g);
	gvn = cuddE(g);
    } else {
	gv = gvn = g;
    }
    // printf("T\n");
    // mpfr_printf("%.16Rf %.16Rf\n", Cudd_V(fv).real, Cudd_V(gv).real);
    // mpfr_printf("%g %g\n", Cudd_V(f), Cudd_V(g));
    T = cuddAddApplyRecur(dd,op,fv,gv);
    if (T == NULL) return(NULL);
    cuddRef(T);

    // printf("E\n");
    // mpfr_printf("%.16Rf %.16Rf\n", Cudd_V(fvn).real, Cudd_V(gvn).real);
    // mpfr_printf("%g %g\n", Cudd_V(f), Cudd_V(g));
    E = cuddAddApplyRecur(dd,op,fvn,gvn);
    if (E == NULL) {
	Cudd_RecursiveDeref(dd,T);
	return(NULL);
    }
    cuddRef(E);

    res = (T == E) ? T : cuddUniqueInter(dd,(int)index,T,E);
    if (res == NULL) {
	Cudd_RecursiveDeref(dd, T);
	Cudd_RecursiveDeref(dd, E);
	return(NULL);
    }
    cuddDeref(T);
    cuddDeref(E);

    /* Store result. */
    cuddCacheInsert2(dd,cacheOp,f,g,res);

    return(res);

} /* end of cuddAddApplyRecur */


/**
  @brief Performs the recursive step of Cudd_addMonadicApply.

  @return a pointer to the result if successful; NULL otherwise.

  @sideeffect None

  @see cuddAddApplyRecur

*/
DdNode *
cuddAddMonadicApplyRecur(
  DdManager * dd,
  DD_MAOP op,
  DdNode * f)
{
    DdNode *res, *ft, *fe, *T, *E;
    unsigned int index;

    /* Check terminal cases. */
    statLine(dd);
    res = (*op)(dd,f);
    if (res != NULL) return(res);

    /* Check cache. */
    res = cuddCacheLookup1(dd,op,f);
    if (res != NULL) return(res);

    checkWhetherToGiveUp(dd);

    /* Recursive step. */
    index = f->index;
    ft = cuddT(f);
    fe = cuddE(f);

    T = cuddAddMonadicApplyRecur(dd,op,ft);
    if (T == NULL) return(NULL);
    cuddRef(T);

    E = cuddAddMonadicApplyRecur(dd,op,fe);
    if (E == NULL) {
	Cudd_RecursiveDeref(dd,T);
	return(NULL);
    }
    cuddRef(E);

    res = (T == E) ? T : cuddUniqueInter(dd,(int)index,T,E);
    if (res == NULL) {
	Cudd_RecursiveDeref(dd, T);
	Cudd_RecursiveDeref(dd, E);
	return(NULL);
    }
    cuddDeref(T);
    cuddDeref(E);

    /* Store result. */
    cuddCacheInsert1(dd,op,f,res);

    return(res);

} /* end of cuddAddMonadicApplyRecur */


/**
  @brief Performs the recursive step of Cudd_addPathCounts.

  @return NULL.

  @sideeffect None

  @see cuddAddApplyRecur

  typedef struct path_info {
    CUDD_VALUE_TYPE weight;
    mpfr_t path_count;
    long int index;
    long int l_index;
    unsigned int l_path_count;
    long int r_index;
    unsigned int r_path_count;
} path_info;

*/
DdNode*
cuddAddPathCountsRecur(
  DdManager * dd,
  DdNode * f,
  int period,
  unsigned int N)
{
    DdNode *ft, *fe, *T, *E, *res;
    unsigned int index;


    /* Check terminal cases. */
    statLine(dd);

    /* Check cache. */
    res = cuddCacheLookup1(dd,Cudd_addPathCounts_dup,f);
    if (res != NULL) return(res);

    if (cuddIsConstant(f)){
      path_info* p = (path_info *)malloc(sizeof(path_info));
      mpfr_init_set(p->weight.real, cuddV(f).real, RND_TYPE);
      mpfr_init_set_si(p->path_count, 1, RND_TYPE);
      p->index = 0;
      p->l_index = -1;
      p->r_index = -1;
      mpfr_init_set_si(p->l_path_count, 0, RND_TYPE);
      mpfr_init_set_si(p->r_path_count, 0, RND_TYPE);

      total_path_info* ps = (total_path_info *)malloc(sizeof(total_path_info));
      ps->info = p;
      ps->size = 1;

      CUDD_VALUE_TYPE value;
      mpfr_init_set(value.real, cuddV(f).real, RND_TYPE);
      DdNode *res = cuddUniqueConst(dd,value);
      res->numPaths = ps;
      mpfr_clear(value.real);

      cuddCacheInsert1(dd,Cudd_addPathCounts_dup,f,res);

      return(res);

    }

    checkWhetherToGiveUp(dd);

    /* Recursive step. */
    index = f->index;
    ft = cuddT(f);
    fe = cuddE(f);

    T = cuddAddPathCountsRecur(dd,ft, period, N);
    if (T == NULL) return(NULL);
    cuddRef(T);

    E = cuddAddPathCountsRecur(dd,fe, period, N);
    if (E == NULL) {
  Cudd_RecursiveDeref(dd,T);
  return(NULL);
    }
    cuddRef(E);

    res = (T == E) ? T : cuddUniqueInter(dd,(int)index,T,E);
    if (res == NULL) {
  Cudd_RecursiveDeref(dd, T);
  Cudd_RecursiveDeref(dd, E);
  return(NULL);
    }

    total_path_info* t_ps = T->numPaths;
    total_path_info* e_ps = E->numPaths;

    path_info* t_p = t_ps->info;
    path_info* e_p = e_ps->info;

    path_info* curr_node_path_info = (path_info *)malloc(sizeof(path_info)*(t_ps->size + e_ps->size));
    unsigned int size = 0, t_i = 0, e_i = 0;
    mpfr_t left_mul, right_mul;

    if (!cuddIsConstant(T)){
      mpfr_init_set_si(left_mul, ((T->index) - index)/period - 1, RND_TYPE);
      mpfr_exp2(left_mul, left_mul, RND_TYPE);
    } else{
      mpfr_init_set_si(left_mul, ((period * N) - index)/period - 1, RND_TYPE);
      mpfr_exp2(left_mul, left_mul, RND_TYPE);
    }

    if (!cuddIsConstant(E)){
      mpfr_init_set_si(right_mul, ((E->index) - index)/period - 1, RND_TYPE);
      mpfr_exp2(right_mul, right_mul, RND_TYPE);
    }else{
      mpfr_init_set_si(right_mul, ((period * N) - index)/period - 1, RND_TYPE);
      mpfr_exp2(right_mul, right_mul, RND_TYPE);
    }

    while (t_i < t_ps->size && e_i < e_ps->size){
      int cmp = mpfr_cmp(t_p[t_i].weight.real, e_p[e_i].weight.real);
      if (cmp == 0){

        path_info p;
        mpfr_init_set(p.weight.real, t_p[t_i].weight.real, RND_TYPE);
        mpfr_init_set(p.path_count, t_p[t_i].path_count, RND_TYPE);
        mpfr_mul(p.path_count, p.path_count, left_mul, RND_TYPE);
        mpfr_t tmp;
        mpfr_init_set(tmp, e_p[e_i].path_count, RND_TYPE);
        mpfr_mul(tmp, tmp, right_mul, RND_TYPE);
        mpfr_add(p.path_count, p.path_count, tmp , RND_TYPE);
        mpfr_clear(tmp);
        p.index = size;
        p.l_index = t_p[t_i].index;
        p.r_index = e_p[e_i].index;
        mpfr_init_set(p.l_path_count, t_p[t_i].path_count, RND_TYPE);
        mpfr_mul(p.l_path_count, p.l_path_count, left_mul, RND_TYPE);
        mpfr_init_set(p.r_path_count, e_p[e_i].path_count, RND_TYPE);
        mpfr_mul(p.r_path_count, p.r_path_count, right_mul, RND_TYPE);

        curr_node_path_info[size] = p;

        t_i++; e_i++;
      } else if (cmp < 0){
        path_info p;
        mpfr_init_set(p.weight.real, t_p[t_i].weight.real, RND_TYPE);
        mpfr_init_set(p.path_count, t_p[t_i].path_count, RND_TYPE);
        mpfr_mul(p.path_count, p.path_count, left_mul, RND_TYPE);
        p.index = size;
        p.l_index = t_p[t_i].index;
        p.r_index = -1;
        mpfr_init_set(p.l_path_count, t_p[t_i].path_count, RND_TYPE);
        mpfr_mul(p.l_path_count, p.l_path_count, left_mul, RND_TYPE);
        mpfr_init_set_si(p.r_path_count, 0, RND_TYPE);

        curr_node_path_info[size] = p;        
        t_i++;
      } else{
        path_info p;
        mpfr_init_set(p.weight.real, e_p[e_i].weight.real, RND_TYPE);
        mpfr_init_set(p.path_count, e_p[e_i].path_count, RND_TYPE);
        mpfr_mul(p.path_count, p.path_count, right_mul, RND_TYPE);
        p.index = size;
        p.l_index = -1;
        p.r_index = e_p[e_i].index;
        mpfr_init_set_si(p.l_path_count, 0, RND_TYPE);
        mpfr_init_set(p.r_path_count, e_p[e_i].path_count, RND_TYPE);
        mpfr_mul(p.r_path_count, p.r_path_count, right_mul, RND_TYPE);

        curr_node_path_info[size] = p;
        e_i++;
      }
      size++;
    }

    while (t_i < t_ps->size){
      path_info p;
        mpfr_init_set(p.weight.real, t_p[t_i].weight.real, RND_TYPE);
        mpfr_init_set(p.path_count, t_p[t_i].path_count, RND_TYPE);
        mpfr_mul(p.path_count, p.path_count, left_mul, RND_TYPE);
        p.index = size;
        p.l_index = t_p[t_i].index;
        p.r_index = -1;
        mpfr_init_set(p.l_path_count, t_p[t_i].path_count, RND_TYPE);
        mpfr_mul(p.l_path_count, p.l_path_count, left_mul, RND_TYPE);
        mpfr_init_set_si(p.r_path_count, 0, RND_TYPE);

        curr_node_path_info[size] = p;             
        t_i++;
        size++;
    }

    while (e_i < e_ps->size){
      path_info p;
        mpfr_init_set(p.weight.real, e_p[e_i].weight.real, RND_TYPE);
        mpfr_init_set(p.path_count, e_p[e_i].path_count, RND_TYPE);
        mpfr_mul(p.path_count, p.path_count, right_mul, RND_TYPE);
        p.index = size;
        p.l_index = -1;
        p.r_index = e_p[e_i].index;
        mpfr_init_set_si(p.l_path_count, 0, RND_TYPE);
        mpfr_init_set(p.r_path_count, e_p[e_i].path_count, RND_TYPE);
        mpfr_mul(p.r_path_count, p.r_path_count, right_mul, RND_TYPE);

        curr_node_path_info[size] = p;
        e_i++;
        size++;
    }

    mpfr_clear(left_mul);
    mpfr_clear(right_mul);

    total_path_info* curr_node_total_path_info = (total_path_info *)malloc(sizeof(total_path_info));
    // if (t_ps->size + e_ps->size > size){
    //   for (unsigned int k = size; k < t_ps->size + e_ps->size; k++){
    //     free (curr_node_path_info + k);
    //   }
    // }
    curr_node_total_path_info->info = curr_node_path_info;
    curr_node_total_path_info->size = size;

    res->numPaths = curr_node_total_path_info;

    cuddDeref(T);
    cuddDeref(E);

    // printf("size: %d index: %d\n", size, index);

    /* Store result. */
    cuddCacheInsert1(dd,Cudd_addPathCounts_dup,f,res);

    return(res);

} /* end of cuddAddMonadicApplyRecur */


/**
  @brief Path Counting and caching of f.

  @return void.

  @sideeffect None

  @see Cudd_addApply Cudd_addLog

*/
DdNode* 
Cudd_addPathCounts(
  DdManager * dd,
  DdNode * f,
  int period,
  unsigned int N)
{
    DdNode *res;

    do {
  dd->reordered = 0;
  res =     cuddAddPathCountsRecur(dd,f, period, N);
    } while (dd->reordered == 1);
    if (dd->errorCode == CUDD_TIMEOUT_EXPIRED && dd->timeoutHandler) {
        dd->timeoutHandler(dd, dd->tohArg);
    }

    return(res);

} /* end of Cudd_addPathCounts */


DdNode* 
Cudd_addPathCounts_dup(
  DdManager * dd,
  DdNode * f)
{
    return f;

} /* end of Cudd_addPathCounts */

/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/
