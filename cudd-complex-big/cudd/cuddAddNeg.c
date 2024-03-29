/**
  @file

  @ingroup cudd

  @brief Function to compute the negation of an %ADD.

  @author Fabio Somenzi, Balakrishna Kumthekar

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
  @brief Computes the additive inverse of an %ADD.

  @return a pointer to the result if successful; NULL otherwise.

  @sideeffect None

  @see Cudd_addCmpl

*/
DdNode *
Cudd_addNegate(
  DdManager * dd,
  DdNode * f)
{
    DdNode *res;

    do {
	dd->reordered = 0;
	res = cuddAddNegateRecur(dd,f);
    } while (dd->reordered == 1);
    if (dd->errorCode == CUDD_TIMEOUT_EXPIRED && dd->timeoutHandler) {
        dd->timeoutHandler(dd, dd->tohArg);
    }
    return(res);

} /* end of Cudd_addNegate */


/**
  @brief Rounds off the discriminants of an %ADD.

  @details The discriminants are rounded off to N digits after the
  decimal.

  @return a pointer to the result %ADD if successful; NULL otherwise.

  @sideeffect None

*/
DdNode *
Cudd_addRoundOff(
  DdManager * dd,
  DdNode * f,
  int  N)
{
    DdNode *res;
    double trunc = pow(10.0,(double)N);

    do {
	dd->reordered = 0;
	res = cuddAddRoundOffRecur(dd,f,trunc);
    } while (dd->reordered == 1);
    if (dd->errorCode == CUDD_TIMEOUT_EXPIRED && dd->timeoutHandler) {
        dd->timeoutHandler(dd, dd->tohArg);
    }
    return(res);

} /* end of Cudd_addRoundOff */


/*---------------------------------------------------------------------------*/
/* Definition of internal functions                                          */
/*---------------------------------------------------------------------------*/


/**
  @brief Implements the recursive step of Cudd_addNegate.

  @return a pointer to the result.

  @sideeffect None

*/
DdNode *
cuddAddNegateRecur(
  DdManager * dd,
  DdNode * f)
{
    DdNode *res,
	    *fv, *fvn,
	    *T, *E;

    statLine(dd);
    /* Check terminal cases. */
    if (cuddIsConstant(f)) {
      CUDD_VALUE_TYPE t_val;
      mpfr_init_set(t_val.real, cuddV(f).real, RND_TYPE);
      mpfr_init_set(t_val.imag, cuddV(f).imag, RND_TYPE);
      mpfr_mul_si(t_val.real,t_val.real,-1, RND_TYPE);
      mpfr_mul_si(t_val.imag,t_val.imag,-1, RND_TYPE);
      res = cuddUniqueConst(dd,t_val);
      mpfr_clear(t_val.real);
      mpfr_clear(t_val.imag);
      return(res);
    }

    /* Check cache */
    res = cuddCacheLookup1(dd,Cudd_addNegate,f);
    if (res != NULL) return(res);

    checkWhetherToGiveUp(dd);

    /* Recursive Step */
    fv = cuddT(f);
    fvn = cuddE(f);
    T = cuddAddNegateRecur(dd,fv);
    if (T == NULL) return(NULL);
    cuddRef(T);

    E = cuddAddNegateRecur(dd,fvn);
    if (E == NULL) {
	Cudd_RecursiveDeref(dd,T);
	return(NULL);
    }
    cuddRef(E);
    res = (T == E) ? T : cuddUniqueInter(dd,(int)f->index,T,E);
    if (res == NULL) {
	Cudd_RecursiveDeref(dd, T);
	Cudd_RecursiveDeref(dd, E);
	return(NULL);
    }
    cuddDeref(T);
    cuddDeref(E);

    /* Store result. */
    cuddCacheInsert1(dd,Cudd_addNegate,f,res);

    return(res);

} /* end of cuddAddNegateRecur */


/**
  @brief Implements the recursive step of Cudd_addRoundOff.

  @return a pointer to the result.

  @sideeffect None

*/
DdNode *
cuddAddRoundOffRecur(
  DdManager * dd,
  DdNode * f,
  double  trunc)
{

    DdNode *res, *fv, *fvn, *T, *E;
    CUDD_VALUE_TYPE n;
    mpfr_init(n.real);
    mpfr_init(n.imag);
    DD_CTFP1 cacheOp;

    statLine(dd);
    if (cuddIsConstant(f)) {
      // n = ceil(cuddV(f)*trunc)/trunc;
      mpfr_mul_d(n.real, cuddV(f).real, trunc, RND_TYPE);
      mpfr_mul_d(n.imag, cuddV(f).imag, trunc, RND_TYPE);
      mpfr_ceil(n.real, n.real);
      mpfr_ceil(n.imag, n.imag);
      mpfr_div_d(n.real, n.real, trunc, RND_TYPE);
      mpfr_div_d(n.imag, n.imag, trunc, RND_TYPE);
      res = cuddUniqueConst(dd,n);
      mpfr_clear(n.real);
      mpfr_clear(n.imag);
      return(res);
    }
    cacheOp = (DD_CTFP1) Cudd_addRoundOff;
    res = cuddCacheLookup1(dd,cacheOp,f);
    if (res != NULL) {
      mpfr_clear(n.real);
      mpfr_clear(n.imag);
	return(res);
    }
    checkWhetherToGiveUp(dd);
    /* Recursive Step */
    fv = cuddT(f);
    fvn = cuddE(f);
    T = cuddAddRoundOffRecur(dd,fv,trunc);
    if (T == NULL) {
      mpfr_clear(n.real);
      mpfr_clear(n.imag);
      return(NULL);
    }
    cuddRef(T);
    E = cuddAddRoundOffRecur(dd,fvn,trunc);
    if (E == NULL) {
      Cudd_RecursiveDeref(dd,T);
      mpfr_clear(n.real);
      mpfr_clear(n.imag);
	    return(NULL);
    }
    cuddRef(E);
    res = (T == E) ? T : cuddUniqueInter(dd,(int)f->index,T,E);
    if (res == NULL) {
      Cudd_RecursiveDeref(dd,T);
      Cudd_RecursiveDeref(dd,E);
      mpfr_clear(n.real);
      mpfr_clear(n.imag);
      return(NULL);
    }
    cuddDeref(T);
    cuddDeref(E);

    /* Store result. */
    cuddCacheInsert1(dd,cacheOp,f,res);
    mpfr_clear(n.real);
    mpfr_clear(n.imag);
    return(res);

} /* end of cuddAddRoundOffRecur */

/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/
