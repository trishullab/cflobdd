/**
  @file

  @ingroup cudd

  @brief Functions to find maximum and minimum in an %ADD and to
  extract the i-th bit.

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
#include "cuddAbsVal.h"

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

static DdNode * addDoIthBit (DdManager *dd, DdNode *f, DdNode *index);

/** \endcond */

/*---------------------------------------------------------------------------*/
/* Definition of exported functions                                          */
/*---------------------------------------------------------------------------*/


// CUDD_VALUE_TYPE getAbsValue(CUDD_VALUE_TYPE v)
// {
//   CUDD_VALUE_TYPE v_abs;
//   mpfr_init(v_abs->real); mpfr_init_set_d(v_abs->imag, 0, RND_TYPE);
//   mpfr_mul(v_abs->real, v.real, v.real, RND_TYPE);
//   mpfr_t tmp; mpfr_init(tmp);
//   mpfr_mul(tmp, v.imag, v.imag, RND_TYPE);
//   mpfr_add(v_abs->real, v_abs->real, tmp, RND_TYPE);
//   mpfr_clear(tmp);
//   return v_abs;
// }

/**
  @brief Finds the maximum discriminant of f.

  @return a pointer to a constant %ADD.

  @sideeffect None

*/
DdNode *
Cudd_addFindMax(
  DdManager * dd,
  DdNode * f)
{
    DdNode *t, *e, *res;

    statLine(dd);
    if (cuddIsConstant(f)) {
	return(f);
    }

    res = cuddCacheLookup1(dd,Cudd_addFindMax,f);
    if (res != NULL) {
	return(res);
    }

    checkWhetherToGiveUp(dd);

    t  = Cudd_addFindMax(dd,cuddT(f));
    if (t == DD_PLUS_INFINITY(dd)) return(t);

    e  = Cudd_addFindMax(dd,cuddE(f));
    CUDD_VALUE_TYPE *t_abs, *e_abs;
    t_abs = (CUDD_VALUE_TYPE *)malloc(sizeof(CUDD_VALUE_TYPE));
    e_abs = (CUDD_VALUE_TYPE *)malloc(sizeof(CUDD_VALUE_TYPE));
    mpfr_init(t_abs->real); mpfr_init(t_abs->imag);
    getAbsValue(cuddV(t), t_abs);
    mpfr_init(e_abs->real); mpfr_init(e_abs->imag);
    getAbsValue(cuddV(e), e_abs);
    res = (mpfr_cmp(t_abs->real,e_abs->real) >= 0) ? t : e;
    mpfr_clear(t_abs->real); mpfr_clear(t_abs->imag);
    mpfr_clear(e_abs->real); mpfr_clear(e_abs->imag);
    free(t_abs); free(e_abs);
    cuddCacheInsert1(dd,Cudd_addFindMax,f,res);

    return(res);

} /* end of Cudd_addFindMax */


/**
  @brief Finds the minimum discriminant of f.

  @return a pointer to a constant %ADD.

  @sideeffect None

*/
DdNode *
Cudd_addFindMin(
  DdManager * dd,
  DdNode * f)
{
    DdNode *t, *e, *res;

    statLine(dd);
    if (cuddIsConstant(f)) {
	return(f);
    }

    res = cuddCacheLookup1(dd,Cudd_addFindMin,f);
    if (res != NULL) {
	return(res);
    }

    checkWhetherToGiveUp(dd);

    t  = Cudd_addFindMin(dd,cuddT(f));
    if (t == DD_MINUS_INFINITY(dd)) return(t);

    e  = Cudd_addFindMin(dd,cuddE(f));
    CUDD_VALUE_TYPE *t_abs, *e_abs;
    t_abs = (CUDD_VALUE_TYPE *)malloc(sizeof(CUDD_VALUE_TYPE));
    e_abs = (CUDD_VALUE_TYPE *)malloc(sizeof(CUDD_VALUE_TYPE));
    mpfr_init(t_abs->real); mpfr_init(t_abs->imag);
    getAbsValue(cuddV(t), t_abs);
    mpfr_init(e_abs->real); mpfr_init(e_abs->imag);
    getAbsValue(cuddV(e), e_abs); 
    res = ((mpfr_cmp(t_abs->real,e_abs->real) <= 0)) ? t : e;
    mpfr_clear(t_abs->real); mpfr_clear(t_abs->imag);
    mpfr_clear(e_abs->real); mpfr_clear(e_abs->imag);
    free(t_abs); free(e_abs);
    cuddCacheInsert1(dd,Cudd_addFindMin,f,res);

    return(res);

} /* end of Cudd_addFindMin */


/**
  @brief Extracts the i-th bit from an %ADD.

  @details Produces an %ADD from another %ADD by replacing all
  discriminants whose i-th bit is equal to 1 with 1, and all other
  discriminants with 0. The i-th bit refers to the integer
  representation of the leaf value. If the value has a fractional
  part, it is ignored. Repeated calls to this procedure allow one to
  transform an integer-valued %ADD into an array of ADDs, one for each
  bit of the leaf values.

  @return a pointer to the resulting %ADD if successful; NULL
  otherwise.

  @sideeffect None

  @see Cudd_addBddIthBit

*/
DdNode *
Cudd_addIthBit(
  DdManager * dd,
  DdNode * f,
  int  bit)
{
    DdNode *res;
    DdNode *index;
    
    /* Use a constant node to remember the bit, so that we can use the
    ** global cache.
    */
    CUDD_VALUE_TYPE t;
    mpfr_init_set_si(t.real, (long int)bit, RND_TYPE);
    mpfr_init(t.imag); 
    mpfr_set_zero(t.imag, RND_TYPE);
    index = cuddUniqueConst(dd, t);
    mpfr_clear(t.real);
    mpfr_clear(t.imag);
    if (index == NULL) return(NULL);
    cuddRef(index);

    do {
	dd->reordered = 0;
	res = addDoIthBit(dd, f, index);
    } while (dd->reordered == 1);

    if (res == NULL) {
	Cudd_RecursiveDeref(dd, index);
        if (dd->errorCode == CUDD_TIMEOUT_EXPIRED && dd->timeoutHandler) {
            dd->timeoutHandler(dd, dd->tohArg);
        }
	return(NULL);
    }
    cuddRef(res);
    Cudd_RecursiveDeref(dd, index);
    cuddDeref(res);
    return(res);

} /* end of Cudd_addIthBit */


/*---------------------------------------------------------------------------*/
/* Definition of internal functions                                          */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/


/**
  @brief Performs the recursive step for Cudd_addIthBit.

  @return a pointer to the %BDD if successful; NULL otherwise.

  @sideeffect None

*/
static DdNode *
addDoIthBit(
  DdManager * dd,
  DdNode * f,
  DdNode * index)
{
    DdNode *res, *T, *E;
    DdNode *fv, *fvn;
    unsigned int mask, value;
    DdHalfWord v;

    statLine(dd);
    /* Check terminal case. */
    if (cuddIsConstant(f)) {
      // TODO: Maybe change this?
	mask = 1U << ((unsigned int) mpfr_get_ui(cuddV(index).real, RND_TYPE));
	value = (unsigned) (int) mpfr_get_ui(cuddV(f).real, RND_TYPE);
	return((value & mask) == 0 ? DD_ZERO(dd) : DD_ONE(dd));
    }

    /* Check cache. */
    res = cuddCacheLookup2(dd,addDoIthBit,f,index);
    if (res != NULL) return(res);

    checkWhetherToGiveUp(dd);

    /* Recursive step. */
    v = f->index;
    fv = cuddT(f); fvn = cuddE(f);

    T = addDoIthBit(dd,fv,index);
    if (T == NULL) return(NULL);
    cuddRef(T);

    E = addDoIthBit(dd,fvn,index);
    if (E == NULL) {
	Cudd_RecursiveDeref(dd, T);
	return(NULL);
    }
    cuddRef(E);

    res = (T == E) ? T : cuddUniqueInter(dd,v,T,E);
    if (res == NULL) {
	Cudd_RecursiveDeref(dd, T);
	Cudd_RecursiveDeref(dd, E);
	return(NULL);
    }
    cuddDeref(T);
    cuddDeref(E);

    /* Store result. */
    cuddCacheInsert2(dd,addDoIthBit,f,index,res);

    return(res);

} /* end of addDoIthBit */

