//
//    Copyright (c) 1999 Thomas W. Reps
//    All Rights Reserved.
//
//    This software is furnished under a license and may be used and
//    copied only in accordance with the terms of such license and the
//    inclusion of the above copyright notice.  This software or any
//    other copies thereof or any derivative works may not be provided
//    or otherwise made available to any other person.  Title to and
//    ownership of the software and any derivative works is retained
//    by Thomas W. Reps.
//
//    THIS IMPLEMENTATION MAY HAVE BUGS, SOME OF WHICH MAY HAVE SERIOUS
//    CONSEQUENCES.  THOMAS W. REPS PROVIDES THIS SOFTWARE IN ITS "AS IS"
//    CONDITION, AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
//    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//    AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
//    THOMAS W. REPS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <cassert>
#include "bool_op.h"

// BoolOp's ----------------------------------------------------------
BoolOp trueOp =          { {true,   true},    //  \a.\b.true
                           {true,   true}
                         };
BoolOp falseOp =         { {false, false},    //  \a.\b.false
                           {false, false}
                         };
BoolOp andOp =           { {false, false},    //  \a.\b.(a && b)
                           {false,  true}
                         };
BoolOp nandOp =          { {true,   true},    //  \a.\b.!(a && b)
                           {true,  false}
                         };
BoolOp orOp  =           { {false,  true},    //  \a.\b.(a || b)
                           {true,   true}
                         };
BoolOp norOp =           { {true,  false},    //  \a.\b.!(a || b)
                           {false, false}
                         };
BoolOp iffOp =           { {true,  false},    //  \a.\b.(a == b)
                           {false,  true}
                         };
BoolOp exclusiveOrOp =   { {false,  true},    //  \a.\b.(a != b)
                           {true,  false}
                         };
BoolOp impliesOp =       { {true,   true},    //  \a.\b.(!a || b)
                           {false,  true}
                         };
BoolOp minusOp =         { {false, false},    //  \a.\b.(a && !b)
                           {true,  false}
                         };
BoolOp quotientOp =      { {true,  false},    //  \a.\b.(a || !b)
                           {true,   true}
                         };
BoolOp notQuotientOp =   { {false,  true},    //  \a.\b.(!a && b)
                           {false, false}
                         };
BoolOp firstOp =         { {false, false},    //  \a.\b.a
                           {true,   true}
                         };
BoolOp notFirstOp =      { {true,   true},    //  \a.\b.!a
                           {false, false}
                         };
BoolOp secondOp =        { {false,  true},    //  \a.\b.b
                           {false,  true}
                         };
BoolOp notSecondOp =     { {true,  false},    //  \a.\b.!b
                           {true,  false}
                         };


int TrueFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return trueOp[a][b];
}

int FalseFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return falseOp[a][b];
}

int AndFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return andOp[a][b];
}

int NandFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return nandOp[a][b];
}

int OrFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return orOp[a][b];
}

int NorFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return norOp[a][b];
}

int IffFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return iffOp[a][b];
}

int ExclusiveOrFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return exclusiveOrOp[a][b];
}

int ImpliesFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return impliesOp[a][b];
}

int MinusFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return minusOp[a][b];
}

int QuotientFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return quotientOp[a][b];
}

int NotQuotientFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return notQuotientOp[a][b];
}

int FirstFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return firstOp[a][b];
}

int NotFirstFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return notFirstOp[a][b];
}

int SecondFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return notSecondOp[a][b];
}

int NotSecondFunc(int a, int b) {
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1);
	return notSecondOp[a][b];
}


// BoolOp3's ----------------------------------------------------------
// \a.\b.\c.(a && b) || (!a && c)
//                           ________________________________________________________
//                                         b
//                                false          true
//                           --------------------------------------------------------
//                                  c               c
//                            false   true    false   true 
//                            _______________________________________________________
BoolOp3 ifThenElseOp =   { { {false,  true}, {false,  true} },    //  | false |  a
                           { {false, false}, { true,  true} }     //  | true  |
                         };

// \a.\b.\c.(b && !a) || (c && !a) || (b && c)
//                          ________________________________________________________
//                                        b
//                               false          true
//                          --------------------------------------------------------
//                                 c               c
//                           false   true    false   true 
//                           _______________________________________________________
BoolOp3 negMajorityOp = { { {false,  true}, {true,  true} },    //  | false |  a
                          { {false, false}, {false, true} }     //  | true  |
                         };

int IfThenElseFunc(int a, int b, int c) {           //  \a.\b.\c.(a && b) || (!a && c)
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1);
	return ifThenElseOp[a][b][c];
}

int NegMajorityFunc(int a, int b, int c) {          //  \a.\b.\c.(b && !a) || (c && !a) || (b && c)
	assert(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1);
	return negMajorityOp[a][b][c];
}


