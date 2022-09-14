#ifndef BOOL_OP_GUARD
#define BOOL_OP_GUARD

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

// BoolOp's ----------------------------------------------------------
typedef bool BoolOp[2][2];
extern BoolOp trueOp;           //  \a.\b.true
extern BoolOp falseOp;          //  \a.\b.false
extern BoolOp andOp;            //  \a.\b.(a && b)
extern BoolOp nandOp;           //  \a.\b.!(a && b)
extern BoolOp orOp;             //  \a.\b.(a || b)
extern BoolOp norOp;            //  \a.\b.!(a || b)
extern BoolOp iffOp;            //  \a.\b.(a == b)
extern BoolOp exclusiveOrOp;    //  \a.\b.(a != b)
extern BoolOp impliesOp;        //  \a.\b.(!a || b)
extern BoolOp minusOp;          //  \a.\b.(a && !b)
extern BoolOp quotientOp;       //  \a.\b.(!b || a)
extern BoolOp notQuotientOp;    //  \a.\b.(b && !a)
extern BoolOp firstOp;          //  \a.\b.a
extern BoolOp notFirstOp;       //  \a.\b.!a
extern BoolOp secondOp;         //  \a.\b.b
extern BoolOp notSecondOp;      //  \a.\b.!b

extern int TrueFunc(int a, int b);           //  \a.\b.true
extern int FalseFunc(int a, int b);          //  \a.\b.false
extern int AndFunc(int a, int b);            //  \a.\b.(a && b)
extern int NandFunc(int a, int b);           //  \a.\b.!(a && b)
extern int OrFunc(int a, int b);             //  \a.\b.(a || b)
extern int NorFunc(int a, int b);            //  \a.\b.!(a || b)
extern int IffFunc(int a, int b);            //  \a.\b.(a == b)
extern int ExclusiveOrFunc(int a, int b);    //  \a.\b.(a != b)
extern int ImpliesFunc(int a, int b);        //  \a.\b.(!a || b)
extern int MinusFunc(int a, int b);          //  \a.\b.(a && !b)
extern int QuotientFunc(int a, int b);       //  \a.\b.(!b || a)
extern int NotQuotientFunc(int a, int b);    //  \a.\b.(b && !a)
extern int FirstFunc(int a, int b);          //  \a.\b.a
extern int NotFirstFunc(int a, int b);       //  \a.\b.!a
extern int SecondFunc(int a, int b);         //  \a.\b.b
extern int NotSecondFunc(int a, int b);      //  \a.\b.!b

// BoolOp3's ----------------------------------------------------------
typedef bool BoolOp3[2][2][2];
extern BoolOp3 ifThenElseOp;           //  \a.\b.\c.(a && b) || (!a && c)
extern BoolOp3 negMajorityOp;          //  \a.\b.\c.(b && !a) || (c && !a) || (b && c)

extern int IfThenElseFunc(int a, int b, int c);           //  \a.\b.\c.(a && b) || (!a && c)
extern int NegMajorityFunc(int a, int b, int c);          //  \a.\b.\c.(b && !a) || (c && !a) || (b && c)

// Some polymorphic non-Boolean functions ------------------------------------------------
template <typename T>
T PlusFunc(T a, T b) {
	return a + b;
}

template <typename T>
T TimesFunc(T a, T b) {
	return a * b;
}

#endif
