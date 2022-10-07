#ifndef MATRIX1234_FOURIER_SEMIRING_GUARD
#define MATRIX1234_FOURIER_SEMIRING_GUARD

//
//    Copyright (c) 2017, 2018 Thomas W. Reps
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

#include <iostream>
#include <fstream>
#include "fourier_semiring.h"
#include "cflobdd_t.h"


namespace CFL_OBDD {

	typedef CFLOBDD_T<fourierSemiring> CFLOBDD_FOURIER;

	namespace Matrix1234Fourier {
		// Initialization routine
		extern void Matrix1234Initializer();

		extern void MatrixPrintRowMajorInterleaved(CFLOBDD_FOURIER c, std::ostream & out);
		extern void MatrixPrintRowMajor(CFLOBDD_FOURIER c, std::ostream & out);


		extern CFLOBDD_FOURIER MkIdRelationInterleaved(unsigned int i); // Representation of identity relation
		extern CFLOBDD_FOURIER MkWalshInterleaved(unsigned int i); // Representation of walsh relation
		extern CFLOBDD_FOURIER MatrixShiftVocs13To24(CFLOBDD_FOURIER c);                // Vocabulary shift in a matrix
		extern CFLOBDD_FOURIER MatrixShiftVocs12To34(CFLOBDD_FOURIER c);                // Vocabulary shift in a matrix
		extern CFLOBDD_FOURIER MatrixShiftVoc43(CFLOBDD_FOURIER c);                     // Vocabulary shift in a matrix
		extern CFLOBDD_FOURIER MatrixShiftVoc42(CFLOBDD_FOURIER c);                     // Vocabulary shift in a matrix
		extern CFLOBDD_FOURIER KroneckerProduct(CFLOBDD_FOURIER m1, CFLOBDD_FOURIER m2);        // Kronecker product (on matrices with room for two extra vocabularies)
		extern CFLOBDD_FOURIER MkDetensorConstraintInterleaved(unsigned int i); // Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal: (W,X,Y,Z) s.t. X==Y with interleaved variables
		extern CFLOBDD_FOURIER MatrixMultiplyV4(CFLOBDD_FOURIER m1, CFLOBDD_FOURIER m2);          // Matrix multiplication
		extern CFLOBDD_FOURIER MatrixMultiplyV4WithInfo(CFLOBDD_FOURIER m1, CFLOBDD_FOURIER m2);          // Matrix multiplication

		// Discrete Fourier Transform, and subroutines
		extern CFLOBDD_FOURIER MkFourierMatrixInterleavedV4(unsigned int i);      // Create representation of the DFT matrix
		extern CFLOBDD_FOURIER MkFourierMatrixInterleavedV4WithInfo(unsigned int i);      // Create representation of the DFT matrix
		extern CFLOBDD_FOURIER MkCFLOBDDMatrixEqVoc14(unsigned int i);          // Create representation of a matrix in which vocabularies 1 and 4 are constrained to be equal: (W,X,Y,Z) s.t. W==Z with interleaved variables
		extern CFLOBDD_FOURIER MkFourierDiagonalComponent(unsigned int i);
		extern CFLOBDD_FOURIER PromoteInterleavedTo12(CFLOBDD_FOURIER c);
		extern CFLOBDD_FOURIER Demote12ToInterleaved(CFLOBDD_FOURIER c);

	}
}

#endif

