#ifndef MATRIX1234_INT_GUARD
#define MATRIX1234_INT_GUARD

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
#include "cflobdd_int.h"

namespace CFL_OBDD {

	namespace Matrix1234Int {

		// Initialization routine
		extern void Matrix1234Initializer();

		extern CFLOBDD MkIdRelationInterleaved(unsigned int i);         // Representation of identity relation
		extern CFLOBDD MkWalshInterleaved(unsigned int i);              // Representation of Walsh matrix
		extern CFLOBDD MkInverseReedMullerInterleaved(unsigned int i);  // Representation of Inverse Reed-Muller matrix

		// Matrix-related operations (on matrices with room for two extra vocabularies) ------------
		extern CFLOBDD MkWalshVoc13(unsigned int i);                    // Create representation of Walsh matrix with room for two extra vocabularies
		extern CFLOBDD MkWalshVoc12(unsigned int i);                    // Create representation of Walsh matrix with room for two extra vocabularies
		extern CFLOBDD MatrixShiftVocs13To24(CFLOBDD c);                // Vocabulary shift in a matrix
		extern CFLOBDD MatrixShiftVocs12To34(CFLOBDD c);                // Vocabulary shift in a matrix
		extern CFLOBDD MatrixShiftVoc43(CFLOBDD c);                     // Vocabulary shift in a matrix
		extern CFLOBDD MatrixShiftVoc42(CFLOBDD c);                     // Vocabulary shift in a matrix
		extern CFLOBDD KroneckerProduct(CFLOBDD m1, CFLOBDD m2);        // Kronecker product (on matrices with room for two extra vocabularies)
		extern CFLOBDD MkDetensorConstraintInterleaved(unsigned int i); // Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal: (W,X,Y,Z) s.t. X==Y with interleaved variables
		// extern CFLOBDD MatrixProjectVoc23(CFLOBDD c);                   // Vocabulary projection
		extern CFLOBDD MatrixDetensor(CFLOBDD k);                       // Detensor of a 4-vocabulary matrix
		extern CFLOBDD MatrixMultiply(CFLOBDD m1, CFLOBDD m2);          // Matrix multiplication
		//extern CFLOBDD MatrixConvertVocs1234To12(CFLOBDD c);			// Convert Matrix into Matrix with additional vocabularies
		extern CFLOBDD MatrixPadWithZeros(CFLOBDD c, unsigned int level);
		extern CFLOBDD MatrixTranspose(CFLOBDD c);

		// Discrete Fourier Transform, and subroutines
		extern CFLOBDD MkFourierMatrixInterleaved(unsigned int i);      // Create representation of the DFT matrix
		extern CFLOBDD MkCFLOBDDMatrixEqVoc14(unsigned int i);          // Create representation of a matrix in which vocabularies 1 and 4 are constrained to be equal: (W,X,Y,Z) s.t. W==Z with interleaved variables
		extern CFLOBDD MkFourierDiagonalComponent(unsigned int i);
		extern CFLOBDD PromoteInterleavedTo12(CFLOBDD c);
		extern CFLOBDD Demote12ToInterleaved(CFLOBDD c);

		extern void MatrixPrintRowMajor(CFLOBDD c, std::ostream & out);
	}
}

#endif

