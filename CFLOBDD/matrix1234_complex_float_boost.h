#ifndef MATRIX1234_COMPLEX_FLOAT_BOOST_GUARD
#define MATRIX1234_COMPLEX_FLOAT_BOOST_GUARD

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
#include <complex>
#include <boost/multiprecision/cpp_complex.hpp>
#include "cflobdd_t.h"
#include "fourier_semiring.h"

namespace mp = boost::multiprecision;

namespace CFL_OBDD {

	typedef mp::cpp_complex_100 BIG_COMPLEX_FLOAT;
	//typedef mp::number<mp::cpp_dec_float<200> > BIG_FLOAT;
	typedef CFLOBDD_T<BIG_COMPLEX_FLOAT> CFLOBDD_COMPLEX_BIG;
	typedef CFLOBDD_T<fourierSemiring> CFLOBDD_FOURIER;

	namespace Matrix1234ComplexFloatBoost {
		// Initialization routine
		extern void Matrix1234Initializer();

		extern void MatrixPrintRowMajorInterleaved(CFLOBDD_COMPLEX_BIG c, std::ostream & out);
		extern void MatrixPrintRowMajor(CFLOBDD_COMPLEX_BIG c, std::ostream & out);


		extern CFLOBDD_COMPLEX_BIG MkIdRelationInterleaved(unsigned int i); // Representation of identity relation
		extern CFLOBDD_COMPLEX_BIG MkWalshInterleaved(unsigned int i);              // Representation of Walsh matrix
		extern CFLOBDD_COMPLEX_BIG MkInverseReedMullerInterleaved(unsigned int i);  // Representation of Inverse Reed-Muller matrix
		extern CFLOBDD_COMPLEX_BIG MkExchangeInterleaved(unsigned int i); // Representation of exchange matrix
		extern CFLOBDD_COMPLEX_BIG MkNegationMatrixInterleaved(unsigned int i);
		extern CFLOBDD_COMPLEX_BIG MkPauliYMatrixInterleaved(unsigned int i);
		extern CFLOBDD_COMPLEX_BIG MkPauliZMatrixInterleaved(unsigned int i);
		extern CFLOBDD_COMPLEX_BIG MkSGateInterleaved(unsigned int i);
		extern CFLOBDD_COMPLEX_BIG MkPhaseShiftGateInterleaved(unsigned int i, double theta);

		extern CFLOBDD_COMPLEX_BIG MkRestrictMatrix(unsigned int level, std::string s); 

		// Matrix-related operations (on matrices with room for two extra vocabularies) ------------
		extern CFLOBDD_COMPLEX_BIG MkWalshVoc13(unsigned int i);                    // Create representation of Walsh matrix with room for two extra vocabularies
		extern CFLOBDD_COMPLEX_BIG MkWalshVoc12(unsigned int i);                    // Create representation of Walsh matrix with room for two extra vocabularies
		extern CFLOBDD_COMPLEX_BIG MatrixShiftVocs13To24(CFLOBDD_COMPLEX_BIG c);                // Vocabulary shift in a matrix
		extern CFLOBDD_COMPLEX_BIG MatrixShiftVocs12To34(CFLOBDD_COMPLEX_BIG c);                // Vocabulary shift in a matrix
		extern CFLOBDD_COMPLEX_BIG MatrixShiftVoc43(CFLOBDD_COMPLEX_BIG c);                     // Vocabulary shift in a matrix
		extern CFLOBDD_COMPLEX_BIG MatrixShiftVoc42(CFLOBDD_COMPLEX_BIG c);                     // Vocabulary shift in a matrix
		extern CFLOBDD_COMPLEX_BIG KroneckerProduct(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2);        // Kronecker product (on matrices with room for two extra vocabularies)
		extern CFLOBDD_COMPLEX_BIG MkDetensorConstraintInterleaved(unsigned int i); // Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal: (W,X,Y,Z) s.t. X==Y with interleaved variables
		// extern CFLOBDD_COMPLEX_BIG MatrixProjectVoc23(CFLOBDD_COMPLEX_BIG c);                   // Vocabulary projection
		extern CFLOBDD_COMPLEX_BIG MatrixDetensor(CFLOBDD_COMPLEX_BIG k);                       // Detensor of a 4-vocabulary matrix
		// extern CFLOBDD_COMPLEX_BIG MatrixMultiply(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2);          // Matrix multiplication
		extern CFLOBDD_COMPLEX_BIG MatrixMultiplyV4(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2);          // Matrix multiplication
		extern CFLOBDD_COMPLEX_BIG MatrixMultiplyV4WithInfo(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2);          // Matrix multiplication

		// Discrete Fourier Transform, and subroutines
		extern CFLOBDD_COMPLEX_BIG MkFourierMatrixInterleaved(unsigned int i);      // Create representation of the DFT matrix
		extern CFLOBDD_COMPLEX_BIG MkFourierMatrixInterleavedV4(unsigned int i);      // Create representation of the DFT matrix
		extern CFLOBDD_COMPLEX_BIG MkFourierMatrixInterleavedV4WithInfo(unsigned int i);      // Create representation of the DFT matrix
		extern CFLOBDD_COMPLEX_BIG MkCFLOBDDMatrixEqVoc14(unsigned int i);          // Create representation of a matrix in which vocabularies 1 and 4 are constrained to be equal: (W,X,Y,Z) s.t. W==Z with interleaved variables
		extern CFLOBDD_COMPLEX_BIG MkFourierDiagonalComponent(unsigned int i);
		extern CFLOBDD_COMPLEX_BIG PromoteInterleavedTo12(CFLOBDD_COMPLEX_BIG c);
		extern CFLOBDD_COMPLEX_BIG Demote12ToInterleaved(CFLOBDD_COMPLEX_BIG c);
		extern CFLOBDD_COMPLEX_BIG ConvertToComplex(CFLOBDD_FOURIER c);
		extern CFLOBDD_COMPLEX_BIG MkCNOTInterleaved(unsigned int i);
		extern CFLOBDD_COMPLEX_BIG MkCPGate(unsigned int i, long int c1, long int c2, double theta);
		extern CFLOBDD_COMPLEX_BIG MkSwapGate(unsigned int i, long int c1, long int c2);
		extern CFLOBDD_COMPLEX_BIG MkiSwapGate(unsigned int i, long int c1, long int c2);
		extern CFLOBDD_COMPLEX_BIG MkCSwapGate(unsigned int i, long int c1, long int x1, long int x2);
		extern CFLOBDD_COMPLEX_BIG MkCNOT(unsigned int level, unsigned int n, long int controller, long int controlled); // Representation of CNOT matrix with index1 as controller and index2 as controlled bits
		extern CFLOBDD_COMPLEX_BIG MkCCNOT(unsigned int level, unsigned int n, long int controller1, long int controller2, long int controlled);

		extern CFLOBDD_COMPLEX_BIG MatrixShiftToAConnection(CFLOBDD_COMPLEX_BIG c);
		extern CFLOBDD_COMPLEX_BIG MatrixShiftToBConnection(CFLOBDD_COMPLEX_BIG c);
		extern CFLOBDD_COMPLEX_BIG KroneckerProduct2Vocs(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2);

	}
}

#endif

