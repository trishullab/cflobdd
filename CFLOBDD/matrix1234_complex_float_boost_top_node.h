#ifndef MATRIX1234_COMPLEX_FLOAT_BOOST_TOP_NODE_GUARD
#define MATRIX1234_COMPLEX_FLOAT_BOOST_TOP_NODE_GUARD

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
#include "matrix1234_complex_float_boost.h"
#include "return_map_T.h"
namespace CFL_OBDD {
	typedef ReturnMapBody<BIG_COMPLEX_FLOAT> ComplexFloatBoostReturnMapBody;
	typedef ReturnMapHandle<BIG_COMPLEX_FLOAT> ComplexFloatBoostReturnMapHandle;
}
#include "connectionT.h"
namespace CFL_OBDD {
	typedef ConnectionT<ComplexFloatBoostReturnMapHandle> ComplexFloatBoostConnection;
}

namespace CFL_OBDD {

	typedef CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT> CFLOBDDTopNodeComplexFloatBoost;
	typedef CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeComplexFloatBoostRefPtr;
	typedef CFLOBDDTopNodeT<fourierSemiring>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeFourierRefPtr;

	namespace Matrix1234ComplexFloatBoost {

		// Initialization routine
		extern void Matrix1234InitializerTop();

		extern void MatrixPrintRowMajorInterleavedTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out);
		extern void MatrixPrintRowMajorTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out);


		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkIdRelationInterleavedTop(unsigned int i); // Representation of identity relation
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkWalshInterleavedTop(unsigned int i); // Representation of Walsh matrix
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkInverseReedMullerInterleavedTop(unsigned int i); // Representation of Inverse Reed-Muller matrix
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkNegationMatrixInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkPauliYMatrixInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkPauliZMatrixInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkSGateInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkPhaseShiftGateInterleavedTop(unsigned int i, double theta);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkExchangeInterleavedTop(unsigned int i); // Representation of exchange matrix

		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkRestrictTop(unsigned int level, std::string s);

		// Matrix-related operations (on matrices with room for two extra vocabularies) ------------
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkWalshVoc13Top(unsigned int i); // Representation of Walsh matrix with room for two additional vocabularies
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkWalshVoc12Top(unsigned int i); // Representation of Walsh matrix with room for two additional vocabularies
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftVocs13To24Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftVocs12To34Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftVoc43Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftVoc42Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkDetensorConstraintInterleavedTop(unsigned int i);
		// extern CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixProjectVoc23Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n); // Vocabulary projection
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixMultiplyV4TopNode(CFLOBDDTopNodeComplexFloatBoostRefPtr c1, CFLOBDDTopNodeComplexFloatBoostRefPtr c2);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixMultiplyV4WithInfoTopNode(CFLOBDDTopNodeComplexFloatBoostRefPtr c1, CFLOBDDTopNodeComplexFloatBoostRefPtr c2);

		// Subroutines for Discrete Fourier Transform
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkCFLOBDDMatrixEqVoc14Top(unsigned int i);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkFourierDiagonalComponentTop(unsigned int i);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr PromoteInterleavedTo12Top(CFLOBDDTopNodeComplexFloatBoostRefPtr c);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr Demote12ToInterleavedTop(CFLOBDDTopNodeComplexFloatBoostRefPtr c);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr ConvertToComplexTop(CFLOBDDTopNodeFourierRefPtr c);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkCNOTInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkCPGateTop(unsigned int level, long int i, long int j, double theta);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkSwapGateTop(unsigned int level, long int i, long int j);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkiSwapGateTop(unsigned int level, long int i, long int j);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkCSwapGateTop(unsigned int level, long int c, long int i, long int j);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkCCNOTTopNode(unsigned int level, unsigned int n, long int controller1, long int controller2, long int controlled);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkMCXTopNode(unsigned int level, unsigned int n, std::vector<long int>& controllers, long int controlled);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkCCPTopNode(unsigned int level, unsigned int n, long int controller1, long int controller2, long int controlled, double theta);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkSXGateTop(unsigned int i);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkSYGateTop(unsigned int i);

		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftToAConnectionTop(CFLOBDDTopNodeComplexFloatBoostRefPtr c);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftToBConnectionTop(CFLOBDDTopNodeComplexFloatBoostRefPtr c);

	}
}

#endif

