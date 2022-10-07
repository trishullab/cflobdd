#ifndef MATRIX1234_FOURIER_TOP_NODE_GUARD
#define MATRIX1234_FOURIER_TOP_NODE_GUARD

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
#include "matrix1234_fourier.h"
#include "return_map_T.h"
namespace CFL_OBDD {
	typedef ReturnMapBody<fourierSemiring> FourierReturnMapBody;
	typedef ReturnMapHandle<fourierSemiring> FourierReturnMapHandle;
}
#include "connectionT.h"
namespace CFL_OBDD {
	typedef ConnectionT<FourierReturnMapHandle> FourierConnection;
}

namespace CFL_OBDD {

	typedef CFLOBDDTopNodeT<fourierSemiring> CFLOBDDTopNodeFourier;
	typedef CFLOBDDTopNodeT<fourierSemiring>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeFourierRefPtr;

	namespace Matrix1234Fourier {

		// Initialization routine
		extern void Matrix1234InitializerTop();

		extern void MatrixPrintRowMajorInterleavedTop(CFLOBDDTopNodeFourierRefPtr n, std::ostream & out);
		extern void MatrixPrintRowMajorTop(CFLOBDDTopNodeFourierRefPtr n, std::ostream & out);


		extern CFLOBDDTopNodeFourierRefPtr MkIdRelationInterleavedTop(unsigned int i); // Representation of identity relation
		extern CFLOBDDTopNodeFourierRefPtr MkWalshInterleavedTop(unsigned int i); // Representation of Walsh matrix

		// Matrix-related operations (on matrices with room for two extra vocabularies) ------------
		extern CFLOBDDTopNodeFourierRefPtr MatrixShiftVocs13To24Top(CFLOBDDTopNodeFourierRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeFourierRefPtr MatrixShiftVocs12To34Top(CFLOBDDTopNodeFourierRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeFourierRefPtr MatrixShiftVoc43Top(CFLOBDDTopNodeFourierRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeFourierRefPtr MatrixShiftVoc42Top(CFLOBDDTopNodeFourierRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeFourierRefPtr MkDetensorConstraintInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeFourierRefPtr MatrixMultiplyV4TopNode(CFLOBDDTopNodeFourierRefPtr c1, CFLOBDDTopNodeFourierRefPtr c2);
		extern CFLOBDDTopNodeFourierRefPtr MatrixMultiplyV4WithInfoTopNode(CFLOBDDTopNodeFourierRefPtr c1, CFLOBDDTopNodeFourierRefPtr c2);

		// Subroutines for Discrete Fourier Transform
		extern CFLOBDDTopNodeFourierRefPtr MkCFLOBDDMatrixEqVoc14Top(unsigned int i);
		extern CFLOBDDTopNodeFourierRefPtr MkFourierDiagonalComponentTop(unsigned int i);
		extern CFLOBDDTopNodeFourierRefPtr MkFourierDiagonalComponent2VocsTop(unsigned int i);
		extern CFLOBDDTopNodeFourierRefPtr PromoteInterleavedTo12Top(CFLOBDDTopNodeFourierRefPtr c);
		extern CFLOBDDTopNodeFourierRefPtr Demote12ToInterleavedTop(CFLOBDDTopNodeFourierRefPtr c);
		extern CFLOBDDTopNodeFourierRefPtr MatrixShiftToAConnectionTop(CFLOBDDTopNodeFourierRefPtr c);
		extern CFLOBDDTopNodeFourierRefPtr MatrixShiftToBConnectionTop(CFLOBDDTopNodeFourierRefPtr c);

	}
}

#endif

