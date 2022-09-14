#ifndef MATRIX1234_TOP_NODE_GUARD
#define MATRIX1234_TOP_NODE_GUARD

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
		extern void Matrix1234InitializerTop();

		extern CFLOBDDTopNodeIntRefPtr MkIdRelationInterleavedTop(unsigned int i); // Representation of identity relation
		extern CFLOBDDTopNodeIntRefPtr MkWalshInterleavedTop(unsigned int i); // Representation of Walsh matrix
		extern CFLOBDDTopNodeIntRefPtr MkInverseReedMullerInterleavedTop(unsigned int i); // Representation of Inverse Reed-Muller matrix

		// Matrix-related operations (on matrices with room for two extra vocabularies) ------------
		extern CFLOBDDTopNodeIntRefPtr MkWalshVoc13Top(unsigned int i); // Representation of Walsh matrix with room for two additional vocabularies
		extern CFLOBDDTopNodeIntRefPtr MkWalshVoc12Top(unsigned int i); // Representation of Walsh matrix with room for two additional vocabularies
		extern CFLOBDDTopNodeIntRefPtr MatrixShiftVocs13To24Top(CFLOBDDTopNodeIntRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeIntRefPtr MatrixShiftVocs12To34Top(CFLOBDDTopNodeIntRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeIntRefPtr MatrixShiftVoc43Top(CFLOBDDTopNodeIntRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeIntRefPtr MatrixShiftVoc42Top(CFLOBDDTopNodeIntRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeIntRefPtr MkDetensorConstraintInterleavedTop(unsigned int i);
		// extern CFLOBDDTopNodeIntRefPtr MatrixProjectVoc23Top(CFLOBDDTopNodeIntRefPtr n); // Vocabulary projection
		//extern CFLOBDDTopNodeIntRefPtr MatrixConvertVocs1234To12Top(CFLOBDDTopNodeIntRefPtr n);
		extern CFLOBDDTopNodeIntRefPtr MatrixTransposeTop(CFLOBDDTopNodeIntRefPtr n);

		// Subroutines for Discrete Fourier Transform
		extern CFLOBDDTopNodeIntRefPtr MkCFLOBDDMatrixEqVoc14Top(unsigned int i);
		extern CFLOBDDTopNodeIntRefPtr MkFourierDiagonalComponentTop(unsigned int i);
		extern CFLOBDDTopNodeIntRefPtr PromoteInterleavedTo12Top(CFLOBDDTopNodeIntRefPtr c);
		extern CFLOBDDTopNodeIntRefPtr Demote12ToInterleavedTop(CFLOBDDTopNodeIntRefPtr c);

		extern void MatrixPrintRowMajorTop(CFLOBDDTopNodeIntRefPtr n, std::ostream & out);
	}
}

#endif

