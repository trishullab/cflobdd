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

#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "matrix1234_node.h"
#include "matrix1234_int_top_node.h"
#include "matrix1234_int.h"
#include "vector_int.h"

namespace CFL_OBDD {

	namespace Matrix1234Int {

		void Matrix1234Initializer()
		{
			Matrix1234InitializerTop();
			return;
		}

		// Create representation of identity relation
		CFLOBDD MkIdRelationInterleaved(unsigned int i)
		{
			return CFLOBDD(MkIdRelationInterleavedTop(i));
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDD MkWalshInterleaved(unsigned int i)
		{
			assert(i <= CFLOBDD::maxLevel);
			return CFLOBDD(MkWalshInterleavedTop(i));
		}

		// Create representation of the Inverse Reed-Muller matrix IRM(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields IRM[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDD MkInverseReedMullerInterleaved(unsigned int i)
		{
			assert(i <= CFLOBDD::maxLevel);
			return CFLOBDD(MkInverseReedMullerInterleavedTop(i));
		}

		// ****************************************************************************
		// Matrix-related operations (on matrices with room for two extra vocabularies)
		// ****************************************************************************

		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for
		// two extra vocabularies
		CFLOBDD MkWalshVoc13(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD::maxLevel);
			return CFLOBDD(MkWalshVoc13Top(i));
		}

		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for
		// two extra vocabularies
		CFLOBDD MkWalshVoc12(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD::maxLevel);
			return CFLOBDD(MkWalshVoc12Top(i));
		}


		// Vocabulary shift in a matrix
		CFLOBDD MatrixShiftVocs13To24(CFLOBDD c)
		{
			return CFLOBDD(MatrixShiftVocs13To24Top(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD MatrixShiftVocs12To34(CFLOBDD c)
		{
			return CFLOBDD(MatrixShiftVocs12To34Top(c.root));
		}

		CFLOBDD PromoteInterleavedTo12(CFLOBDD c)
		{
			return CFLOBDD(PromoteInterleavedTo12Top(c.root));
		}

		extern CFLOBDD Demote12ToInterleaved(CFLOBDD c)
		{
			return CFLOBDD(Demote12ToInterleavedTop(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD MatrixShiftVoc43(CFLOBDD c)
		{
			return CFLOBDD(MatrixShiftVoc43Top(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD MatrixShiftVoc42(CFLOBDD c)
		{
			return CFLOBDD(MatrixShiftVoc42Top(c.root));
		}

		// Return the Kronecker product of two matrices
		CFLOBDD KroneckerProduct(CFLOBDD m1, CFLOBDD m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD m2_12To34 = MatrixShiftVocs12To34(m2);
			//m1.PrintYield(&std::cout);
			//std::cout << std::endl;
			//std::cout << std::endl;
			//m2_12To34.PrintYield(&std::cout);
			//std::cout << std::endl;
			//std::cout << std::endl;
			return m1 * m2_12To34;
		}

		//
		// MatrixDetensor
		//
		// Return the result of detensoring a 4-vocabulary matrix
		//
		// CFLOBDD MatrixDetensor(CFLOBDD k)
		// {
		// 	assert(k.root->level >= 2);

		// 	CFLOBDD e = MkDetensorConstraintInterleaved(k.root->level);
		// 	CFLOBDD m = k * e;
		// 	CFLOBDD p = MatrixProjectVoc23(m);
		// 	CFLOBDD ans = MatrixShiftVoc42(p);
		// 	return ans;
		// }

		//
		// MatrixMultiply
		//
		// Return the matrix product of m1 and m2
		//
		// CFLOBDD MatrixMultiply(CFLOBDD m1, CFLOBDD m2)
		// {
		// 	assert(m1.root->level == m2.root->level);
		// 	assert(m1.root->level >= 2);

		// 	CFLOBDD k = KroneckerProduct(m1, m2);
		// 	CFLOBDD ans = MatrixDetensor(k);
		// 	return ans;
		// }

		// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal.
		// (W,X,Y,Z) s.t. X==Y with interleaved variables
		CFLOBDD MkDetensorConstraintInterleaved(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD::maxLevel);
			return CFLOBDD(MkDetensorConstraintInterleavedTop(i));
		}

		// Create representation of a matrix in which vocabularies 1 and $ are constrained to be equal.
		// (W,X,Y,Z) s.t. W==Z with interleaved variables
		CFLOBDD MkCFLOBDDMatrixEqVoc14(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD::maxLevel);
			return CFLOBDD(MkCFLOBDDMatrixEqVoc14Top(i));
		}

		// // Vocabulary projection
		// CFLOBDD MatrixProjectVoc23(CFLOBDD c)
		// {
		// 	return CFLOBDD(MatrixProjectVoc23Top(c.root));
		// }

		// Convert Matrix to Matrix with additional vocabularies
		/*CFLOBDD MatrixConvertVocs1234To12(CFLOBDD c)
		{
			return CFLOBDD(MatrixConvertVocs1234To12Top(c.root));
		}*/

		CFLOBDD MatrixPadWithZeros(CFLOBDD c, unsigned int level)
		{
			assert(c.root->level < level);
			CFLOBDD m(c);
			for (unsigned int i = c.root->level; i < level; i++)
			{
				CFLOBDD tempV = VectorInt::MkBasisVector(m.root->level-1, 0);
				CFLOBDD tempM = VectorInt::VectorToMatrixInterleaved(tempV);
				tempM = Matrix1234Int::PromoteInterleavedTo12(tempM);
				CFLOBDD m_1234 = Matrix1234Int::PromoteInterleavedTo12(m);
				m = Matrix1234Int::KroneckerProduct(tempM, m_1234);
			}

			assert(m.root->level == level);
			return m;
		}

		CFLOBDD MatrixTranspose(CFLOBDD c)
		{
			return CFLOBDD(MatrixTransposeTop(c.root));
		}


		// Create representation of the diagonal matrix that is the second component of the Kronecker-product-based
		// construction of the Fourier matrix
		CFLOBDD MkFourierDiagonalComponent(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD::maxLevel);
			return CFLOBDD(MkFourierDiagonalComponentTop(i));
		}

		// Create representation of the matrix for the Discrete Fourier Transform (DFT)
		// The answer matrix is a level-i object that uses interleaved vocabularies
		// CFLOBDD MkFourierMatrixInterleaved(unsigned int i)
		// {
		// 	assert(1 <= i && i < (CFLOBDD::maxLevel - 1));  // Need two extra levels: +1 for Kronecker Product; +1 for matrix multiplication

		// 	if (i == 1) {
		// 		return MkWalshInterleaved(1);                                                          // Level 1, interleaved
		// 	}
		// 	else {  // Need to promote all matrices to vocs 1,2 to allow final matrix multiplication
		// 		CFLOBDD Id = PromoteInterleavedTo12(MkIdRelationInterleaved(i - 1));                   // Level i, Vocs 1,2

		// 		//std::cout << "---Id------------------------------------------------------------" << std::endl;
		// 		//std::cout << Id << std::endl;
		// 		//Id.PrintYield(&std::cout);
		// 		//std::cout << "----------------------------------------------------------------" << std::endl;

		// 		CFLOBDD FourierIMinusOne = PromoteInterleavedTo12(MkFourierMatrixInterleaved(i - 1));  // Level i, Vocs 1,2

		// 		//std::cout << "---FourierIMinusOne------------------------------------------------------------" << std::endl;
		// 		//std::cout << FourierIMinusOne << std::endl;
		// 		//FourierIMinusOne.PrintYield(&std::cout);
		// 		//std::cout << "----------------------------------------------------------------" << std::endl;

		// 		CFLOBDD term1 = PromoteInterleavedTo12(KroneckerProduct(FourierIMinusOne, Id));        // Level i+1, Vocs 1,2

		// 		//std::cout << "---term1------------------------------------------------------------" << std::endl;
		// 		//std::cout << term1 << std::endl;
		// 		//term1.PrintYield(&std::cout);
		// 		//std::cout << "----------------------------------------------------------------" << std::endl;

		// 		CFLOBDD term2 = PromoteInterleavedTo12(MkFourierDiagonalComponent(i));                 // Level i+1, Vocs 1,2

		// 		//std::cout << "---term2------------------------------------------------------------" << std::endl;
		// 		//std::cout << term2 << std::endl;
		// 		//term2.PrintYield(&std::cout);
		// 		//std::cout << "----------------------------------------------------------------" << std::endl;

		// 		CFLOBDD term3 = PromoteInterleavedTo12(KroneckerProduct(Id, FourierIMinusOne));        // Level i+1, Vocs 1,2

		// 		//std::cout << "---term3------------------------------------------------------------" << std::endl;
		// 		//std::cout << term3 << std::endl;
		// 		//term3.PrintYield(&std::cout);
		// 		//std::cout << "----------------------------------------------------------------" << std::endl;

		// 		CFLOBDD term4 = PromoteInterleavedTo12(MkCFLOBDDMatrixEqVoc14(i));                     // Level i+1, Vocs 1,2

		// 		//std::cout << "---term4------------------------------------------------------------" << std::endl;
		// 		//std::cout << term4 << std::endl;
		// 		//term4.PrintYield(&std::cout);
		// 		//std::cout << "----------------------------------------------------------------" << std::endl;

		// 		CFLOBDD temp12 = MatrixMultiply(term1, term2);                                          // Level i+1, Vocs 1,2

		// 		//std::cout << "---temp12------------------------------------------------------------" << std::endl;
		// 		//std::cout << temp12 << std::endl;
		// 		//temp12.PrintYield(&std::cout);
		// 		//std::cout << "----------------------------------------------------------------" << std::endl;

		// 		CFLOBDD temp34 = MatrixMultiply(term3, term4);                                          // Level i+1, Vocs 1,2

		// 		//std::cout << "---temp34------------------------------------------------------------" << std::endl;
		// 		//std::cout << temp34 << std::endl;
		// 		//temp34.PrintYield(&std::cout);
		// 		//std::cout << "----------------------------------------------------------------" << std::endl;

		// 		CFLOBDD temp1234 = MatrixMultiply(temp12, temp34);                                          // Level i+1, Vocs 1,2

		// 		//std::cout << "---temp1234------------------------------------------------------------" << std::endl;
		// 		//std::cout << temp1234 << std::endl;
		// 		//temp1234.PrintYield(&std::cout);
		// 		//std::cout << "----------------------------------------------------------------" << std::endl;

		// 		CFLOBDD tempInterleaved = Demote12ToInterleaved(temp1234);                                        // Level i, interleaved
		// 		return tempInterleaved;
		// 	}
		// }

		void MatrixPrintRowMajor(CFLOBDD c, std::ostream & out)
		{
			MatrixPrintRowMajorTop(c.root, out);
			return;
		}
	}
}

