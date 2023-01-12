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
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "matrix1234_node.h"
#include "matrix1234_complex_float_boost_top_node.h"
#include "matrix1234_complex_float_boost.h"

namespace CFL_OBDD {

	namespace Matrix1234ComplexFloatBoost {

		void Matrix1234Initializer()
		{
			Matrix1234InitializerTop();
			return;
		}

		void MatrixPrintRowMajorInterleaved(CFLOBDD_COMPLEX_BIG c, std::ostream & out)
		{
			MatrixPrintRowMajorInterleavedTop(c.root, out);
			return;
		}

		void MatrixPrintRowMajor(CFLOBDD_COMPLEX_BIG c, std::ostream & out)
		{
			MatrixPrintRowMajorTop(c.root, out);
			return;
		}

		// Create representation of identity relation
		CFLOBDD_COMPLEX_BIG MkIdRelationInterleaved(unsigned int i)
		{
			return CFLOBDD_COMPLEX_BIG(MkIdRelationInterleavedTop(i));
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDD_COMPLEX_BIG MkWalshInterleaved(unsigned int i)
		{
			assert(i <= CFLOBDD_COMPLEX_BIG::maxLevel);
			return CFLOBDD_COMPLEX_BIG(MkWalshInterleavedTop(i));
		}

		CFLOBDD_COMPLEX_BIG MkPauliYMatrixInterleaved(unsigned int i)
		{
			assert(i <= CFLOBDD_COMPLEX_BIG::maxLevel);
			return CFLOBDD_COMPLEX_BIG(MkPauliYMatrixInterleavedTop(i));
		}

		CFLOBDD_COMPLEX_BIG MkPauliZMatrixInterleaved(unsigned int i)
		{
			assert(i <= CFLOBDD_COMPLEX_BIG::maxLevel);
			return CFLOBDD_COMPLEX_BIG(MkPauliZMatrixInterleavedTop(i));
		}

		CFLOBDD_COMPLEX_BIG MkSGateInterleaved(unsigned int i)
		{
			assert(i <= CFLOBDD_COMPLEX_BIG::maxLevel);
			return CFLOBDD_COMPLEX_BIG(MkSGateInterleavedTop(i));
		}

		CFLOBDD_COMPLEX_BIG MkPhaseShiftGateInterleaved(unsigned int i, double theta)
		{
			assert(i == 1);
			return CFLOBDD_COMPLEX_BIG(MkPhaseShiftGateInterleavedTop(i, theta));
		}

		CFLOBDD_COMPLEX_BIG MkCNOTInterleaved(unsigned int i)
		{
			return CFLOBDD_COMPLEX_BIG(Matrix1234ComplexFloatBoost::MkCNOTInterleavedTop(i));
		}

		// Create representation of the Inverse Reed-Muller matrix IRM(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields IRM[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDD_COMPLEX_BIG MkInverseReedMullerInterleaved(unsigned int i)
		{
			assert(i <= CFLOBDD_COMPLEX_BIG::maxLevel);
			return CFLOBDD_COMPLEX_BIG(MkInverseReedMullerInterleavedTop(i));
		}

		// ****************************************************************************
		// Matrix-related operations (on matrices with room for two extra vocabularies)
		// ****************************************************************************

		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for
		// two extra vocabularies
		CFLOBDD_COMPLEX_BIG MkWalshVoc13(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD_COMPLEX_BIG::maxLevel);
			return CFLOBDD_COMPLEX_BIG(MkWalshVoc13Top(i));
		}

		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for
		// two extra vocabularies
		CFLOBDD_COMPLEX_BIG MkWalshVoc12(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD_COMPLEX_BIG::maxLevel);
			return CFLOBDD_COMPLEX_BIG(MkWalshVoc12Top(i));
		}


		// Vocabulary shift in a matrix
		CFLOBDD_COMPLEX_BIG MatrixShiftVocs13To24(CFLOBDD_COMPLEX_BIG c)
		{
			return CFLOBDD_COMPLEX_BIG(MatrixShiftVocs13To24Top(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_COMPLEX_BIG MatrixShiftVocs12To34(CFLOBDD_COMPLEX_BIG c)
		{
			return CFLOBDD_COMPLEX_BIG(MatrixShiftVocs12To34Top(c.root));
		}

		CFLOBDD_COMPLEX_BIG PromoteInterleavedTo12(CFLOBDD_COMPLEX_BIG c)
		{
			return CFLOBDD_COMPLEX_BIG(PromoteInterleavedTo12Top(c.root));
		}

		extern CFLOBDD_COMPLEX_BIG Demote12ToInterleaved(CFLOBDD_COMPLEX_BIG c)
		{
			return CFLOBDD_COMPLEX_BIG(Demote12ToInterleavedTop(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_COMPLEX_BIG MatrixShiftVoc43(CFLOBDD_COMPLEX_BIG c)
		{
			return CFLOBDD_COMPLEX_BIG(MatrixShiftVoc43Top(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_COMPLEX_BIG MatrixShiftVoc42(CFLOBDD_COMPLEX_BIG c)
		{
			return CFLOBDD_COMPLEX_BIG(MatrixShiftVoc42Top(c.root));
		}

		// Return the Kronecker product of two matrices
		CFLOBDD_COMPLEX_BIG KroneckerProduct(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_COMPLEX_BIG m2_12To34 = MatrixShiftVocs12To34(m2);

			//std::cout << "---Kronecker Product------------------------------------------------------------" << std::endl;
			//std::cout << "---m1 ------------------------------------------------------------" << std::endl;
			//MatrixPrintRowMajor(m1, std::cout);
			//std::cout << "---m2 ------------------------------------------------------------" << std::endl;
			//MatrixPrintRowMajor(m2, std::cout);
			//std::cout << "---m2_12To34 ------------------------------------------------------------" << std::endl;
			//MatrixPrintRowMajor(m2_12To34, std::cout);
			//CFLOBDD_COMPLEX_BIG temp = m1 * m2_12To34;
			//std::cout << "---temp ------------------------------------------------------------" << std::endl;
			//MatrixPrintRowMajor(temp, std::cout);
			//std::cout << "----------------------------------------------------------------" << std::endl;


			return m1 * m2_12To34;
		}

		//
		// MatrixDetensor
		//
		// Return the result of detensoring a 4-vocabulary matrix
		//
		// CFLOBDD_COMPLEX_BIG MatrixDetensor(CFLOBDD_COMPLEX_BIG k)
		// {
		// 	assert(k.root->level >= 2);

		// 	CFLOBDD_COMPLEX_BIG e = MkDetensorConstraintInterleaved(k.root->level);
		// 	CFLOBDD_COMPLEX_BIG m = k * e;
		// 	CFLOBDD_COMPLEX_BIG p = MatrixProjectVoc23(m);
		// 	CFLOBDD_COMPLEX_BIG ans = MatrixShiftVoc42(p);
		// 	return ans;
		// }

		//
		// MatrixMultiply
		//
		// Return the matrix product of m1 and m2
		//
		// CFLOBDD_COMPLEX_BIG MatrixMultiply(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2)
		// {
		// 	assert(m1.root->level == m2.root->level);
		// 	assert(m1.root->level >= 2);

		// 	CFLOBDD_COMPLEX_BIG k = KroneckerProduct(m1, m2);
		// 	CFLOBDD_COMPLEX_BIG ans = MatrixDetensor(k);
		// 	return ans;
		// }

		CFLOBDD_COMPLEX_BIG MatrixMultiplyV4(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2)
		{
			assert(m1.root->level == m2.root->level);
			assert(m1.root->level >= 2);
			return CFLOBDD_COMPLEX_BIG(MatrixMultiplyV4TopNode(m1.root, m2.root));
		}

		CFLOBDD_COMPLEX_BIG MatrixMultiplyV4WithInfo(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2)
		{
			assert(m1.root->level == m2.root->level);
			assert(m1.root->level >= 2);
			return CFLOBDD_COMPLEX_BIG(MatrixMultiplyV4WithInfoTopNode(m1.root, m2.root));
		}

		// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal.
		// (W,X,Y,Z) s.t. X==Y with interleaved variables
		CFLOBDD_COMPLEX_BIG MkDetensorConstraintInterleaved(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD_COMPLEX_BIG::maxLevel);
			return CFLOBDD_COMPLEX_BIG(MkDetensorConstraintInterleavedTop(i));
		}

		// Create representation of a matrix in which vocabularies 1 and $ are constrained to be equal.
		// (W,X,Y,Z) s.t. W==Z with interleaved variables
		CFLOBDD_COMPLEX_BIG MkCFLOBDDMatrixEqVoc14(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD_COMPLEX_BIG::maxLevel);
			return CFLOBDD_COMPLEX_BIG(MkCFLOBDDMatrixEqVoc14Top(i));
		}

		// Vocabulary projection
		// CFLOBDD_COMPLEX_BIG MatrixProjectVoc23(CFLOBDD_COMPLEX_BIG c)
		// {
		// 	return CFLOBDD_COMPLEX_BIG(MatrixProjectVoc23Top(c.root));
		// }


		// Create representation of the diagonal matrix that is the second component of the Kronecker-product-based
		// construction of the Fourier matrix
		CFLOBDD_COMPLEX_BIG MkFourierDiagonalComponent(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD_COMPLEX_BIG::maxLevel);
			return CFLOBDD_COMPLEX_BIG(MkFourierDiagonalComponentTop(i));
		}

		// Create representation of the matrix for the Discrete Fourier Transform (DFT)
		// The answer matrix is a level-i object that uses interleaved vocabularies
		CFLOBDD_COMPLEX_BIG MkFourierMatrixInterleaved(unsigned int i)
		{
			assert(1 <= i && i < (CFLOBDD_COMPLEX_BIG::maxLevel - 1));  // Need two extra levels: +1 for Kronecker Product; +1 for matrix multiplication

			if (i == 1) {
				return MkWalshInterleaved(1);                                                          // Level 1, interleaved
			}
			else {  // Need to promote all matrices to vocs 1,2 to allow final matrix multiplication
				CFLOBDD_COMPLEX_BIG Id = PromoteInterleavedTo12(MkIdRelationInterleaved(i - 1));                   // Level i, Vocs 1,2

				//std::cout << "---Id------------------------------------------------------------" << std::endl;
				//std::cout << Id << std::endl;
				//Id.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG Id_Demoted = Demote12ToInterleaved(Id);
				/*std::cout << "---Id_Demoted------------------------------------------------------------" << std::endl;
				Matrix1234ComplexDouble::MatrixPrintRowMajor(Id_Demoted, std::cout);
				std::cout << Id_Demoted << std::endl;
				std::cout << "----------------------------------------------------------------" << std::endl;*/

				CFLOBDD_COMPLEX_BIG FourierIMinusOne = PromoteInterleavedTo12(MkFourierMatrixInterleaved(i - 1));  // Level i, Vocs 1,2

				//std::cout << "---FourierIMinusOne------------------------------------------------------------" << std::endl;
				//std::cout << FourierIMinusOne << std::endl;
				//FourierIMinusOne.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG FourierIMinusOne_Demoted = Demote12ToInterleaved(FourierIMinusOne);
				/*std::cout << "---FourierIMinusOne_Demoted------------------------------------------------------------" << std::endl;
				Matrix1234ComplexDouble::MatrixPrintRowMajor(FourierIMinusOne_Demoted, std::cout);
				std::cout << FourierIMinusOne_Demoted << std::endl;
				std::cout << "----------------------------------------------------------------" << std::endl;*/

				CFLOBDD_COMPLEX_BIG term1 = PromoteInterleavedTo12(KroneckerProduct(FourierIMinusOne, Id));        // Level i+1, Vocs 1,2

				//std::cout << "---term1------------------------------------------------------------" << std::endl;
				//std::cout << term1 << std::endl;
				//term1.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG term1_Demoted = Demote12ToInterleaved(term1);
				/*std::cout << "---term1_Demoted------------------------------------------------------------" << std::endl;
				Matrix1234ComplexDouble::MatrixPrintRowMajor(term1_Demoted, std::cout);
				std::cout << term1_Demoted << std::endl;
				std::cout << "----------------------------------------------------------------" << std::endl;*/

				CFLOBDD_COMPLEX_BIG term2 = PromoteInterleavedTo12(MkFourierDiagonalComponent(i));                 // Level i+1, Vocs 1,2

				//std::cout << "---term2------------------------------------------------------------" << std::endl;
				//std::cout << term2 << std::endl;
				//term2.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG term2_Demoted = Demote12ToInterleaved(term2);
				/*std::cout << "---term2_Demoted------------------------------------------------------------" << std::endl;
				Matrix1234ComplexDouble::MatrixPrintRowMajor(term2_Demoted, std::cout);
				std::cout << term2_Demoted << std::endl;
				std::cout << "----------------------------------------------------------------" << std::endl;*/

				CFLOBDD_COMPLEX_BIG term3 = PromoteInterleavedTo12(KroneckerProduct(Id, FourierIMinusOne));        // Level i+1, Vocs 1,2

				//std::cout << "---term3------------------------------------------------------------" << std::endl;
				//std::cout << term3 << std::endl;
				//term3.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG term3_Demoted = Demote12ToInterleaved(term3);
				/*std::cout << "---term3_Demoted------------------------------------------------------------" << std::endl;
				Matrix1234ComplexDouble::MatrixPrintRowMajor(term3_Demoted, std::cout);
				std::cout << term3_Demoted << std::endl;
				std::cout << "----------------------------------------------------------------" << std::endl;*/

				CFLOBDD_COMPLEX_BIG term4 = PromoteInterleavedTo12(MkCFLOBDDMatrixEqVoc14(i) * MkDetensorConstraintInterleaved(i));                     // Level i+1, Vocs 1,2

				//std::cout << "---term4------------------------------------------------------------" << std::endl;
				//std::cout << term4 << std::endl;
				//term4.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG term4_Demoted = Demote12ToInterleaved(term4);
				/*std::cout << "---term4_Demoted------------------------------------------------------------" << std::endl;
				Matrix1234ComplexDouble::MatrixPrintRowMajor(term4_Demoted, std::cout);
				std::cout << term4_Demoted << std::endl;
				std::cout << "----------------------------------------------------------------" << std::endl;*/

				CFLOBDD_COMPLEX_BIG temp12 = MatrixMultiplyV4(term1, term2);                                          // Level i+1, Vocs 1,2

				//std::cout << "---temp12------------------------------------------------------------" << std::endl;
				//std::cout << temp12 << std::endl;
				//temp12.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG temp12_Demoted = Demote12ToInterleaved(temp12);
				/*std::cout << "---temp12_Demoted------------------------------------------------------------" << std::endl;
				Matrix1234ComplexDouble::MatrixPrintRowMajor(temp12_Demoted, std::cout);
				std::cout << temp12_Demoted << std::endl;
				std::cout << "----------------------------------------------------------------" << std::endl;*/

				CFLOBDD_COMPLEX_BIG temp34 = MatrixMultiplyV4(term3, term4);                                          // Level i+1, Vocs 1,2

				//std::cout << "---temp34------------------------------------------------------------" << std::endl;
				//std::cout << temp34 << std::endl;
				//temp34.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG temp34_Demoted = Demote12ToInterleaved(temp34);
				/*std::cout << "---temp34_Demoted------------------------------------------------------------" << std::endl;
				Matrix1234ComplexDouble::MatrixPrintRowMajor(temp34_Demoted, std::cout);
				std::cout << temp34_Demoted << std::endl;
				std::cout << "----------------------------------------------------------------" << std::endl;*/

				CFLOBDD_COMPLEX_BIG temp1234 = MatrixMultiplyV4(temp12, temp34);                                          // Level i+1, Vocs 1,2

				//std::cout << "---temp1234------------------------------------------------------------" << std::endl;
				//std::cout << temp1234 << std::endl;
				//temp1234.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG tempInterleaved = Demote12ToInterleaved(temp1234);                                        // Level i, interleaved
				/*std::cout << "---tempInterleaved------------------------------------------------------------" << std::endl;
				Matrix1234ComplexDouble::MatrixPrintRowMajor(tempInterleaved, std::cout);
				std::cout << tempInterleaved << std::endl;
				std::cout << "----------------------------------------------------------------" << std::endl;*/

				return tempInterleaved;
			}
		}

		CFLOBDD_COMPLEX_BIG ConvertToComplex(CFLOBDD_FOURIER c)
		{
			return CFLOBDD_COMPLEX_BIG(ConvertToComplexTop(c.root));
		}


		CFLOBDD_COMPLEX_BIG MkCPGate(unsigned int i, long c1, long c2, double theta)
		{
			return CFLOBDD_COMPLEX_BIG(MkCPGateTop(i, c1, c2, theta));
		}

		CFLOBDD_COMPLEX_BIG MkExchangeInterleaved(unsigned int i)
		{
			return CFLOBDD_COMPLEX_BIG(MkExchangeInterleavedTop(i));
		}

		CFLOBDD_COMPLEX_BIG MkNegationMatrixInterleaved(unsigned int i)
		{
			return CFLOBDD_COMPLEX_BIG(MkNegationMatrixInterleavedTop(i));
		}

		CFLOBDD_COMPLEX_BIG MkCNOT(unsigned int level, unsigned int n, long int controller, long int controlled){
			return CFLOBDD_COMPLEX_BIG(MkCNOTTopNode(level, n, controller, controlled));
		}

		CFLOBDD_COMPLEX_BIG MkCCNOT(unsigned int level, unsigned int n, long int controller1, long int controller2, long int controlled){
			return CFLOBDD_COMPLEX_BIG(MkCCNOTTopNode(level, n, controller1, controller2, controlled));
		}
		
		CFLOBDD_COMPLEX_BIG MkSwapGate(unsigned int i, long c1, long c2)
		{
			return CFLOBDD_COMPLEX_BIG(MkSwapGateTop(i, c1, c2));
		}

		CFLOBDD_COMPLEX_BIG MkiSwapGate(unsigned int i, long c1, long c2)
		{
			return CFLOBDD_COMPLEX_BIG(MkiSwapGateTop(i, c1, c2));
		}

		CFLOBDD_COMPLEX_BIG MkCSwapGate(unsigned int i, long int c1, long int x1, long int x2)
		{
			return CFLOBDD_COMPLEX_BIG(MkCSwapGateTop(i, c1, x1, x2));
		}

		CFLOBDD_COMPLEX_BIG MatrixShiftToAConnection(CFLOBDD_COMPLEX_BIG c)
		{
			return CFLOBDD_COMPLEX_BIG(MatrixShiftToAConnectionTop(c.root));
		}

		CFLOBDD_COMPLEX_BIG MatrixShiftToBConnection(CFLOBDD_COMPLEX_BIG c)
		{
			return CFLOBDD_COMPLEX_BIG(MatrixShiftToBConnectionTop(c.root));
		}

		// Return the Kronecker product of two matrices
		CFLOBDD_COMPLEX_BIG KroneckerProduct2Vocs(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_COMPLEX_BIG m1_A = MatrixShiftToAConnection(m1);
			CFLOBDD_COMPLEX_BIG m2_B = MatrixShiftToBConnection(m2);
			CFLOBDD_COMPLEX_BIG c = m1_A * m2_B;
			return c;
		}

		CFLOBDD_COMPLEX_BIG MkFourierMatrixInterleavedV4(unsigned int i)
		{
			if (i >= 5){
				DisposeOfPairProductCache();
				InitPairProductCache();
			}
			assert(1 <= i && i < (CFLOBDD_COMPLEX_BIG::maxLevel - 1));  // Need two extra levels: +1 for Kronecker Product; +1 for matrix multiplication

			if (i == 1) {
				return MkWalshInterleaved(1);                                                          // Level 1, interleaved
			}
			else {  // Need to promote all matrices to vocs 1,2 to allow final matrix multiplication


				CFLOBDD_COMPLEX_BIG FourierIMinusOne = PromoteInterleavedTo12(MkFourierMatrixInterleavedV4(i - 1));  // Level i, Vocs 1,2

				std::cout << "---FourierIMinusOne------------------------------------------------------------" << std::endl;
				//std::cout << FourierIMinusOne << std::endl;
				//FourierIMinusOne.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				std::cout << "---Level i: " << i << " ------------------------------------------------------------" << std::endl;


				CFLOBDD_COMPLEX_BIG Id = PromoteInterleavedTo12(MkIdRelationInterleaved(i - 1));                   // Level i-1, Vocs 1

				std::cout << "---Id------------------------------------------------------------" << std::endl;
				//std::cout << Id << std::endl;
				//Id.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;


				CFLOBDD_COMPLEX_BIG term1 = KroneckerProduct(FourierIMinusOne, Id);        // Level i+1, Vocs 1,2

				std::cout << "---term1------------------------------------------------------------" << std::endl;
				//std::cout << term1 << std::endl;
				//term1.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG term2 = MkFourierDiagonalComponent(i);                 // Level i+1, Vocs 1,2

				std::cout << "---term2------------------------------------------------------------" << std::endl;
				//std::cout << term2 << std::endl;
				//term2.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG term3 = KroneckerProduct(Id, FourierIMinusOne);        // Level i+1, Vocs 1,2

				std::cout << "---term3------------------------------------------------------------" << std::endl;
				//std::cout << term3 << std::endl;
				//term3.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG term4 = MkCFLOBDDMatrixEqVoc14(i) * MkDetensorConstraintInterleaved(i);                     // Level i+1, Vocs 1,2

				std::cout << "---term4------------------------------------------------------------" << std::endl;
				//std::cout << term4 << std::endl;
				//term4.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;
				//printMemory();

				CFLOBDD_COMPLEX_BIG temp12 = MatrixMultiplyV4(term1, term2);                                          // Level i+1, Vocs 1,2

				std::cout << "---temp12------------------------------------------------------------" << std::endl;
				//std::cout << temp12 << std::endl;
				//temp12.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG temp34 = MatrixMultiplyV4(term3, term4);                                          // Level i+1, Vocs 1,2
				
				std::cout << "---temp34------------------------------------------------------------" << std::endl;
				//std::cout << temp34 << std::endl;
				//temp34.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;


				CFLOBDD_COMPLEX_BIG temp1234 = MatrixMultiplyV4(temp12, temp34);                                          // Level i+1, Vocs 1,2

				std::cout << "---temp1234------------------------------------------------------------" << std::endl;
				//std::cout << temp1234 << std::endl;
				//temp1234.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				return temp1234;
			}
		}

		CFLOBDD_COMPLEX_BIG MkFourierMatrixInterleavedV4WithInfo(unsigned int i)
		{
			if (i >= 5){
				DisposeOfPairProductCache();
				InitPairProductCache();
			}
			assert(1 <= i && i < (CFLOBDD_COMPLEX_BIG::maxLevel - 1));  // Need two extra levels: +1 for Kronecker Product; +1 for matrix multiplication

			if (i == 1) {
				return MkWalshInterleaved(1);                                                          // Level 1, interleaved
			}
			else {  // Need to promote all matrices to vocs 1,2 to allow final matrix multiplication


				CFLOBDD_COMPLEX_BIG FourierIMinusOne = PromoteInterleavedTo12(MkFourierMatrixInterleavedV4(i - 1));  // Level i, Vocs 1,2

				std::cout << "---FourierIMinusOne------------------------------------------------------------" << std::endl;
				//std::cout << FourierIMinusOne << std::endl;
				//FourierIMinusOne.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				std::cout << "---Level i: " << i << " ------------------------------------------------------------" << std::endl;


				CFLOBDD_COMPLEX_BIG Id = PromoteInterleavedTo12(MkIdRelationInterleaved(i - 1));                   // Level i-1, Vocs 1

				std::cout << "---Id------------------------------------------------------------" << std::endl;
				//std::cout << Id << std::endl;
				//Id.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;


				CFLOBDD_COMPLEX_BIG term1 = KroneckerProduct(FourierIMinusOne, Id);        // Level i+1, Vocs 1,2

				std::cout << "---term1------------------------------------------------------------" << std::endl;
				//std::cout << term1 << std::endl;
				//term1.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG term2 = MkFourierDiagonalComponent(i);                 // Level i+1, Vocs 1,2

				std::cout << "---term2------------------------------------------------------------" << std::endl;
				//std::cout << term2 << std::endl;
				//term2.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG term3 = KroneckerProduct(Id, FourierIMinusOne);        // Level i+1, Vocs 1,2

				std::cout << "---term3------------------------------------------------------------" << std::endl;
				//std::cout << term3 << std::endl;
				//term3.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG term4 = MkCFLOBDDMatrixEqVoc14(i) * MkDetensorConstraintInterleaved(i);                     // Level i+1, Vocs 1,2

				std::cout << "---term4------------------------------------------------------------" << std::endl;
				//std::cout << term4 << std::endl;
				//term4.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;
				//printMemory();

				CFLOBDD_COMPLEX_BIG temp12 = MatrixMultiplyV4WithInfo(term1, term2);                                          // Level i+1, Vocs 1,2

				std::cout << "---temp12------------------------------------------------------------" << std::endl;
				//std::cout << temp12 << std::endl;
				//temp12.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_COMPLEX_BIG temp34 = MatrixMultiplyV4WithInfo(term3, term4);                                          // Level i+1, Vocs 1,2

				std::cout << "---temp34------------------------------------------------------------" << std::endl;
				//std::cout << temp34 << std::endl;
				//temp34.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;



				CFLOBDD_COMPLEX_BIG temp1234 = MatrixMultiplyV4WithInfo(temp12, temp34);                                          // Level i+1, Vocs 1,2

				std::cout << "---temp1234------------------------------------------------------------" << std::endl;
				//std::cout << temp1234 << std::endl;
				//temp1234.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				return temp1234;
			}
		}

	}
}

