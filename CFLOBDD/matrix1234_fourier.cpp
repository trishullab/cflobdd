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
#include "matrix1234_fourier_top_node.h"
#include "matrix1234_fourier.h"

namespace CFL_OBDD {

	namespace Matrix1234Fourier {

		void Matrix1234Initializer()
		{
			Matrix1234InitializerTop();
			return;
		}

		void MatrixPrintRowMajorInterleaved(CFLOBDD_FOURIER c, std::ostream & out)
		{
			MatrixPrintRowMajorInterleavedTop(c.root, out);
			return;
		}

		void MatrixPrintRowMajor(CFLOBDD_FOURIER c, std::ostream & out)
		{
			MatrixPrintRowMajorTop(c.root, out);
			return;
		}

		// Create representation of identity relation
		CFLOBDD_FOURIER MkIdRelationInterleaved(unsigned int i)
		{
			return CFLOBDD_FOURIER(MkIdRelationInterleavedTop(i));
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDD_FOURIER MkWalshInterleaved(unsigned int i)
		{
			assert(i <= CFLOBDD_FOURIER::maxLevel);
			return CFLOBDD_FOURIER(MkWalshInterleavedTop(i));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_FOURIER MatrixShiftVocs13To24(CFLOBDD_FOURIER c)
		{
			return CFLOBDD_FOURIER(MatrixShiftVocs13To24Top(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_FOURIER MatrixShiftVocs12To34(CFLOBDD_FOURIER c)
		{
			return CFLOBDD_FOURIER(MatrixShiftVocs12To34Top(c.root));
		}

		CFLOBDD_FOURIER PromoteInterleavedTo12(CFLOBDD_FOURIER c)
		{
			return CFLOBDD_FOURIER(PromoteInterleavedTo12Top(c.root));
		}

		extern CFLOBDD_FOURIER Demote12ToInterleaved(CFLOBDD_FOURIER c)
		{
			return CFLOBDD_FOURIER(Demote12ToInterleavedTop(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_FOURIER MatrixShiftVoc43(CFLOBDD_FOURIER c)
		{
			return CFLOBDD_FOURIER(MatrixShiftVoc43Top(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_FOURIER MatrixShiftVoc42(CFLOBDD_FOURIER c)
		{
			return CFLOBDD_FOURIER(MatrixShiftVoc42Top(c.root));
		}

		// Return the Kronecker product of two matrices
		CFLOBDD_FOURIER KroneckerProduct(CFLOBDD_FOURIER m1, CFLOBDD_FOURIER m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_FOURIER m2_12To34 = MatrixShiftVocs12To34(m2);

			//std::cout << "---Kronecker Product------------------------------------------------------------" << std::endl;
			//std::cout << "---m1 ------------------------------------------------------------" << std::endl;
			//MatrixPrintRowMajor(m1, std::cout);
			//std::cout << "---m2 ------------------------------------------------------------" << std::endl;
			//MatrixPrintRowMajor(m2, std::cout);
			//std::cout << "---m2_12To34 ------------------------------------------------------------" << std::endl;
			//MatrixPrintRowMajor(m2_12To34, std::cout);
			//CFLOBDD_FOURIER temp = m1 * m2_12To34;
			//std::cout << "---temp ------------------------------------------------------------" << std::endl;
			//MatrixPrintRowMajor(temp, std::cout);
			//std::cout << "----------------------------------------------------------------" << std::endl;


			return m1 * m2_12To34;
		}


		CFLOBDD_FOURIER MatrixMultiplyV4(CFLOBDD_FOURIER m1, CFLOBDD_FOURIER m2)
		{
			assert(m1.root->level == m2.root->level);
			assert(m1.root->level >= 2);
			return CFLOBDD_FOURIER(MatrixMultiplyV4TopNode(m1.root, m2.root));
		}

		CFLOBDD_FOURIER MatrixMultiplyV4WithInfo(CFLOBDD_FOURIER m1, CFLOBDD_FOURIER m2)
		{
			assert(m1.root->level == m2.root->level);
			assert(m1.root->level >= 2);
			return CFLOBDD_FOURIER(MatrixMultiplyV4WithInfoTopNode(m1.root, m2.root));
		}

		// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal.
		// (W,X,Y,Z) s.t. X==Y with interleaved variables
		CFLOBDD_FOURIER MkDetensorConstraintInterleaved(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD_FOURIER::maxLevel);
			return CFLOBDD_FOURIER(MkDetensorConstraintInterleavedTop(i));
		}

		// Create representation of a matrix in which vocabularies 1 and $ are constrained to be equal.
		// (W,X,Y,Z) s.t. W==Z with interleaved variables
		CFLOBDD_FOURIER MkCFLOBDDMatrixEqVoc14(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD_FOURIER::maxLevel);
			return CFLOBDD_FOURIER(MkCFLOBDDMatrixEqVoc14Top(i));
		}

		CFLOBDD_FOURIER MatrixShiftToAConnection(CFLOBDD_FOURIER c)
		{
			return CFLOBDD_FOURIER(MatrixShiftToAConnectionTop(c.root));
		}

		CFLOBDD_FOURIER MatrixShiftToBConnection(CFLOBDD_FOURIER c)
		{
			return CFLOBDD_FOURIER(MatrixShiftToBConnectionTop(c.root));
		}


		// Create representation of the diagonal matrix that is the second component of the Kronecker-product-based
		// construction of the Fourier matrix
		CFLOBDD_FOURIER MkFourierDiagonalComponent(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD_FOURIER::maxLevel);
			return CFLOBDD_FOURIER(MkFourierDiagonalComponentTop(i));
		}

		CFLOBDD_FOURIER MkFourierDiagonalComponent2Vocs(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD_FOURIER::maxLevel);
			return CFLOBDD_FOURIER(MkFourierDiagonalComponent2VocsTop(i));
		}

		// Return the Kronecker product of two matrices
		CFLOBDD_FOURIER KroneckerProduct2Vocs(CFLOBDD_FOURIER m1, CFLOBDD_FOURIER m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_FOURIER m1_A = MatrixShiftToAConnection(m1);
			CFLOBDD_FOURIER m2_B = MatrixShiftToBConnection(m2);
			CFLOBDD_FOURIER c = m1_A * m2_B;
			return c;
		}

		CFLOBDD_FOURIER MkFourierMatrixInterleavedV4(unsigned int i)
		{
			if (i >= 5){
				DisposeOfPairProductCache();
				InitPairProductCache();
			}
			assert(1 <= i && i < (CFLOBDD_FOURIER::maxLevel - 1));  // Need two extra levels: +1 for Kronecker Product; +1 for matrix multiplication

			if (i == 1) {
				return MkWalshInterleaved(1);                                                          // Level 1, interleaved
			}
			else {  // Need to promote all matrices to vocs 1,2 to allow final matrix multiplication


				CFLOBDD_FOURIER FourierIMinusOne = PromoteInterleavedTo12(MkFourierMatrixInterleavedV4(i - 1));  // Level i, Vocs 1,2

				std::cout << "---FourierIMinusOne------------------------------------------------------------" << std::endl;
				//std::cout << FourierIMinusOne << std::endl;
				//FourierIMinusOne.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				std::cout << "---Level i: " << i << " ------------------------------------------------------------" << std::endl;


				CFLOBDD_FOURIER Id = PromoteInterleavedTo12(MkIdRelationInterleaved(i - 1));                   // Level i-1, Vocs 1

				std::cout << "---Id------------------------------------------------------------" << std::endl;
				//std::cout << Id << std::endl;
				//Id.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;


				CFLOBDD_FOURIER term1 = KroneckerProduct(FourierIMinusOne, Id);        // Level i+1, Vocs 1,2

				std::cout << "---term1------------------------------------------------------------" << std::endl;
				//std::cout << term1 << std::endl;
				//term1.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_FOURIER term2 = MkFourierDiagonalComponent(i);                 // Level i+1, Vocs 1,2

				std::cout << "---term2------------------------------------------------------------" << std::endl;
				//std::cout << term2 << std::endl;
				//term2.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_FOURIER term3 = KroneckerProduct(Id, FourierIMinusOne);        // Level i+1, Vocs 1,2

				std::cout << "---term3------------------------------------------------------------" << std::endl;
				//std::cout << term3 << std::endl;
				//term3.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_FOURIER term4 = MkCFLOBDDMatrixEqVoc14(i) * MkDetensorConstraintInterleaved(i);                     // Level i+1, Vocs 1,2

				std::cout << "---term4------------------------------------------------------------" << std::endl;
				//std::cout << term4 << std::endl;
				//term4.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;
				//printMemory();

				CFLOBDD_FOURIER temp12 = MatrixMultiplyV4(term1, term2);                                          // Level i+1, Vocs 1,2

				//Matrix1234Fourier::MatrixPrintRowMajorInterleaved(term1, std::cout);
				//Matrix1234Fourier::MatrixPrintRowMajorInterleaved(term2, std::cout);
				std::cout << "---temp12------------------------------------------------------------" << std::endl;
				//std::cout << temp12 << std::endl;
				//temp12.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_FOURIER temp34 = MatrixMultiplyV4(term3, term4);                                          // Level i+1, Vocs 1,2
				
				std::cout << "---temp34------------------------------------------------------------" << std::endl;
				//std::cout << temp34 << std::endl;
				//temp34.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;


				CFLOBDD_FOURIER temp1234 = MatrixMultiplyV4(temp12, temp34);                                          // Level i+1, Vocs 1,2

				std::cout << "---temp1234------------------------------------------------------------" << std::endl;
				//std::cout << temp1234 << std::endl;
				//temp1234.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				return temp1234;
			}
		}

		CFLOBDD_FOURIER MkFourierMatrixInterleavedV4WithInfo(unsigned int i)
		{
			/*if (i >= 5){
				DisposeOfPairProductCache();
				InitPairProductCache();
			}*/
			assert(1 <= i && i < (CFLOBDD_FOURIER::maxLevel - 1));  // Need two extra levels: +1 for Kronecker Product; +1 for matrix multiplication

			if (i == 1) {
				return MkWalshInterleaved(1);                                                          // Level 1, interleaved
			}
			else {  // Need to promote all matrices to vocs 1,2 to allow final matrix multiplication


				CFLOBDD_FOURIER FourierIMinusOne = PromoteInterleavedTo12(MkFourierMatrixInterleavedV4WithInfo(i - 1));  // Level i, Vocs 1,2

				std::cout << "---FourierIMinusOne------------------------------------------------------------" << std::endl;
				//std::cout << FourierIMinusOne << std::endl;
				//FourierIMinusOne.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				std::cout << "---Level i: " << i << " ------------------------------------------------------------" << std::endl;


				CFLOBDD_FOURIER Id = PromoteInterleavedTo12(MkIdRelationInterleaved(i - 1));                   // Level i-1, Vocs 1

				std::cout << "---Id------------------------------------------------------------" << std::endl;
				//std::cout << Id << std::endl;
				//Id.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;


				CFLOBDD_FOURIER term1 = KroneckerProduct(FourierIMinusOne, Id);        // Level i+1, Vocs 1,2

				std::cout << "---term1------------------------------------------------------------" << std::endl;
				//std::cout << term1 << std::endl;
				//term1.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_FOURIER term2 = MkFourierDiagonalComponent(i);                 // Level i+1, Vocs 1,2
				Matrix1234Fourier::MatrixPrintRowMajorInterleaved(term2, std::cout);

				std::cout << "---term2------------------------------------------------------------" << std::endl;
				//std::cout << term2 << std::endl;
				//term2.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_FOURIER term3 = KroneckerProduct(Id, FourierIMinusOne);        // Level i+1, Vocs 1,2

				std::cout << "---term3------------------------------------------------------------" << std::endl;
				//std::cout << term3 << std::endl;
				//term3.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_FOURIER term4 = MkCFLOBDDMatrixEqVoc14(i) * MkDetensorConstraintInterleaved(i);                     // Level i+1, Vocs 1,2
				//Matrix1234Fourier::MatrixPrintRowMajorInterleaved(term4, std::cout);

				std::cout << "---term4------------------------------------------------------------" << std::endl;
				//std::cout << term4 << std::endl;
				//term4.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;
				//printMemory();

				CFLOBDD_FOURIER temp12 = MatrixMultiplyV4WithInfo(term1, term2);                                          // Level i+1, Vocs 1,2
				
				std::cout << "---temp12------------------------------------------------------------" << std::endl;
				//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
				//temp12.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
				//std::cout << "Step i : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
				//std::cout << temp12 << std::endl;
				//temp12.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				CFLOBDD_FOURIER temp34 = MatrixMultiplyV4WithInfo(term3, term4);                                          // Level i+1, Vocs 1,2

				std::cout << "---temp34------------------------------------------------------------" << std::endl;
				//temp34.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
				//std::cout << "Step i : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
				//std::cout << temp34 << std::endl;
				//temp34.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;



				CFLOBDD_FOURIER temp1234 = MatrixMultiplyV4(temp12, temp34);                                          // Level i+1, Vocs 1,2

				std::cout << "---temp1234------------------------------------------------------------" << std::endl;
				//temp1234.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
				//std::cout << "Step i : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
				//std::cout << temp1234 << std::endl;
				//temp1234.PrintYield(&std::cout);
				//std::cout << "----------------------------------------------------------------" << std::endl;

				return temp1234;
			}
		}

	}
}

