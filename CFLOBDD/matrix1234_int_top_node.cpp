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

namespace CFL_OBDD {

	namespace Matrix1234Int {

		void Matrix1234InitializerTop()
		{
			Matrix1234InitializerNode();
			return;
		}

		// Create representation of identity relation (with interleaved variable order).
		// That is, input (x0,y0,x1,y1,...,xN,yN) yield Id[(x0,x1,...,xN)][(y0,y1,...,yN)]
		// which equals 1 iff xi == yi, for 0 <= i <= N.
		CFLOBDDTopNodeIntRefPtr MkIdRelationInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeIntRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			CFLOBDDReturnMapHandle m10;

			tempHandle = MkIdRelationInterleavedNode(i);
			m10.AddToEnd(1);
			m10.AddToEnd(0);
			m10.Canonicalize();
			v = new CFLOBDDTopNode(tempHandle, m10);
			return v;
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDDTopNodeIntRefPtr MkWalshInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeIntRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			CFLOBDDReturnMapHandle m;

			assert(i <= CFLOBDDTopNode::maxLevel);

			tempHandle = MkWalshInterleavedNode(i);
			m.AddToEnd(1);
			m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNode(tempHandle, m);
			return v;
		}

		// Create representation of the Inverse Reed-Muller matrix IRM(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields IRM[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDDTopNodeIntRefPtr MkInverseReedMullerInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeIntRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			CFLOBDDReturnMapHandle m;

			assert(i <= CFLOBDDTopNode::maxLevel);

			tempHandle = MkInverseReedMullerInterleavedNode(i);
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNode(tempHandle, m);
			return v;
		}

		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for two
		// additional vocabularies
		CFLOBDDTopNodeIntRefPtr MkWalshVoc13Top(unsigned int i)
		{
			CFLOBDDTopNodeIntRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			CFLOBDDReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNode::maxLevel);

			tempHandle = MkWalshVoc13Node(i);
			m.AddToEnd(1);
			m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNode(tempHandle, m);
			return v;
		}

		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for two
		// additional vocabularies
		CFLOBDDTopNodeIntRefPtr MkWalshVoc12Top(unsigned int i)
		{
			CFLOBDDTopNodeIntRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			CFLOBDDReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNode::maxLevel);

			tempHandle = MkWalshVoc12Node(i);
			m.AddToEnd(1);
			m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNode(tempHandle, m);
			return v;
		}


		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeIntRefPtr MatrixShiftVocs13To24Top(CFLOBDDTopNodeIntRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVocs13To24Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeIntRefPtr MatrixShiftVocs12To34Top(CFLOBDDTopNodeIntRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVocs12To34Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		/*CFLOBDDTopNodeIntRefPtr MatrixConvertVocs1234To12Top(CFLOBDDTopNodeIntRefPtr n)
		{
			assert(n->level >= 1);
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;
			CFLOBDDNodeHandle tempHandle = MatrixConvertVocs1234To12Node(memoTable, n->rootConnection.entryPointHandle);
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}*/

		CFLOBDDTopNodeIntRefPtr PromoteInterleavedTo12Top(CFLOBDDTopNodeIntRefPtr n)
		{
			assert(1 <= n->rootConnection.entryPointHandle->handleContents->level);
			assert(n->rootConnection.entryPointHandle->handleContents->level < CFLOBDDTopNode::maxLevel);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = PromoteInterleavedTo12Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeIntRefPtr Demote12ToInterleavedTop(CFLOBDDTopNodeIntRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = Demote12ToInterleavedNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}


		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeIntRefPtr MatrixShiftVoc43Top(CFLOBDDTopNodeIntRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVoc43Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeIntRefPtr MatrixShiftVoc42Top(CFLOBDDTopNodeIntRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVoc42Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal:
		// (W,X,Y,Z) s.t. X==Y with interleaved variables
		CFLOBDDTopNodeIntRefPtr MkDetensorConstraintInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeIntRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			CFLOBDDReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNode::maxLevel);

			tempHandle = MkDetensorConstraintInterleavedNode(i);
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			v = new CFLOBDDTopNode(tempHandle, m);
			return v;
		}

		// Create representation of a matrix in which vocabularies 1 and 4 are constrained to be equal:
		// (W,X,Y,Z) s.t. W==Z with interleaved variables
		CFLOBDDTopNodeIntRefPtr MkCFLOBDDMatrixEqVoc14Top(unsigned int i)
		{
			CFLOBDDTopNodeIntRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			CFLOBDDReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNode::maxLevel);

			tempHandle = MkCFLOBDDMatrixEqVoc14Node(i);
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			v = new CFLOBDDTopNode(tempHandle, m);
			return v;
		}

		// Vocabulary projection
		// CFLOBDDTopNodeIntRefPtr MatrixProjectVoc23Top(CFLOBDDTopNodeIntRefPtr n)
		// {
		// 	assert(n->rootConnection.entryPointHandle.handleContents->level >= 2);

		// 	CFLOBDDLinearMapMemoTableRefPtr memoTable = new CFLOBDDLinearMapMemoTable;

		// 	CFLOBDDTopNodeLinearMapRefPtr temp = MatrixProjectVoc23Node(memoTable, n->rootConnection.entryPointHandle, TopLevelVisit);

		// 	// Begin construction of the appropriate CFLOBDDReturnMapHandle by applying each LinearMap in temp's
		// 	// CFLOBDDLinearMapTupleHandle to n's CFLOBDDReturnMapHandle
		// 	CFLOBDDReturnMapHandle rmh;
		// 	for (unsigned int i = 0; i < temp->rootConnection.returnMapHandle.Size(); i++) {
		// 		rmh.AddToEnd(ApplyLinearMap(temp->rootConnection.returnMapHandle[i], n->rootConnection.returnMapHandle));
		// 	}

		// 	CFLOBDDNodeHandle tempHandle = temp->rootConnection.entryPointHandle;

		// 	// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
		// 	ReductionMapHandle inducedReductionMapHandle;
		// 	CFLOBDDReturnMapHandle inducedReturnMap;
		// 	rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
		// 	//CFLOBDDNodeHandle::InitReduceCache();
		// 	CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
		// 	//CFLOBDDNodeHandle::DisposeOfReduceCache();

		// 	// Create and return CFLOBDDTopNode
		// 	return(new CFLOBDDTopNode(reduced_tempHandle, inducedReturnMap));
		// }

		CFLOBDDTopNodeIntRefPtr MkFourierDiagonalComponentTop(unsigned int i)
		{
			CFLOBDDTopNodeIntRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			CFLOBDDReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNode::maxLevel);

			tempHandle = MkFourierDiagonalComponentNode(i);
			unsigned long int n = 1UL << (1 << (i - 2));
			unsigned long int nSq = n * n;
			assert(nSq != 0);   // Check for overflow
			for (unsigned long int j = 0; j < n; j++) {
				for (unsigned long int k = 0; k < n; k++) {
					m.AddToEnd((j*k) % n + 1);    // Eventually w^(j*k), where w is the (canonical?) n-th root of unity
					if (j == 0 && k == 0)
						m.AddToEnd(-1000);     // Eventually m.AddToEnd(0)
				}
			}
			m.Canonicalize();

			// Perform reduction on tempHandle, with respect to the common elements that m maps together
			ReductionMapHandle inducedReductionMapHandle;
			CFLOBDDReturnMapHandle inducedReturnMap;
			m.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//CFLOBDDNodeHandle::InitReduceCache();
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			//CFLOBDDNodeHandle::DisposeOfReduceCache();

			// Create and return CFLOBDDTopNode
			return(new CFLOBDDTopNode(reduced_tempHandle, inducedReturnMap));
		}

		CFLOBDDTopNodeIntRefPtr MatrixTransposeTop(CFLOBDDTopNodeIntRefPtr n)
		{
			std::unordered_map<CFLOBDDNodeHandle, std::pair<CFLOBDDNodeHandle, CFLOBDDReturnMapHandle>, 
				CFLOBDDNodeHandle::CFLOBDDNodeHandle_Hash> hashMap;
			auto pc = MatrixTransposeNode(hashMap, *(n->rootConnection.entryPointHandle));
			CFLOBDDNodeHandle temp = pc.first;
			ReductionMapHandle reductionMapHandle;
			CFLOBDDReturnMapHandle v;
			for (unsigned int i = 0; i < pc.second.Size(); i++){
				reductionMapHandle.AddToEnd(i);
				v.AddToEnd(n->rootConnection.returnMapHandle[pc.second[i]]);
			}
			reductionMapHandle.Canonicalize();
			v.Canonicalize();
			temp = temp.Reduce(reductionMapHandle, v.Size(), true);
			return new CFLOBDDTopNode(temp, v);
		}

		//
		// Print the matrix in row-major order
		// For each entry, the assignment created has Booleans for the bits in root-to-leaf order.
		// The loop that fills the assignment runs least-significant bit to most-significant bit,
		// interleaving row and column bits.  The most-significant bit (of the row value) corresponds
		// to the root; the least-significant bit (of the col value) corresponds to the leaf.
		//
		void MatrixPrintRowMajorTop(CFLOBDDTopNodeIntRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 2 && level <= 4) {
				unsigned int indexBits = 1 << (level - 1);
				unsigned int totalBits = 2 * indexBits;
				unsigned long int rows = 1UL << indexBits;
				unsigned long int cols = rows;
				SH_OBDD::Assignment a(totalBits);
				for (unsigned long int i = 0UL; i < rows; i++) { // Vocs 1 and 3
					// Fill in the even positions of a 
					unsigned long int maskVoc1 = 1UL << (indexBits / 2);
					unsigned long int maskVoc3 = 1UL;
					for (int k = (indexBits / 2 - 1); k >= 0; k--) {
						a[4 * k] = (i & maskVoc1);
						a[4 * k + 2] = (i & maskVoc3);
						maskVoc1 = maskVoc1 << 1;
						maskVoc3 = maskVoc3 << 1;
					}
					for (unsigned long int j = 0UL; j < cols; j++) {  // Vocs 2 and 4
						// Fill in the odd positions of a
						unsigned long int maskVoc2 = 1UL << (indexBits / 2);
						unsigned long int maskVoc4 = 1UL;
						for (int k = (indexBits / 2 - 1); k >= 0; k--) {
							a[4 * k + 1] = (j & maskVoc2);
							a[4 * k + 3] = (j & maskVoc4);
							maskVoc2 = maskVoc2 << 1;
							maskVoc4 = maskVoc4 << 1;
						}
						//out << "i = " << i << ", j = " << j << ", a = ";
						//a.print(out);
						//out << std::endl;
						int b;
						b = n->EvaluateIteratively(a);
						out << b << " ";
					}
					out << std::endl;
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [2 .. 4]" << std::endl;
			}
		}

	}
}

