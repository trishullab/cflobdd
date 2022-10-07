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
#include "matrix1234_node.h"
#include "matrix1234_fourier_top_node.h"
#include "assignment.h"

namespace CFL_OBDD {

	namespace Matrix1234Fourier {

		void Matrix1234InitializerTop()
		{
			Matrix1234InitializerNode();
			return;
		}

		// MatrixPrintRowMajorInterleaved
		//
		// Given a CFLOBDD n that represents a 2-vocabulary matrix M_n(W,X) with the
		// vocabulary interleaving (W |x| X), print M_n in row-major order as if M_n were the matrix M(W, X).
		//
		// For each entry, the assignment created has Booleans for the bits in root-to-leaf order.
		// The loop that fills the assignment runs least-significant bit to most-significant bit,
		// interleaving row and column bits.  The most-significant bit (of the row value) corresponds
		// to the root; the least-significant bit (of the col value) corresponds to the leaf.
		//
		void MatrixPrintRowMajorInterleavedTop(CFLOBDDTopNodeFourierRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 1 && level <= 4) {
				unsigned int indexBits = 1 << (level - 1);
				unsigned int totalBits = 2 * indexBits;
				unsigned long int rows = 1UL << indexBits;
				unsigned long int cols = rows;
				SH_OBDD::Assignment a(totalBits);
				// For each row, column position (i,j) in matrix M_n, set a to the assignment
				// such that evaluating CFLOBDD n w.r.t. a yields the value of M_n[i,j]
				for (unsigned long int i = 0UL; i < rows; i++) {  // Voc1
					// Fill in the even positions of a
					//  o the least-significant bits of a (= a[large]) are obtained from the least-signicant bits of i
					//  o the most-significant bits of a (= a[small]) are obtained from the most-significant bits of i
					unsigned long int mask1 = 1UL;
					for (int k = indexBits - 1; k >= 0; k--) {  // least-significant to most-significant bits
						a[2 * k] = (i & mask1);
						mask1 = mask1 << 1;  // Prepare mask1 for the next-most-significant bit
					}
					for (unsigned long int j = 0UL; j < cols; j++) {  // Voc2
						// Fill in the odd positions of a
						//  o the least-significant bits of a (= a[large]) are obtained from the least-signicant bits of j
						//  o the most-significant bits of a (= a[small]) are obtained from the most-significant bits of j
						unsigned long int mask2 = 1UL;
						for (int k = indexBits - 1; k >= 0; k--) {  // least-significant to most-significant bits
							a[2 * k + 1] = (j & mask2);
							mask2 = mask2 << 1;  // Prepare mask2 for the next-most-significant bit
						}
						//out << "i = " << i << ", j = " << j << ", a = ";
						//a.print(out);
						//out << std::endl;
						fourierSemiring b;
						b = n->EvaluateIteratively(a);
						//out << b << " ";
						if (!(b.GetRingSize() == 1 && b.GetVal() == 0))
							out << i << " " << j << " " << b << std::endl;
					}
					//out << std::endl;
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [1 .. 4]" << std::endl;
			}
		}

		// MatrixPrintRowMajor
		//
		// Given a CFLOBDD n that represents a 4-vocabulary matrix M_n(W,X,Y,Z) with the
		// vocabulary interleaving (W |x| X) |x| (Y |x| Z) = (W |x| Y |x| X |x| Z),
		// print M_n in row-major order as if M_n were the matrix M(W || X, Y || Z).
		//
		// For each entry, the assignment created has Booleans for the bits in root-to-leaf order.
		// The loop that fills the assignment runs least-significant bit to most-significant bit,
		// interleaving row and column bits.  The most-significant bit (of the row value) corresponds
		// to the root; the least-significant bit (of the col value) corresponds to the leaf.
		//
		void MatrixPrintRowMajorTop(CFLOBDDTopNodeFourierRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 2 && level <= 4) {
				unsigned int indexBits = 1 << (level - 1);
				unsigned int totalBits = 2 * indexBits;
				unsigned long int rows = 1UL << indexBits;
				unsigned long int cols = rows;
				SH_OBDD::Assignment a(totalBits);
				// For each row, column position (i,j) in matrix M_n, set a to the assignment
				// such that evaluating CFLOBDD n w.r.t. a yields the value of M_n[i,j]
				for (unsigned long int i = 0UL; i < rows; i++) { // Vocs 1 and 3
					// Fill in the even positions of a:
					//  o Voc1 bits are obtained from the top (most-significant) half of i
					//  o Voc3 bits are obtained from the bottom (least-significant) half of i
					// For each vocabulary,
					//  o the least-significant bits of a (= a[large]) are obtained from the least-signicant bits
					//    of the appropriate half of i
					//  o the most-significant bits of a (= a[small]) are obtained from the most-significant bits
					//    of the appropriate half of i
					unsigned long int maskVoc1 = 1UL << (indexBits / 2);
					unsigned long int maskVoc3 = 1UL;
					for (int k = (indexBits / 2 - 1); k >= 0; k--) {  // least-significant to most-significant bits
						a[4 * k] = (i & maskVoc1);
						a[4 * k + 2] = (i & maskVoc3);
						maskVoc1 = maskVoc1 << 1;  // Prepare mask1 for the next-most-significant bit
						maskVoc3 = maskVoc3 << 1;  // Prepare mask3 for the next-most-significant bit
					}
					for (unsigned long int j = 0UL; j < cols; j++) {  // Vocs 2 and 4
						// Fill in the odd positions of a:
						//  o Voc2 bits are obtained from the top (most-significant) half of j
						//  o Voc4 bits are obtained from the bottom (least-significant) half of j
						// For each vocabulary,
						//  o the least-significant bits of a (= a[large]) are obtained from the least-signicant bits
						//    of the appropriate half of j
						//  o the most-significant bits of a (= a[small]) are obtained from the most-significant bits
						//    of the appropriate half of j
						unsigned long int maskVoc2 = 1UL << (indexBits / 2);
						unsigned long int maskVoc4 = 1UL;
						for (int k = (indexBits / 2 - 1); k >= 0; k--) {  // least-significant to most-significant bits
							a[4 * k + 1] = (j & maskVoc2);
							a[4 * k + 3] = (j & maskVoc4);
							maskVoc2 = maskVoc2 << 1;  // Prepare mask2 for the next-most-significant bit
							maskVoc4 = maskVoc4 << 1;  // Prepare mask4 for the next-most-significant bit
						}
						//out << "i = " << i << ", j = " << j << ", a = ";
						//a.print(out);
						//out << std::endl;
						fourierSemiring b;
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

		// Create representation of identity relation (with interleaved variable order).
		// That is, input (x0,y0,x1,y1,...,xN,yN) yield Id[(x0,x1,...,xN)][(y0,y1,...,yN)]
		// which equals 1 iff xi == yi, for 0 <= i <= N.
		CFLOBDDTopNodeFourierRefPtr MkIdRelationInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeFourierRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FourierReturnMapHandle m10;

			tempHandle = MkIdRelationInterleavedNode(i);
			m10.AddToEnd(fourierSemiring(1, 1));
			m10.AddToEnd(fourierSemiring(0, 1));
			m10.Canonicalize();
			v = new CFLOBDDTopNodeFourier(tempHandle, m10);
			return v;
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDDTopNodeFourierRefPtr MkWalshInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeFourierRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FourierReturnMapHandle m;

			assert(i <= CFLOBDDTopNodeFourier::maxLevel);

			tempHandle = MkWalshInterleavedNode(i);
			m.AddToEnd(fourierSemiring(0, 2));
			m.AddToEnd(fourierSemiring(1, 2));
			m.Canonicalize();
			v = new CFLOBDDTopNodeFourier(tempHandle, m);
			return v;
		}

		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeFourierRefPtr MatrixShiftVocs13To24Top(CFLOBDDTopNodeFourierRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVocs13To24Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFourierRefPtr v = new CFLOBDDTopNodeFourier(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeFourierRefPtr MatrixShiftVocs12To34Top(CFLOBDDTopNodeFourierRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVocs12To34Node(memoTable,*(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFourierRefPtr v = new CFLOBDDTopNodeFourier(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeFourierRefPtr PromoteInterleavedTo12Top(CFLOBDDTopNodeFourierRefPtr n)
		{
			assert(1 <= n->rootConnection.entryPointHandle->handleContents->level);
			assert(n->rootConnection.entryPointHandle->handleContents->level < CFLOBDDTopNodeFourier::maxLevel);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = PromoteInterleavedTo12Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFourierRefPtr v = new CFLOBDDTopNodeFourier(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeFourierRefPtr Demote12ToInterleavedTop(CFLOBDDTopNodeFourierRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = Demote12ToInterleavedNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFourierRefPtr v = new CFLOBDDTopNodeFourier(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}


		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeFourierRefPtr MatrixShiftVoc43Top(CFLOBDDTopNodeFourierRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVoc43Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFourierRefPtr v = new CFLOBDDTopNodeFourier(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeFourierRefPtr MatrixShiftVoc42Top(CFLOBDDTopNodeFourierRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVoc42Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFourierRefPtr v = new CFLOBDDTopNodeFourier(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal:
		// (W,X,Y,Z) s.t. X==Y with interleaved variables
		CFLOBDDTopNodeFourierRefPtr MkDetensorConstraintInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeFourierRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FourierReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeFourier::maxLevel);

			tempHandle = MkDetensorConstraintInterleavedNode(i);
			m.AddToEnd(fourierSemiring(1, 1));
			m.AddToEnd(fourierSemiring(0, 1));
			m.Canonicalize();
			v = new CFLOBDDTopNodeFourier(tempHandle, m);
			return v;
		}

		// Create representation of a matrix in which vocabularies 1 and 4 are constrained to be equal:
		// (W,X,Y,Z) s.t. W==Z with interleaved variables
		CFLOBDDTopNodeFourierRefPtr MkCFLOBDDMatrixEqVoc14Top(unsigned int i)
		{
			CFLOBDDTopNodeFourierRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FourierReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeFourier::maxLevel);

			tempHandle = MkCFLOBDDMatrixEqVoc14Node(i);
			m.AddToEnd(fourierSemiring(1, 1));
			m.AddToEnd(fourierSemiring(0, 1));
			m.Canonicalize();
			v = new CFLOBDDTopNodeFourier(tempHandle, m);
			return v;
		}

		CFLOBDDTopNodeFourierRefPtr MkFourierDiagonalComponentTop(unsigned int i)
		{
			CFLOBDDTopNodeFourierRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FourierReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeFourier::maxLevel);

			tempHandle = MkFourierDiagonalComponentNode(i);
			unsigned long int n = 1UL << (1 << (i - 2));
			unsigned long int log_n = 1 << (i - 2);
			unsigned long int nSq = 1UL << (1 << (i - 1));

			// Create a listing of the diagonal entries for a matrix with the interleaved order
			// (Voc1 || Voc4) |x| (Voc2 || Voc3)
			int *diagonalEntries = new (std::nothrow) int[nSq];
			if (diagonalEntries == nullptr) {
				std::cout << "memory allocation failed in MkFourierDiagonalComponentNode" << std::endl;
				abort();
			}
			for (unsigned long int j = 0; j < n; j++) {
				for (unsigned long int k = 0; k < n; k++) {
					diagonalEntries[j*n + k] = (j * k) % nSq;
				}
			}
			// Create ReturnMapHandle m in a way that respects interleaved order
			// Voc1 |x| Voc2 |x| Voc3 |x| Voc4.
			std::unordered_map<int, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
			for (unsigned long int j = 0; j < nSq; j++) {
				// Set adjustedIndex to the unshuffle of the bits of j
				unsigned long int adjustedIndex = 0;
				unsigned long int maskVoc34 = 1UL;
				unsigned long int maskVoc12 = 2UL;
				for (unsigned long int k = 0; k < log_n; k++) {
					if (j & maskVoc34)
						adjustedIndex += (1UL << k);
					if (j & maskVoc12)
						adjustedIndex += (1UL << (log_n + k));
					maskVoc34 = maskVoc34 << 2;
					maskVoc12 = maskVoc12 << 2;
				}
				auto it = reductionMap.find(diagonalEntries[adjustedIndex]);
				if (it == reductionMap.end()){
					m.AddToEnd(fourierSemiring(diagonalEntries[adjustedIndex], nSq));
					reductionMap.emplace(diagonalEntries[adjustedIndex], m.Size() - 1);
					reductionMapHandle.AddToEnd(m.Size() - 1);
				}
				else{
					reductionMapHandle.AddToEnd(reductionMap[diagonalEntries[adjustedIndex]]);
				}
				if (j == 0){
					m.AddToEnd(fourierSemiring(0, 1));
					reductionMap.emplace(-1, m.Size() - 1);
					reductionMapHandle.AddToEnd(m.Size() - 1);
				}
			}
			delete[] diagonalEntries;
			m.Canonicalize();
			reductionMapHandle.Canonicalize();

			// Perform reduction on tempHandle, with respect to the common elements that m maps together
			//ReductionMapHandle inducedReductionMapHandle;
			//FourierReturnMapHandle inducedReturnMap;
			//m.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//     CFLOBDDNodeHandle::InitReduceCache();
			//CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, m.Size(), true);
			//     CFLOBDDNodeHandle::DisposeOfReduceCache();

			// Create and return CFLOBDDTopNodeFourier
			//return(new CFLOBDDTopNodeFourier(reduced_tempHandle, inducedReturnMap));
			return(new CFLOBDDTopNodeFourier(reduced_tempHandle, m));
		}

		CFLOBDDTopNodeFourierRefPtr MkFourierDiagonalComponent2VocsTop(unsigned int i)
		{
			CFLOBDDTopNodeFourierRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FourierReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeFourier::maxLevel);

			tempHandle = MkFourierDiagonalComponentNode(i);
			unsigned long int n = 1UL << (1 << (i - 2));
			unsigned long int log_n = 1 << (i - 2);
			unsigned long int nSq = 1UL << (1 << (i - 1));

			// Create a listing of the diagonal entries for a matrix with the interleaved order
			// (Voc1 || Voc4) |x| (Voc2 || Voc3)
			int *diagonalEntries = new (std::nothrow) int[nSq];
			if (diagonalEntries == nullptr) {
				std::cout << "memory allocation failed in MkFourierDiagonalComponentNode" << std::endl;
				abort();
			}
			for (unsigned long int j = 0; j < n; j++) {
				for (unsigned long int k = 0; k < n; k++) {
					diagonalEntries[j*n + k] = (j * k) % nSq;
				}
			}
			
			std::unordered_map<int, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
			for (unsigned long int j = 0; j < nSq; j++) {
				auto it = reductionMap.find(diagonalEntries[j]);
				if (it == reductionMap.end()){
					m.AddToEnd(fourierSemiring(diagonalEntries[j], nSq));
					reductionMap.emplace(diagonalEntries[j], m.Size() - 1);
					reductionMapHandle.AddToEnd(m.Size() - 1);
				}
				else{
					reductionMapHandle.AddToEnd(it->second);
				}
				if (j == 0){
					m.AddToEnd(fourierSemiring(0, 1));
					reductionMap.emplace(-1, m.Size() - 1);
					reductionMapHandle.AddToEnd(m.Size() - 1);
				}
			}
			delete[] diagonalEntries;
			m.Canonicalize();
			reductionMapHandle.Canonicalize();

			// Perform reduction on tempHandle, with respect to the common elements that m maps together
			//ReductionMapHandle inducedReductionMapHandle;
			//FourierReturnMapHandle inducedReturnMap;
			//m.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//     CFLOBDDNodeHandle::InitReduceCache();
			//CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, m.Size(), true);
			//     CFLOBDDNodeHandle::DisposeOfReduceCache();

			// Create and return CFLOBDDTopNodeFourier
			//return(new CFLOBDDTopNodeFourier(reduced_tempHandle, inducedReturnMap));
			return(new CFLOBDDTopNodeFourier(reduced_tempHandle, m));
		}

		CFLOBDDTopNodeFourierRefPtr MatrixShiftToAConnectionTop(CFLOBDDTopNodeFourierRefPtr c)
		{
			CFLOBDDTopNodeFourierRefPtr v;
			CFLOBDDNodeHandle tempHandle;

			tempHandle = MatrixShiftToAConnectionNode(*(c->rootConnection.entryPointHandle));
			v = new CFLOBDDTopNodeFourier(tempHandle, c->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeFourierRefPtr MatrixShiftToBConnectionTop(CFLOBDDTopNodeFourierRefPtr c)
		{
			CFLOBDDTopNodeFourierRefPtr v;
			CFLOBDDNodeHandle tempHandle;

			tempHandle = MatrixShiftToBConnectionNode(*(c->rootConnection.entryPointHandle));
			v = new CFLOBDDTopNodeFourier(tempHandle, c->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeFourierRefPtr MatrixMultiplyV4TopNode(CFLOBDDTopNodeFourierRefPtr c1, CFLOBDDTopNodeFourierRefPtr c2)
		{
			std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash> hashMap;
			//if (c1->level >= 5)
				//clearMultMap();
			CFLOBDDTopNodeMatMultMapRefPtr c = MatrixMultiplyV4Node(hashMap, *(c1->rootConnection.entryPointHandle), *(c2->rootConnection.entryPointHandle));
			FourierReturnMapHandle v;
			std::unordered_map<fourierSemiring, unsigned int, fourierSemiring::fourierSemiring_hash> reductionMap;
			ReductionMapHandle reductionMapHandle;
			for (unsigned int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++){
				MatMultMapHandle r = c->rootConnection.returnMapHandle[i];
				fourierSemiring val;
				for (auto &j : r.mapContents->map){
					unsigned int index1 = j.first.first;
					unsigned int index2 = j.first.second;
					auto factor = j.second.convert_to<unsigned long long int>();
					val = val + (factor * (c1->rootConnection.returnMapHandle[index1] * c2->rootConnection.returnMapHandle[index2]));
					//std::cout << "index1: " << index1 << " index2: " << index2 << " factor: " << factor << " val: " << val << std::endl;
				}
				if (reductionMap.find(val) == reductionMap.end()){
					v.AddToEnd(val);
					reductionMap.insert(std::make_pair(val, v.Size() - 1));
					reductionMapHandle.AddToEnd(v.Size() - 1);
				}
				else{
					reductionMapHandle.AddToEnd(reductionMap[val]);
				}
			}

			v.Canonicalize();
			reductionMapHandle.Canonicalize();
			//std::cout << v << std::endl;
			CFLOBDDNodeHandle tempHandle = *(c->rootConnection.entryPointHandle);
			// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
			/*ReductionMapHandle inducedReductionMapHandle;
			FloatBoostReturnMapHandle inducedReturnMap;
			v.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());*/
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, v.Size(), true);
			// Create and return CFLOBDDTopNode
			//return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, inducedReturnMap));
			return(new CFLOBDDTopNodeFourier(reduced_tempHandle, v));
		}

		CFLOBDDTopNodeFourierRefPtr MatrixMultiplyV4WithInfoTopNode(CFLOBDDTopNodeFourierRefPtr c1, CFLOBDDTopNodeFourierRefPtr c2)
		{
			//if (c1->level >= 5)
				//clearMultMap();
			std::unordered_map<ZeroValNodeInfo, ZeroIndicesMapHandle, ZeroValNodeInfo::ZeroValNodeInfoHash> hashMap;
			int c1_zero_index = -1, c2_zero_index = -1;
			fourierSemiring zero;
			c1_zero_index = c1->rootConnection.returnMapHandle.LookupInv(zero);
			c2_zero_index = c2->rootConnection.returnMapHandle.LookupInv(zero);
			CFLOBDDTopNodeMatMultMapRefPtr c = MatrixMultiplyV4WithInfoNode
				(hashMap, *(c1->rootConnection.entryPointHandle), *(c2->rootConnection.entryPointHandle),
				c1_zero_index, c2_zero_index);
			FourierReturnMapHandle v;
			std::unordered_map<fourierSemiring, unsigned int, fourierSemiring::fourierSemiring_hash> reductionMap;
			ReductionMapHandle reductionMapHandle;
			for (unsigned int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++){
				MatMultMapHandle r = c->rootConnection.returnMapHandle[i];
				fourierSemiring val;
				for (auto &j : r.mapContents->map){
					unsigned int index1 = j.first.first;
					unsigned int index2 = j.first.second;
					if (index1 != -1 && index2 != -1){
						auto factor = j.second.convert_to<unsigned long long int>();
						val = val + (factor * (c1->rootConnection.returnMapHandle[index1] * c2->rootConnection.returnMapHandle[index2]));
					}
				}
				if (reductionMap.find(val) == reductionMap.end()){
					v.AddToEnd(val);
					reductionMap.insert(std::make_pair(val, v.Size() - 1));
					reductionMapHandle.AddToEnd(v.Size() - 1);
				}
				else{
					reductionMapHandle.AddToEnd(reductionMap[val]);
				}
			}

			v.Canonicalize();
			reductionMapHandle.Canonicalize();
			//std::cout << v << std::endl;
			CFLOBDDNodeHandle tempHandle = *(c->rootConnection.entryPointHandle);
			// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
			/*ReductionMapHandle inducedReductionMapHandle;
			FloatBoostReturnMapHandle inducedReturnMap;
			v.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());*/
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, v.Size(), true);
			// Create and return CFLOBDDTopNode
			//return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, inducedReturnMap));
			return(new CFLOBDDTopNodeFourier(reduced_tempHandle, v));
		}
	}
}

